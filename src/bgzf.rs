//! Opening plain, gzip and BGZF files, decompressing BGZF blocks in parallel.

use std::cell::RefCell;
use std::fs::File;
use std::io::{self, BufRead, BufReader, Read, Seek, SeekFrom};
use std::path::Path;

use anyhow::{Context, Result};
use flate2::read::MultiGzDecoder;
use rayon::prelude::*;
use zflate::{CompressionLevel, Compressor, Decompressor};

/// BGZF blocks decompressed together; 256 blocks are at most 16 MiB of output.
const BLOCKS_PER_BATCH: usize = 256;

/// Fixed part of a gzip member header, up to and including XLEN.
const GZIP_HEADER_LEN: usize = 12;

thread_local! {
    static DECOMPRESSOR: RefCell<Decompressor> = RefCell::new(Decompressor::new());
}

/// Opens `path` for reading, transparently decompressing gzip and BGZF.
///
/// # Errors
///
/// Fails when the file cannot be opened or its first bytes cannot be read.
pub fn open(path: &Path) -> Result<Box<dyn Read + Send>> {
    let mut file = File::open(path).with_context(|| format!("opening {}", path.display()))?;
    let mut head = [0u8; 16];
    let n = read_up_to(&mut file, &mut head)?;
    file.seek(SeekFrom::Start(0))?;
    let head = &head[..n];
    if is_bgzf(head) {
        Ok(Box::new(BgzfReader::new(file)))
    } else if head.starts_with(&[0x1f, 0x8b]) {
        let buffered = BufReader::with_capacity(1 << 20, file);
        Ok(Box::new(BufReader::with_capacity(
            1 << 20,
            MultiGzDecoder::new(buffered),
        )))
    } else {
        Ok(Box::new(BufReader::with_capacity(1 << 20, file)))
    }
}

/// Opens standard input, transparently decompressing gzip and BGZF.
///
/// # Errors
///
/// Fails when standard input cannot be read.
pub fn open_stdin() -> Result<Box<dyn Read + Send>> {
    let mut stdin = BufReader::with_capacity(1 << 20, io::stdin());
    let head = stdin.fill_buf().context("reading standard input")?;
    if is_bgzf(head) {
        Ok(Box::new(BgzfReader::new(stdin)))
    } else if head.starts_with(&[0x1f, 0x8b]) {
        Ok(Box::new(BufReader::with_capacity(
            1 << 20,
            MultiGzDecoder::new(stdin),
        )))
    } else {
        Ok(Box::new(stdin))
    }
}

/// Uncompressed bytes per BGZF block written, below the format's 64 KiB limit.
const WRITE_BLOCK: usize = 0xff00;

/// The 28-byte empty block that ends a BGZF file.
pub const EOF_BLOCK: [u8; 28] = [
    0x1f, 0x8b, 8, 4, 0, 0, 0, 0, 0, 0xff, 6, 0, b'B', b'C', 2, 0, 0x1b, 0, 3, 0, 0, 0, 0, 0, 0, 0,
    0, 0,
];

thread_local! {
    static COMPRESSOR: RefCell<Compressor> =
        RefCell::new(Compressor::new(CompressionLevel::new(5).expect("level 5 exists")));
}

/// Appends `data` to `out` as BGZF blocks (without the end-of-file block).
pub fn compress(data: &[u8], out: &mut Vec<u8>) {
    let mut buffer = vec![0u8; zflate::deflate_compress_bound(WRITE_BLOCK)];
    for chunk in data.chunks(WRITE_BLOCK) {
        let len = COMPRESSOR
            .with(|c| c.borrow_mut().deflate_compress(chunk, &mut buffer))
            .expect("a buffer of the compress bound always fits");
        let total = 18 + len + 8;
        out.extend_from_slice(&[
            0x1f, 0x8b, 8, 4, 0, 0, 0, 0, 0, 0xff, 6, 0, b'B', b'C', 2, 0,
        ]);
        out.extend_from_slice(
            &u16::try_from(total - 1)
                .expect("a block fits in 64 KiB")
                .to_le_bytes(),
        );
        out.extend_from_slice(&buffer[..len]);
        out.extend_from_slice(&zflate::crc32(0, chunk).to_le_bytes());
        out.extend_from_slice(
            &u32::try_from(chunk.len())
                .expect("block is small")
                .to_le_bytes(),
        );
    }
}

fn read_up_to(reader: &mut impl Read, buf: &mut [u8]) -> io::Result<usize> {
    let mut filled = 0;
    while filled < buf.len() {
        match reader.read(&mut buf[filled..])? {
            0 => break,
            n => filled += n,
        }
    }
    Ok(filled)
}

/// A gzip member with the FEXTRA flag whose first extra subfield is BGZF's `BC`.
fn is_bgzf(head: &[u8]) -> bool {
    head.len() >= 16
        && head[..4] == [0x1f, 0x8b, 8, 4]
        && u16::from_le_bytes([head[10], head[11]]) >= 6
        && head[12..14] == *b"BC"
}

/// Reads a BGZF stream, inflating batches of blocks on the rayon pool.
///
/// A plain gzip member after BGZF ones (files concatenated with `cat`) is legal gzip: from
/// there on the rest of the stream is read as ordinary multi-member gzip.
pub struct BgzfReader<R> {
    /// `None` once the stream has been handed to `plain`.
    inner: Option<R>,
    out: Vec<u8>,
    pos: usize,
    done: bool,
    /// Whether the last block read was empty, as the end-of-file block is.
    last_empty: bool,
    /// The bytes of a plain gzip member's start, read before it was recognized.
    pending_plain: Option<Vec<u8>>,
    plain: Option<MultiGzDecoder<io::Chain<io::Cursor<Vec<u8>>, R>>>,
}

impl<R> std::fmt::Debug for BgzfReader<R> {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.debug_struct("BgzfReader")
            .field("pos", &self.pos)
            .field("done", &self.done)
            .field("plain", &self.plain.is_some())
            .finish_non_exhaustive()
    }
}

/// What the stream holds next.
enum Next {
    Block(Vec<u8>),
    End,
    /// A gzip member that is not BGZF, with the bytes of it already read.
    Plain(Vec<u8>),
}

impl<R: Read> BgzfReader<R> {
    pub fn new(inner: R) -> Self {
        Self {
            inner: Some(inner),
            out: Vec::new(),
            pos: 0,
            done: false,
            last_empty: false,
            pending_plain: None,
            plain: None,
        }
    }

    fn inner(&mut self) -> &mut R {
        self.inner
            .as_mut()
            .expect("read only before switching to plain gzip")
    }

    /// Reads one block's bytes after the fixed header.
    fn next_raw_block(&mut self) -> io::Result<Next> {
        let mut header = [0u8; GZIP_HEADER_LEN];
        let last_empty = self.last_empty;
        match read_up_to(self.inner(), &mut header)? {
            // Writers end BGZF with an empty block. A file cut at a block boundary (a
            // killed basecaller, an interrupted copy) lacks it and would otherwise read as
            // complete, with some of its reads missing.
            0 if last_empty => return Ok(Next::End),
            0 => {
                return Err(io::Error::new(
                    io::ErrorKind::UnexpectedEof,
                    "BGZF file ends without its end-of-file block: it is truncated (an interrupted write or copy)",
                ))
            }
            GZIP_HEADER_LEN => {}
            n => return Err(invalid(format!("truncated BGZF header ({n} bytes)"))),
        }
        if header[..3] != [0x1f, 0x8b, 8] {
            return Err(invalid("not a BGZF block".to_owned()));
        }
        if header[3] & 4 == 0 {
            return Ok(Next::Plain(header.to_vec()));
        }
        let xlen = usize::from(u16::from_le_bytes([header[10], header[11]]));
        let mut extra = vec![0u8; xlen];
        self.inner().read_exact(&mut extra)?;
        let Some(block_size) = bsize(&extra) else {
            // A gzip member with extra fields but no BGZF size: plain gzip too.
            let mut seen = header.to_vec();
            seen.extend_from_slice(&extra);
            return Ok(Next::Plain(seen));
        };
        // The block holds the header, the extra field, the data and an 8-byte trailer.
        let rest = (block_size + 1)
            .checked_sub(GZIP_HEADER_LEN + xlen)
            .filter(|&rest| rest >= 8)
            .ok_or_else(|| invalid(format!("invalid BGZF block size {}", block_size + 1)))?;
        let mut raw = vec![0u8; rest];
        self.inner().read_exact(&mut raw)?;
        self.last_empty = raw[rest - 4..] == [0, 0, 0, 0];
        Ok(Next::Block(raw))
    }

    fn fill(&mut self) -> io::Result<()> {
        let mut raw_blocks = Vec::with_capacity(BLOCKS_PER_BATCH);
        while raw_blocks.len() < BLOCKS_PER_BATCH {
            match self.next_raw_block()? {
                Next::Block(raw) => raw_blocks.push(raw),
                Next::End => {
                    self.done = true;
                    break;
                }
                Next::Plain(seen) => {
                    self.pending_plain = Some(seen);
                    break;
                }
            }
        }
        let inflated: Vec<Vec<u8>> = raw_blocks
            .par_iter()
            .map(|raw| inflate(raw))
            .collect::<io::Result<_>>()?;
        self.out.clear();
        self.pos = 0;
        for block in inflated {
            self.out.extend_from_slice(&block);
        }
        Ok(())
    }
}

impl<R: Read> Read for BgzfReader<R> {
    fn read(&mut self, buf: &mut [u8]) -> io::Result<usize> {
        while self.pos == self.out.len() {
            if let Some(plain) = &mut self.plain {
                return plain.read(buf);
            }
            if let Some(seen) = self.pending_plain.take() {
                let inner = self.inner.take().expect("not yet switched");
                self.plain = Some(MultiGzDecoder::new(io::Cursor::new(seen).chain(inner)));
                continue;
            }
            if self.done {
                return Ok(0);
            }
            self.fill()?;
        }
        let n = buf.len().min(self.out.len() - self.pos);
        buf[..n].copy_from_slice(&self.out[self.pos..self.pos + n]);
        self.pos += n;
        Ok(n)
    }
}

/// BSIZE from the extra subfields: the total block size minus one.
fn bsize(mut extra: &[u8]) -> Option<usize> {
    while extra.len() >= 4 {
        let len = usize::from(u16::from_le_bytes([extra[2], extra[3]]));
        let data = extra.get(4..4 + len)?;
        if extra[..2] == *b"BC" && len == 2 {
            return Some(usize::from(u16::from_le_bytes([data[0], data[1]])));
        }
        extra = &extra[4 + len..];
    }
    None
}

/// Inflates a block's compressed data plus CRC32 and ISIZE trailer.
fn inflate(raw: &[u8]) -> io::Result<Vec<u8>> {
    let Some(split) = raw.len().checked_sub(8) else {
        return Err(invalid("BGZF block too short".to_owned()));
    };
    let (data, trailer) = raw.split_at(split);
    let crc = u32::from_le_bytes([trailer[0], trailer[1], trailer[2], trailer[3]]);
    let size = u32::from_le_bytes([trailer[4], trailer[5], trailer[6], trailer[7]]) as usize;
    let mut out = vec![0u8; size];
    if size > 0 {
        let written = DECOMPRESSOR
            .with(|d| d.borrow_mut().deflate_decompress(data, &mut out))
            .map_err(|e| invalid(format!("inflating BGZF block: {e}")))?;
        if written != size {
            return Err(invalid(format!(
                "BGZF block inflated to {written} bytes, expected {size}"
            )));
        }
    }
    if zflate::crc32(0, &out) != crc {
        return Err(invalid("BGZF block CRC mismatch".to_owned()));
    }
    Ok(out)
}

fn invalid(message: String) -> io::Error {
    io::Error::new(io::ErrorKind::InvalidData, message)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn refuses_streams_without_the_eof_block() {
        let mut out = Vec::new();
        compress(&[7u8; 200_000], &mut out);
        let mut back = Vec::new();
        let err = BgzfReader::new(&out[..])
            .read_to_end(&mut back)
            .expect_err("cut at a block boundary");
        assert!(err.to_string().contains("end-of-file block"), "{err}");
        out.extend_from_slice(&EOF_BLOCK);
        back.clear();
        BgzfReader::new(&out[..])
            .read_to_end(&mut back)
            .expect("complete stream");
        assert_eq!(back.len(), 200_000);
    }

    #[test]
    fn reads_plain_gzip_members_after_bgzf() {
        let mut data = bgzf(b"@a\nACGT\n+\nIIII\n", 64);
        let mut plain = Vec::new();
        {
            let mut gz = flate2::write::GzEncoder::new(&mut plain, flate2::Compression::default());
            std::io::Write::write_all(&mut gz, b"@b\nGG\n+\nII\n").expect("compress");
            gz.finish().expect("finish");
        }
        data.extend_from_slice(&plain);
        let mut back = String::new();
        BgzfReader::new(&data[..])
            .read_to_string(&mut back)
            .expect("BGZF then gzip");
        assert_eq!(back, "@a\nACGT\n+\nIIII\n@b\nGG\n+\nII\n");
    }

    #[test]
    fn rejects_impossible_block_sizes() {
        let mut out = bgzf(b"some data", 10);
        out[16] = 5; // BSIZE smaller than the block's header
        out[17] = 0;
        let err = BgzfReader::new(&out[..])
            .read_to_end(&mut Vec::new())
            .expect_err("corrupt block");
        assert!(err.to_string().contains("block size"), "{err}");
    }
    use std::io::Write;

    /// Writes `data` as BGZF blocks of at most `chunk` bytes plus the EOF block.
    pub(crate) fn bgzf(data: &[u8], chunk: usize) -> Vec<u8> {
        let mut out = Vec::new();
        let mut write_block = |payload: &[u8]| {
            let mut encoder =
                flate2::write::DeflateEncoder::new(Vec::new(), flate2::Compression::fast());
            encoder.write_all(payload).expect("in-memory write");
            let cdata = encoder.finish().expect("in-memory write");
            let total = 18 + cdata.len() + 8;
            out.extend_from_slice(&[
                0x1f, 0x8b, 8, 4, 0, 0, 0, 0, 0, 0xff, 6, 0, b'B', b'C', 2, 0,
            ]);
            out.extend_from_slice(&u16::try_from(total - 1).expect("block fits").to_le_bytes());
            out.extend_from_slice(&cdata);
            out.extend_from_slice(&zflate::crc32(0, payload).to_le_bytes());
            out.extend_from_slice(
                &u32::try_from(payload.len())
                    .expect("block fits")
                    .to_le_bytes(),
            );
        };
        for piece in data.chunks(chunk) {
            write_block(piece);
        }
        write_block(&[]);
        out
    }

    #[test]
    fn round_trips_many_blocks() {
        let data: Vec<u8> = (0..1_000_000u32).map(|i| (i % 251) as u8).collect();
        let compressed = bgzf(&data, 60_000);
        let mut reader = BgzfReader::new(&compressed[..]);
        let mut back = Vec::new();
        reader.read_to_end(&mut back).expect("valid stream");
        assert_eq!(back, data);
    }

    #[test]
    fn compresses_readable_bgzf() {
        let data: Vec<u8> = (0..300_000u32).flat_map(|i| i.to_le_bytes()).collect();
        let mut out = Vec::new();
        compress(&data, &mut out);
        out.extend_from_slice(&EOF_BLOCK);
        assert!(is_bgzf(&out[..16]));
        let mut back = Vec::new();
        BgzfReader::new(&out[..])
            .read_to_end(&mut back)
            .expect("valid BGZF");
        assert_eq!(back, data);
        let mut gz = Vec::new();
        MultiGzDecoder::new(&out[..])
            .read_to_end(&mut gz)
            .expect("valid gzip");
        assert_eq!(gz, data);
    }

    #[test]
    fn detects_bgzf_header() {
        assert!(is_bgzf(&bgzf(b"x", 10)[..16]));
        assert!(!is_bgzf(b"@read\nACGT\n+\n!!!!"));
    }
}
