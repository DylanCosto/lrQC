//! FASTQ and FASTA record splitting and decoding, including the metadata basecallers put
//! in read headers.

use std::io::BufRead;

use anyhow::{ensure, Result};

use crate::barcode::LabelCache;
use crate::read::{DuplexStatus, Molecule, Origin, Platform, Read, ASCII_TO_CODE};
use crate::source::Source;
use crate::time;

/// Appends the next FASTQ record to `buf`: usually four lines, but sequence and qualities
/// may each be wrapped over several lines, as the format allows. Returns `false` at end of
/// input.
///
/// # Errors
///
/// Fails on I/O errors and on a record cut short.
pub fn read_fastq_record(reader: &mut impl BufRead, buf: &mut Vec<u8>) -> Result<bool> {
    let start = buf.len();
    // Skip blank lines between records.
    loop {
        let n = reader.read_until(b'\n', buf)?;
        if n == 0 {
            return Ok(false);
        }
        if buf[start..].iter().any(|b| !b.is_ascii_whitespace()) {
            break;
        }
        buf.truncate(start);
    }
    ensure!(buf[start] == b'@', "FASTQ record does not start with '@'");
    let content = |line: &[u8]| {
        let line = line.strip_suffix(b"\n").unwrap_or(line);
        line.strip_suffix(b"\r").unwrap_or(line).len()
    };
    // Sequence lines, up to the '+' line. The first line is sequence even when empty.
    let mut seq_len = 0;
    let mut first = true;
    loop {
        let line_start = buf.len();
        ensure!(reader.read_until(b'\n', buf)? > 0, "truncated FASTQ record");
        if !first && buf[line_start] == b'+' {
            break;
        }
        first = false;
        seq_len += content(&buf[line_start..]);
    }
    // Quality lines, until they cover the sequence. A quality line may start with '@'.
    let mut qual_len = 0;
    loop {
        let line_start = buf.len();
        ensure!(reader.read_until(b'\n', buf)? > 0, "truncated FASTQ record");
        qual_len += content(&buf[line_start..]);
        if qual_len >= seq_len {
            break;
        }
    }
    Ok(true)
}

/// Streams FASTA records, which span a header line and any number of sequence lines.
#[derive(Debug, Default)]
pub struct FastaSplitter {
    pending: Vec<u8>,
}

impl FastaSplitter {
    /// Appends the next FASTA record to `buf`. Returns `false` at end of input.
    ///
    /// # Errors
    ///
    /// Fails on I/O errors and on input that does not start with '>'.
    pub fn next_record(&mut self, reader: &mut impl BufRead, buf: &mut Vec<u8>) -> Result<bool> {
        let mut line = std::mem::take(&mut self.pending);
        // Blank lines before the first record, or at the end, hold no record.
        while line.iter().all(u8::is_ascii_whitespace) {
            line.clear();
            if reader.read_until(b'\n', &mut line)? == 0 {
                return Ok(false);
            }
        }
        ensure!(
            line.first() == Some(&b'>'),
            "FASTA record does not start with '>'"
        );
        buf.extend_from_slice(&line);
        loop {
            line.clear();
            if reader.read_until(b'\n', &mut line)? == 0 {
                return Ok(true);
            }
            if line.first() == Some(&b'>') {
                self.pending = line;
                return Ok(true);
            }
            buf.extend_from_slice(&line);
        }
    }
}

/// Decodes a FASTQ (four lines) or FASTA (header plus sequence lines) record.
///
/// # Errors
///
/// Fails when a FASTQ record's quality and sequence lengths differ.
pub fn decode(
    record: &[u8],
    source: &Source,
    cache: &mut LabelCache,
    read: &mut Read,
) -> Result<()> {
    read.clear();
    read.platform = source.platform;
    read.barcode = source.barcode.clone();
    read.filter = source.filter;
    read.run = source.run.clone();

    let mut lines = record
        .split(|&b| b == b'\n')
        .map(|l| l.strip_suffix(b"\r").unwrap_or(l));
    let header = lines.next().unwrap_or_default();
    let header = &header[1.min(header.len())..];
    parse_header(header, cache, read);

    if record.first() == Some(&b'@') {
        // Sequence lines until the '+' line (usually one), then quality lines.
        let mut first = true;
        for line in lines.by_ref() {
            if !first && line.first() == Some(&b'+') {
                break;
            }
            first = false;
            read.bases
                .extend(line.iter().map(|&b| ASCII_TO_CODE[usize::from(b)]));
        }
        for line in lines {
            read.quals
                .extend(line.iter().map(|&q| q.saturating_sub(33)));
        }
        ensure!(
            read.quals.len() == read.bases.len(),
            "FASTQ record has {} bases but {} qualities",
            read.bases.len(),
            read.quals.len()
        );
    } else {
        for line in lines {
            read.bases
                .extend(line.iter().map(|&b| ASCII_TO_CODE[usize::from(b)]));
        }
    }
    read.len = read.bases.len() as u64;
    read.sanitize();
    if read.platform == Platform::Unknown {
        read.platform = guess_platform(header, read);
    }
    if read.platform == Platform::PacBio {
        read.reported_q = None;
    }
    Ok(())
}

/// Reads the name and the metadata after it: dorado's SAM-style tags (`qs:f:18.2`,
/// `BC:Z:barcode01`) and guppy/MinKNOW's `key=value` pairs (`barcode=barcode01`, `ch=12`).
fn parse_header(header: &[u8], cache: &mut LabelCache, read: &mut Read) {
    let mut fields = header
        .split(|b| b.is_ascii_whitespace())
        .filter(|f| !f.is_empty());
    let name = fields.next().unwrap_or_default();
    read.name.extend_from_slice(name);
    read.name_hash = crate::read::hash_name(name);
    read.zmw = crate::read::zmw_from_name(name);
    if let Some(movie) = crate::read::movie_from_name(name) {
        read.run = Some(cache.intern_verbatim(movie));
    }
    let mut tag_barcode = None;
    let mut header_barcode = None;
    let mut group_barcode = None;
    for field in fields {
        if field.len() > 5 && field[2] == b':' && field[4] == b':' {
            let value = &field[5..];
            match &field[..2] {
                b"BC" if !crate::read::is_barcode_sequence(value) => {
                    tag_barcode = Some(cache.intern(value));
                }
                b"qs" => read.reported_q = parse(value),
                b"ch" => read.channel = parse(value),
                b"st" => read.start_secs = time::parse_iso8601(value),
                b"du" => read.duration_secs = parse(value),
                b"rq" => read.accuracy = parse(value),
                b"np" => read.passes = parse(value),
                b"er" => read.end_reason = Some(cache.intern_verbatim(value)),
                b"dx" => read.duplex = parse(value).and_then(DuplexStatus::from_tag),
                b"pi" => read.origin = Origin::SplitFromParent,
                b"mx" => read.mux = parse(value),
                b"sm" => read.signal_shift = parse(value),
                b"sd" => read.signal_scale = parse(value),
                b"ns" => read.signal_samples = parse(value),
                b"pt" => read.poly_a = parse(value),
                b"zm" => read.zmw = parse(value),
                b"RG" => {
                    group_barcode = cache.read_group_barcode(value);
                    read.molecule = Molecule::from_model_name(&String::from_utf8_lossy(value));
                    // dorado's read groups start with the run ID.
                    let run = value.split(|&b| b == b'_').next().unwrap_or(value);
                    read.run = Some(cache.intern_verbatim(&run[..run.len().min(8)]));
                }
                b"MM" | b"Mm" => read.mod_positions.extend_from_slice(value),
                b"ML" | b"Ml" => {
                    // "C,12,255,..." in SAM text form.
                    let values = value.split(|&b| b == b',').skip(1);
                    read.mod_probs.extend(values.filter_map(parse::<u8>));
                }
                _ => {}
            }
        } else if let Some((key, value)) = split_once(field, b'=') {
            match key {
                b"barcode" => header_barcode = Some(cache.intern(value)),
                b"ch" => read.channel = parse(value),
                b"start_time" => read.start_secs = time::parse_iso8601(value),
                b"np" => read.passes = parse(value),
                b"rq" => read.accuracy = parse(value),
                b"mux" => read.mux = parse(value),
                b"runid" => read.run = Some(cache.intern_verbatim(&value[..value.len().min(8)])),
                // MinKNOW names the model (`rna004_130bps_sup@v5.1.0`), dorado the read group.
                b"basecall_model_version_id" | b"basecall_model" | b"model_version_id" => {
                    read.molecule = Molecule::from_model_name(&String::from_utf8_lossy(value));
                }
                _ => {}
            }
        }
    }
    // A BC tag or barcode= names the barcode directly; dorado FASTQ otherwise names it
    // only at the end of the read group.
    if let Some(label) = tag_barcode.or(header_barcode).or(group_barcode) {
        read.barcode = Some(label);
    }
}

fn split_once(bytes: &[u8], delimiter: u8) -> Option<(&[u8], &[u8])> {
    let at = memchr::memchr(delimiter, bytes)?;
    Some((&bytes[..at], &bytes[at + 1..]))
}

fn parse<T: std::str::FromStr>(bytes: &[u8]) -> Option<T> {
    std::str::from_utf8(bytes).ok()?.parse().ok()
}

/// Platform from the read name and header metadata.
fn guess_platform(header: &[u8], read: &Read) -> Platform {
    let name = header
        .split(|b| b.is_ascii_whitespace())
        .next()
        .unwrap_or_default();
    let pacbio_name = name.starts_with(b"m") && name.iter().filter(|&&b| b == b'/').count() >= 2;
    if pacbio_name || read.accuracy.is_some() || read.passes.is_some() {
        Platform::PacBio
    } else if read.channel.is_some()
        || read.start_secs.is_some()
        || read.reported_q.is_some()
        || is_uuid(name)
    {
        Platform::Ont
    } else {
        Platform::Unknown
    }
}

/// ONT read IDs are UUIDs: 8-4-4-4-12 hex digits.
fn is_uuid(name: &[u8]) -> bool {
    name.len() == 36
        && name.iter().enumerate().all(|(i, &b)| match i {
            8 | 13 | 18 | 23 => b == b'-',
            _ => b.is_ascii_hexdigit(),
        })
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::io::Cursor;

    #[test]
    fn splits_and_decodes_fastq() {
        let input = b"@r1 qs:f:18.5\tBC:Z:SQK-RBK114-24_barcode03 ch:i:7\nACGTN\n+\n!+5?I\n\n@r2\nAC\n+\nII\n";
        let mut reader = Cursor::new(&input[..]);
        let mut buf = Vec::new();
        assert!(read_fastq_record(&mut reader, &mut buf).expect("valid"));
        let mut read = Read::default();
        decode(
            &buf,
            &Source::default(),
            &mut LabelCache::default(),
            &mut read,
        )
        .expect("valid");
        assert_eq!(read.bases, vec![0, 1, 2, 3, 4]);
        assert_eq!(read.quals, vec![0, 10, 20, 30, 40]);
        assert_eq!(read.barcode.as_deref(), Some("barcode03"));
        assert_eq!(read.reported_q, Some(18.5));
        assert_eq!(read.channel, Some(7));
        assert_eq!(read.platform, Platform::Ont);
        buf.clear();
        assert!(read_fastq_record(&mut reader, &mut buf).expect("valid"));
        buf.clear();
        assert!(!read_fastq_record(&mut reader, &mut buf).expect("valid"));
    }

    #[test]
    fn reads_barcodes_from_dorado_read_groups() {
        let record = b"@75d661a6-68e6-4770-9200-df57dff9d3cb\tqs:f:19.56\tRG:Z:7a4eec12-4b56_dna_r10.4.1_e8.2_400bps_sup@v5.2.0_barcode03\nAC\n+\nII\n";
        let mut read = Read::default();
        decode(
            record,
            &Source::default(),
            &mut LabelCache::default(),
            &mut read,
        )
        .expect("valid");
        assert_eq!(read.barcode.as_deref(), Some("barcode03"));
        // A BC tag wins over the read group.
        let tagged =
            b"@r\tBC:Z:SQK-RBK114-24_barcode07\tRG:Z:x_dna_sup@v5.2.0_barcode03\nAC\n+\nII\n";
        decode(
            tagged,
            &Source::default(),
            &mut LabelCache::default(),
            &mut read,
        )
        .expect("valid");
        assert_eq!(read.barcode.as_deref(), Some("barcode07"));
    }

    #[test]
    fn reads_guppy_headers() {
        let record = b"@feba0342-9b5e-408e-a8f0-395a471477f6 runid=abc ch=12 start_time=2024-01-01T00:00:00Z barcode=barcode11\nAC\n+\nII\n";
        let mut read = Read::default();
        decode(
            record,
            &Source::default(),
            &mut LabelCache::default(),
            &mut read,
        )
        .expect("valid");
        assert_eq!(read.barcode.as_deref(), Some("barcode11"));
        assert_eq!(read.channel, Some(12));
        assert_eq!(read.platform, Platform::Ont);
    }

    #[test]
    fn splits_multiline_fasta() {
        let input = b">m84011/1/ccs\nACG\nTT\n>m84011/2/ccs\nGG\n";
        let mut reader = Cursor::new(&input[..]);
        let mut splitter = FastaSplitter::default();
        let mut buf = Vec::new();
        assert!(splitter.next_record(&mut reader, &mut buf).expect("valid"));
        let mut read = Read::default();
        decode(
            &buf,
            &Source::default(),
            &mut LabelCache::default(),
            &mut read,
        )
        .expect("valid");
        assert_eq!(read.len, 5);
        assert_eq!(read.platform, Platform::PacBio);
        buf.clear();
        assert!(splitter.next_record(&mut reader, &mut buf).expect("valid"));
        assert_eq!(buf, b">m84011/2/ccs\nGG\n");
        buf.clear();
        assert!(!splitter.next_record(&mut reader, &mut buf).expect("valid"));
    }
}
