//! BAM header and record decoding, into [`Read`].

use std::collections::HashMap;
use std::io::{self, Read as IoRead};

use std::sync::Arc;

use anyhow::{bail, ensure, Context, Result};

use crate::barcode::{self, LabelCache};
use crate::read::{
    self as reads, Alignment, DuplexStatus, Filter, Kinetics, Mapped, Molecule, Origin, Platform,
    Read, RecordKind, BASE_N,
};
use crate::source::Source;
use crate::time;

/// BAM's 4-bit base encoding, `=ACMGRSVTWYHKDBN`, as base codes.
const NIBBLE_TO_CODE: [u8; 16] = [
    BASE_N, 0, 1, BASE_N, 2, BASE_N, BASE_N, BASE_N, 3, BASE_N, BASE_N, BASE_N, BASE_N, BASE_N,
    BASE_N, BASE_N,
];

const FLAG_UNMAPPED: u16 = 0x4;
const FLAG_REVERSE: u16 = 0x10;
const FLAG_SECONDARY: u16 = 0x100;
const FLAG_QC_FAIL: u16 = 0x200;
const FLAG_SUPPLEMENTARY: u16 = 0x800;

/// Fixed-length part of a record, from `refID` through `tlen`.
const FIXED_LEN: usize = 32;

/// The parts of a BAM header lrqc uses.
#[derive(Debug, Clone, Default)]
pub struct Header {
    pub read_groups: Vec<ReadGroup>,
    /// `@PG` lines, verbatim.
    pub programs: Vec<String>,
    pub references: usize,
}

/// One `@RG` line.
#[derive(Debug, Clone, Default)]
pub struct ReadGroup {
    pub id: String,
    pub fields: Vec<(String, String)>,
    pub platform: Platform,
    /// Barcode the read group stands for (dorado's `..._barcode01` IDs, PacBio's `/0--0`).
    pub barcode: Option<String>,
}

/// Reads the BAM magic, header text and reference list.
///
/// # Errors
///
/// Fails on a stream that is not BAM or is truncated.
pub fn read_header(reader: &mut impl IoRead, barcode_names: &[String]) -> Result<Header> {
    let mut magic = [0u8; 4];
    reader.read_exact(&mut magic).context("reading BAM magic")?;
    ensure!(&magic == b"BAM\x01", "not a BAM file");
    let text_len = read_u32(reader)? as usize;
    let mut text = vec![0u8; text_len];
    reader.read_exact(&mut text)?;
    let references = read_u32(reader)? as usize;
    for _ in 0..references {
        let name_len = read_u32(reader)? as usize;
        let mut skip = vec![0u8; name_len + 4];
        reader.read_exact(&mut skip)?;
    }
    let text = String::from_utf8_lossy(&text);
    let mut header = Header {
        references,
        ..Header::default()
    };
    for line in text.lines() {
        if let Some(rest) = line.strip_prefix("@RG\t") {
            header
                .read_groups
                .push(parse_read_group(rest, barcode_names));
        } else if line.starts_with("@PG\t") {
            header.programs.push(line.to_owned());
        }
    }
    Ok(header)
}

fn parse_read_group(rest: &str, barcode_names: &[String]) -> ReadGroup {
    let fields: Vec<(String, String)> = rest
        .split('\t')
        .filter_map(|f| f.split_once(':'))
        .map(|(k, v)| (k.to_owned(), v.to_owned()))
        .collect();
    let field = |key: &str| {
        fields
            .iter()
            .find(|(k, _)| k == key)
            .map(|(_, v)| v.as_str())
    };
    let id = field("ID").unwrap_or_default().to_owned();
    let platform = match field("PL").map(str::to_ascii_uppercase).as_deref() {
        Some("ONT") | Some("NANOPORE") => Platform::Ont,
        Some("PACBIO") => Platform::PacBio,
        _ => Platform::Unknown,
    };
    let barcode =
        pacbio_group_barcode(&id, barcode_names).or_else(|| barcode::from_read_group(&id));
    ReadGroup {
        id,
        fields,
        platform,
        barcode,
    }
}

/// SMRT Link names demultiplexed read groups `<movie hash>/<forward>--<reverse>`.
fn pacbio_group_barcode(id: &str, barcode_names: &[String]) -> Option<String> {
    let (_, pair) = id.rsplit_once('/')?;
    let (forward, reverse) = pair.split_once("--")?;
    Some(barcode::from_indices(
        forward.parse().ok()?,
        reverse.parse().ok()?,
        barcode_names,
    ))
}

/// Reads one record's bytes (after its `block_size`) into `buf`.
///
/// Returns `false` at a clean end of stream.
///
/// # Errors
///
/// Fails on a truncated record.
pub fn read_record(reader: &mut impl IoRead, buf: &mut Vec<u8>) -> Result<bool> {
    let mut size = [0u8; 4];
    let mut filled = 0;
    while filled < 4 {
        match reader.read(&mut size[filled..])? {
            0 if filled == 0 => return Ok(false),
            0 => bail!("truncated BAM record length"),
            n => filled += n,
        }
    }
    let size = u32::from_le_bytes(size) as usize;
    let start = buf.len();
    buf.resize(start + size, 0);
    reader
        .read_exact(&mut buf[start..])
        .context("truncated BAM record")?;
    Ok(true)
}

fn read_u32(reader: &mut impl IoRead) -> io::Result<u32> {
    let mut bytes = [0u8; 4];
    reader.read_exact(&mut bytes)?;
    Ok(u32::from_le_bytes(bytes))
}

/// Decodes a record into `read`.
///
/// # Errors
///
/// Fails on a record whose lengths point past its end.
pub fn decode(data: &[u8], source: &Source, cache: &mut LabelCache, read: &mut Read) -> Result<()> {
    read.clear();
    ensure!(
        data.len() >= FIXED_LEN,
        "BAM record shorter than its fixed fields"
    );
    let name_len = usize::from(data[8]);
    let mapq = data[9];
    let cigar_ops = usize::from(u16::from_le_bytes([data[12], data[13]]));
    let flag = u16::from_le_bytes([data[14], data[15]]);
    let seq_len = u32::from_le_bytes([data[16], data[17], data[18], data[19]]) as usize;

    let name_end = FIXED_LEN + name_len;
    let cigar_end = name_end + 4 * cigar_ops;
    let seq_end = cigar_end + seq_len.div_ceil(2);
    let qual_end = seq_end + seq_len;
    ensure!(
        data.len() >= qual_end,
        "BAM record shorter than its sequence"
    );
    let name = &data[FIXED_LEN..name_end];
    let name = name.strip_suffix(&[0]).unwrap_or(name);

    read.name.extend_from_slice(name);
    read.name_hash = reads::hash_name(name);
    read.kind = if flag & FLAG_SECONDARY != 0 {
        RecordKind::Secondary
    } else if flag & FLAG_SUPPLEMENTARY != 0 {
        RecordKind::Supplementary
    } else {
        RecordKind::Primary
    };
    read.len = seq_len as u64;
    read.platform = source.platform;
    read.barcode = source.barcode.clone();
    read.filter = source.filter;
    read.run = source.run.clone();
    read.zmw = reads::zmw_from_name(name);
    if let Some(movie) = reads::movie_from_name(name) {
        read.run = Some(cache.intern_verbatim(movie));
    }
    if flag & FLAG_QC_FAIL != 0 {
        read.filter = Some(Filter::Fail);
    }

    let tags = &data[qual_end..];
    let mut cigar = Cigar::from_bytes(&data[name_end..cigar_end]);
    scan_tags(tags, source, cache, read, &mut cigar)?;

    if read.kind != RecordKind::Primary {
        return Ok(());
    }

    let packed = &data[cigar_end..seq_end];
    read.bases.reserve(seq_len);
    for &byte in packed {
        read.bases.push(NIBBLE_TO_CODE[usize::from(byte >> 4)]);
        read.bases.push(NIBBLE_TO_CODE[usize::from(byte & 0xf)]);
    }
    read.bases.truncate(seq_len);
    let quals = &data[seq_end..qual_end];
    if quals.first().is_some_and(|&q| q != 0xff) {
        read.quals.extend_from_slice(quals);
    }

    read.alignment = if flag & FLAG_UNMAPPED != 0 {
        if source.aligned {
            Alignment::Unmapped
        } else {
            Alignment::NotAligned
        }
    } else {
        let mut mapped = cigar.summarize();
        mapped.mapq = mapq;
        mapped.edit_distance = cigar.edit_distance;
        mapped.split = cigar.split;
        if flag & FLAG_REVERSE != 0 {
            read.reverse_complement();
        }
        Alignment::Mapped(mapped)
    };
    Ok(())
}

/// CIGAR operations, from the record or from the `CG` tag of long CIGARs.
#[derive(Debug, Default)]
struct Cigar<'rec> {
    ops: &'rec [u8],
    edit_distance: Option<u32>,
    split: bool,
}

impl<'rec> Cigar<'rec> {
    fn from_bytes(ops: &'rec [u8]) -> Self {
        Self {
            ops,
            ..Self::default()
        }
    }

    fn summarize(&self) -> Mapped {
        let mut mapped = Mapped::default();
        for &op in self.ops.as_chunks::<4>().0 {
            let op = u32::from_le_bytes(op);
            let len = u64::from(op >> 4);
            match op & 0xf {
                0 | 7 | 8 => mapped.aligned_columns += len,
                1 => {
                    mapped.inserted += len;
                    mapped.insertion_runs += 1;
                }
                2 => {
                    mapped.deleted += len;
                    mapped.deletion_runs += 1;
                }
                _ => {}
            }
        }
        mapped
    }
}

/// A tag value, borrowing from the record.
#[derive(Debug)]
enum Value<'rec> {
    Int(i64),
    Float(f32),
    Text(&'rec [u8]),
    Array {
        subtype: u8,
        count: usize,
        bytes: &'rec [u8],
    },
}

impl Value<'_> {
    fn as_int(&self) -> Option<i64> {
        match *self {
            Value::Int(v) => Some(v),
            _ => None,
        }
    }

    fn as_float(&self) -> Option<f32> {
        match *self {
            Value::Float(v) => Some(v),
            // Precision loss is irrelevant at the scale of quality scores.
            Value::Int(v) => Some(v as f32),
            _ => None,
        }
    }
}

fn scan_tags<'rec>(
    mut tags: &'rec [u8],
    source: &Source,
    cache: &mut LabelCache,
    read: &mut Read,
    cigar: &mut Cigar<'rec>,
) -> Result<()> {
    let mut tag_barcode = None;
    let mut group_barcode = None;
    while !tags.is_empty() {
        ensure!(tags.len() >= 3, "truncated BAM tag");
        let tag = [tags[0], tags[1]];
        let (value, rest) = parse_value(tags[2], &tags[3..])?;
        tags = rest;
        match (&tag, value) {
            (b"BC", Value::Text(text)) if !crate::read::is_barcode_sequence(text) => {
                tag_barcode = Some(cache.intern(text));
            }
            (
                b"bc",
                Value::Array {
                    subtype: b'S',
                    count: 2,
                    bytes,
                },
            ) => {
                let forward = u16::from_le_bytes([bytes[0], bytes[1]]);
                let reverse = u16::from_le_bytes([bytes[2], bytes[3]]);
                let name = barcode::from_indices(forward, reverse, &source.barcode_names);
                tag_barcode = Some(cache.intern(name.as_bytes()));
            }
            (b"RG", Value::Text(id)) => {
                if let Some(group) = source.read_groups.get(id) {
                    group_barcode.clone_from(&group.barcode);
                    if group.run.is_some() {
                        read.run.clone_from(&group.run);
                    }
                    read.molecule = group.molecule;
                    if source.platform == Platform::Unknown {
                        read.platform = group.platform;
                    }
                } else {
                    // A read group the header does not list (a BAM rebuilt from FASTQ, a
                    // header stripped by a tool): read its ID as the FASTQ reader does.
                    group_barcode = cache.read_group_barcode(id).map(|b| b.to_string());
                    let run = id.split(|&b| b == b'_').next().unwrap_or(id);
                    if run.len() < id.len() {
                        read.run = Some(cache.intern_verbatim(&run[..run.len().min(8)]));
                    }
                    read.molecule = Molecule::from_model_name(&String::from_utf8_lossy(id));
                }
            }
            (b"bq", v) => read.barcode_score = v.as_int().and_then(|v| u32::try_from(v).ok()),
            (b"qs", Value::Float(q)) => read.reported_q = Some(q),
            (b"rq", v) => read.accuracy = v.as_float(),
            (b"np", v) => read.passes = v.as_int().and_then(|v| u32::try_from(v).ok()),
            (b"ec", v) => read.effective_coverage = v.as_float(),
            (b"ch", v) => read.channel = v.as_int().and_then(|v| u32::try_from(v).ok()),
            (b"st", Value::Text(text)) => read.start_secs = time::parse_iso8601(text),
            (b"du", v) => read.duration_secs = v.as_float(),
            (b"er", Value::Text(reason)) => read.end_reason = Some(cache.intern_verbatim(reason)),
            (b"dx", v) => read.duplex = v.as_int().and_then(DuplexStatus::from_tag),
            (b"pi", Value::Text(_)) => read.origin = Origin::SplitFromParent,
            (b"mx", v) => read.mux = v.as_int().and_then(|v| u8::try_from(v).ok()),
            (b"sm", v) => read.signal_shift = v.as_float(),
            (b"sd", v) => read.signal_scale = v.as_float(),
            (b"ns", v) => read.signal_samples = v.as_int().and_then(|v| u64::try_from(v).ok()),
            (b"pt", v) => read.poly_a = v.as_int().and_then(|v| u32::try_from(v).ok()),
            (b"zm", v) => read.zmw = v.as_int().and_then(|v| u64::try_from(v).ok()),
            (b"fi" | b"ri" | b"ip" | b"pw" | b"fp" | b"rp", Value::Array { .. }) => {
                read.kinetics = Kinetics::Present;
            }
            (b"MM" | b"Mm", Value::Text(text)) => {
                read.mod_positions.clear();
                read.mod_positions.extend_from_slice(text);
            }
            (
                b"ML" | b"Ml",
                Value::Array {
                    subtype: b'C',
                    bytes,
                    ..
                },
            ) => {
                read.mod_probs.clear();
                read.mod_probs.extend_from_slice(bytes);
            }
            (b"NM", v) => cigar.edit_distance = v.as_int().and_then(|v| u32::try_from(v).ok()),
            (b"SA", Value::Text(_)) => cigar.split = true,
            (
                b"CG",
                Value::Array {
                    subtype: b'I',
                    bytes,
                    ..
                },
            ) => cigar.ops = bytes,
            _ => {}
        }
    }
    read.sanitize();
    if read.platform == Platform::Unknown {
        read.platform = guess_platform(read);
    }
    // PacBio's `qs` is a query start, not a quality.
    if read.platform == Platform::PacBio {
        read.reported_q = None;
    }
    if let Some(label) = tag_barcode.or(group_barcode.map(|name| cache.intern(name.as_bytes()))) {
        // A PacBio index pair without --barcode-names names the barcode less well than a
        // lima output file name (`*.bc2001--bc2001.bam`) does.
        let index_only = label.starts_with(barcode::INDEX_PREFIX);
        if !(index_only && source.barcode.is_some()) {
            read.barcode = Some(label);
        }
    }
    Ok(())
}

/// Platform from per-read tags, for files whose header does not say.
fn guess_platform(read: &Read) -> Platform {
    if read.accuracy.is_some() || read.passes.is_some() || read.kinetics == Kinetics::Present {
        Platform::PacBio
    } else if read.channel.is_some() || read.start_secs.is_some() || read.duration_secs.is_some() {
        Platform::Ont
    } else {
        Platform::Unknown
    }
}

/// Parses one tag value of type `kind` from the front of `bytes`, returning the rest.
fn parse_value(kind: u8, bytes: &[u8]) -> Result<(Value<'_>, &[u8])> {
    fn take(bytes: &[u8], n: usize) -> Result<(&[u8], &[u8])> {
        ensure!(bytes.len() >= n, "truncated BAM tag value");
        Ok(bytes.split_at(n))
    }
    Ok(match kind {
        b'A' | b'c' | b'C' => {
            let (v, rest) = take(bytes, 1)?;
            let v = if kind == b'c' {
                i64::from(v[0] as i8)
            } else {
                i64::from(v[0])
            };
            (Value::Int(v), rest)
        }
        b's' | b'S' => {
            let (v, rest) = take(bytes, 2)?;
            let raw = [v[0], v[1]];
            let v = if kind == b's' {
                i64::from(i16::from_le_bytes(raw))
            } else {
                i64::from(u16::from_le_bytes(raw))
            };
            (Value::Int(v), rest)
        }
        b'i' | b'I' => {
            let (v, rest) = take(bytes, 4)?;
            let raw = [v[0], v[1], v[2], v[3]];
            let v = if kind == b'i' {
                i64::from(i32::from_le_bytes(raw))
            } else {
                i64::from(u32::from_le_bytes(raw))
            };
            (Value::Int(v), rest)
        }
        b'f' => {
            let (v, rest) = take(bytes, 4)?;
            (
                Value::Float(f32::from_le_bytes([v[0], v[1], v[2], v[3]])),
                rest,
            )
        }
        b'Z' | b'H' => {
            let end = memchr::memchr(0, bytes).context("unterminated BAM string tag")?;
            (Value::Text(&bytes[..end]), &bytes[end + 1..])
        }
        b'B' => {
            let (head, rest) = take(bytes, 5)?;
            let subtype = head[0];
            let count = u32::from_le_bytes([head[1], head[2], head[3], head[4]]) as usize;
            let width = match subtype {
                b'c' | b'C' => 1,
                b's' | b'S' => 2,
                b'i' | b'I' | b'f' => 4,
                other => bail!("unknown BAM array type {:?}", char::from(other)),
            };
            let (array, rest) = take(rest, count * width)?;
            (
                Value::Array {
                    subtype,
                    count,
                    bytes: array,
                },
                rest,
            )
        }
        other => bail!("unknown BAM tag type {:?}", char::from(other)),
    })
}

/// Read groups by ID, for looking up each record's `RG`.
pub fn read_group_map(header: &Header) -> HashMap<Vec<u8>, crate::source::GroupInfo> {
    header
        .read_groups
        .iter()
        .map(|group| {
            let field = |key: &str| {
                group
                    .fields
                    .iter()
                    .find(|(k, _)| k == key)
                    .map(|(_, v)| v.as_str())
            };
            // dorado's IDs start with the run ID; PacBio's are hashes, with the movie in PU.
            let run = field("PU").map(str::to_owned).or_else(|| {
                let (run, _) = group.id.split_once('_')?;
                Some(run.chars().take(8).collect())
            });
            let info = crate::source::GroupInfo {
                barcode: group.barcode.clone(),
                platform: group.platform,
                run: run.map(|r| Arc::from(r.as_str())),
                molecule: Molecule::from_model_name(&format!(
                    "{} {}",
                    group.id,
                    field("DS").unwrap_or_default()
                )),
            };
            (group.id.as_bytes().to_vec(), info)
        })
        .collect()
}

/// The platform all read groups and programs agree on, if any.
pub fn header_platform(header: &Header) -> Platform {
    let mut platforms = header
        .read_groups
        .iter()
        .map(|g| g.platform)
        .filter(|&p| p != Platform::Unknown);
    if let Some(first) = platforms.next() {
        if platforms.all(|p| p == first) {
            return first;
        }
        return Platform::Unknown;
    }
    let programs = header.programs.join("\n").to_ascii_lowercase();
    if ["id:dorado", "id:basecaller", "id:guppy", "pn:dorado"]
        .iter()
        .any(|p| programs.contains(p))
    {
        Platform::Ont
    } else if ["id:ccs", "id:lima", "id:primrose", "id:jasmine", "pn:ccs"]
        .iter()
        .any(|p| programs.contains(p))
    {
        Platform::PacBio
    } else {
        Platform::Unknown
    }
}

#[cfg(test)]
pub(crate) mod tests {
    use super::*;

    /// Encodes one BAM record (without its `block_size`).
    pub(crate) fn record(
        name: &str,
        flag: u16,
        seq: &[u8],
        qual: &[u8],
        cigar: &[(u32, u32)],
        tags: &[u8],
    ) -> Vec<u8> {
        let mut out = Vec::new();
        let ref_id: i32 = if flag & FLAG_UNMAPPED != 0 { -1 } else { 0 };
        out.extend_from_slice(&ref_id.to_le_bytes());
        out.extend_from_slice(&0i32.to_le_bytes());
        out.push(u8::try_from(name.len() + 1).expect("short name"));
        out.push(60);
        out.extend_from_slice(&0u16.to_le_bytes());
        out.extend_from_slice(
            &u16::try_from(cigar.len())
                .expect("short cigar")
                .to_le_bytes(),
        );
        out.extend_from_slice(&flag.to_le_bytes());
        out.extend_from_slice(&u32::try_from(seq.len()).expect("short seq").to_le_bytes());
        out.extend_from_slice(&(-1i32).to_le_bytes());
        out.extend_from_slice(&(-1i32).to_le_bytes());
        out.extend_from_slice(&0i32.to_le_bytes());
        out.extend_from_slice(name.as_bytes());
        out.push(0);
        for &(len, op) in cigar {
            out.extend_from_slice(&(len << 4 | op).to_le_bytes());
        }
        let code = |b: u8| -> u8 {
            match b {
                b'A' => 1,
                b'C' => 2,
                b'G' => 4,
                b'T' => 8,
                _ => 15,
            }
        };
        for pair in seq.chunks(2) {
            let hi = code(pair[0]);
            let lo = pair.get(1).map_or(0, |&b| code(b));
            out.push(hi << 4 | lo);
        }
        if qual.is_empty() {
            out.extend(std::iter::repeat_n(0xff, seq.len()));
        } else {
            out.extend_from_slice(qual);
        }
        out.extend_from_slice(tags);
        out
    }

    pub(crate) fn text_tag(tag: &[u8; 2], value: &str) -> Vec<u8> {
        let mut out = tag.to_vec();
        out.push(b'Z');
        out.extend_from_slice(value.as_bytes());
        out.push(0);
        out
    }

    #[test]
    fn decodes_reverse_strand_record() {
        let mut tags = text_tag(b"BC", "SQK-NBD114-24_barcode05");
        tags.extend_from_slice(b"NMC\x03");
        let data = record(
            "r1",
            FLAG_REVERSE,
            b"AACGT",
            &[10, 20, 30, 40, 50],
            &[(5, 0)],
            &tags,
        );
        let source = Source::for_tests(Platform::Ont);
        let mut read = Read::default();
        decode(&data, &source, &mut LabelCache::default(), &mut read).expect("valid record");
        assert_eq!(read.barcode.as_deref(), Some("barcode05"));
        // ACGTT reverse-complemented back to sequencing orientation.
        assert_eq!(read.bases, vec![0, 1, 2, 3, 3]);
        assert_eq!(read.quals, vec![50, 40, 30, 20, 10]);
        let Alignment::Mapped(mapped) = read.alignment else {
            panic!("expected a mapped record, got {:?}", read.alignment);
        };
        assert_eq!(mapped.aligned_columns, 5);
        assert_eq!(mapped.edit_distance, Some(3));
    }

    #[test]
    fn names_pacbio_barcode_pairs() {
        let mut tags = b"bcBS\x02\x00\x00\x00\x01\x00\x01\x00".to_vec();
        tags.extend_from_slice(b"rqf");
        tags.extend_from_slice(&0.999f32.to_le_bytes());
        tags.extend_from_slice(b"qsi\x19\x00\x00\x00");
        let data = record("m84/1/ccs", FLAG_UNMAPPED, b"ACGT", &[], &[], &tags);
        let mut source = Source::for_tests(Platform::Unknown);
        source.barcode_names = Arc::new(vec!["bc2001".into(), "bc2002".into()]);
        let mut read = Read::default();
        decode(&data, &source, &mut LabelCache::default(), &mut read).expect("valid record");
        assert_eq!(read.barcode.as_deref(), Some("bc2002--bc2002"));
        assert_eq!(read.platform, Platform::PacBio);
        assert_eq!(read.reported_q, None, "PacBio qs is a query start");
        assert!(read.quals.is_empty());
    }

    #[test]
    fn reads_group_barcodes() {
        let group = parse_read_group(
            "ID:abc_dna_r10.4.1_e8.2_400bps_sup@v5.0.0_SQK-NBD114-24_barcode12\tPL:ONT",
            &[],
        );
        assert_eq!(group.barcode.as_deref(), Some("barcode12"));
        assert_eq!(group.platform, Platform::Ont);
        let group = parse_read_group("ID:e4f3b1d2/3--3\tPL:PACBIO", &[]);
        assert_eq!(group.barcode.as_deref(), Some("index3--index3"));
        assert_eq!(
            parse_read_group("ID:70eb6e24_dna_r10.4.1_e8.2_400bps_sup@v5.0.0", &[]).barcode,
            None
        );
    }
}
