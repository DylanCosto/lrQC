//! Barcode names: normalizing basecaller and demultiplexer labels, and finding them in paths.

use std::sync::{Arc, LazyLock};

use regex::Regex;

/// Name of the group for reads a demultiplexer could not assign.
pub const UNCLASSIFIED: &str = "unclassified";

/// Prefix of PacBio barcode pairs named by index, when the barcode names are unknown.
pub const INDEX_PREFIX: &str = "index";

/// Name of the group for reads with no barcode information at all.
pub const NO_BARCODE: &str = "no_barcode";

/// ONT `barcodeNN`, PacBio `bcNNNN` or `bcNNNN--bcNNNN`, and the unassigned spellings,
/// delimited by anything but a letter or digit.
static BARCODE: LazyLock<Regex> = LazyLock::new(|| {
    Regex::new(
        r"(?i)(?:^|[^a-z0-9])(barcode\d+|bc\d{2,5}(?:--bc\d{2,5})?|unclassified|unassigned|unbarcoded|default--default)(?:$|[^a-z0-9])",
    )
    .expect("barcode pattern is valid")
});

/// Normalizes a barcode label from a record tag or read header.
///
/// `SQK-NBD114-24_barcode01` becomes `barcode01`, the unassigned spellings become
/// [`UNCLASSIFIED`], and labels with no recognizable barcode are kept as they are.
pub fn normalize(label: &str) -> String {
    let label = label.trim();
    match find(label) {
        Some(name) => name,
        None => label.to_owned(),
    }
}

/// The barcode named in `text`, if any.
fn find(text: &str) -> Option<String> {
    let name = BARCODE
        .captures(text)?
        .get(1)?
        .as_str()
        .to_ascii_lowercase();
    Some(match name.as_str() {
        "unassigned" | "unbarcoded" | "default--default" => UNCLASSIFIED.to_owned(),
        _ => name,
    })
}

/// The barcode named by a path: the file name first, then its directories from the
/// innermost out, as in MinKNOW's `fastq_pass/barcode01/` and lima's `*.bc1001--bc1001.bam`.
#[cfg(test)]
fn from_path(path: &std::path::Path) -> Option<String> {
    let parts: Vec<&str> = path.iter().rev().filter_map(|c| c.to_str()).collect();
    from_path_parts(&parts, &[])
}

/// The barcode named by a path split as `input` splits it, innermost first: `inner` (what
/// the user named and everything below it) may name a barcode anywhere in a component
/// (`PAY1_pass_barcode01_abc_0.fastq.gz`); `outer` (directories above it) only by a
/// component that is exactly a barcode (`barcode01/`), so a `proj_bc01/` folder higher up
/// does not label every read.
pub fn from_path_parts(inner: &[&str], outer: &[&str]) -> Option<String> {
    inner.iter().copied().find_map(find).or_else(|| {
        outer.iter().copied().find_map(|component| {
            find(component).filter(|name| {
                name.eq_ignore_ascii_case(component)
                    || (name == UNCLASSIFIED
                        && ["unassigned", "unbarcoded", "default--default"]
                            .iter()
                            .any(|a| a.eq_ignore_ascii_case(component)))
            })
        })
    })
}

/// Name for a PacBio `bc` tag's pair of barcode indices.
///
/// With `names` (the barcode FASTA lima was given) the pair is named as SMRT Link names it,
/// `bc1001--bc1001`; without it, by index. Out-of-range indices mean unassigned.
pub fn from_indices(forward: u16, reverse: u16, names: &[String]) -> String {
    if forward == u16::MAX || reverse == u16::MAX {
        return UNCLASSIFIED.to_owned();
    }
    let name = |index: u16| names.get(usize::from(index)).cloned();
    match (name(forward), name(reverse)) {
        (Some(f), Some(r)) => format!("{f}--{r}"),
        _ if names.is_empty() => format!("{INDEX_PREFIX}{forward}--{INDEX_PREFIX}{reverse}"),
        _ => UNCLASSIFIED.to_owned(),
    }
}

/// The barcode a dorado read group names at its end (`<run>_<model>_barcode03`), if any.
pub fn from_read_group(id: &str) -> Option<String> {
    let name = normalize(id.rsplit('_').next().unwrap_or_default());
    (name.starts_with("barcode") || name == UNCLASSIFIED).then_some(name)
}

/// Caches normalized labels, so each distinct raw label is parsed once per worker.
#[derive(Debug, Default)]
pub struct LabelCache {
    labels: std::collections::HashMap<Vec<u8>, Arc<str>>,
    read_groups: std::collections::HashMap<Vec<u8>, Option<Arc<str>>>,
}

impl LabelCache {
    /// The barcode named by a read group (see [`from_read_group`]), parsed once per group.
    pub fn read_group_barcode(&mut self, id: &[u8]) -> Option<Arc<str>> {
        if let Some(barcode) = self.read_groups.get(id) {
            return barcode.clone();
        }
        let barcode: Option<Arc<str>> =
            from_read_group(&String::from_utf8_lossy(id)).map(Arc::from);
        self.read_groups.insert(id.to_vec(), barcode.clone());
        barcode
    }

    /// The shared form of a raw string, not normalized: run names and the like.
    pub fn intern_verbatim(&mut self, raw: &[u8]) -> Arc<str> {
        let key: Vec<u8> = [b"\0".as_slice(), raw].concat();
        if let Some(label) = self.labels.get(&key) {
            return Arc::clone(label);
        }
        let label: Arc<str> = String::from_utf8_lossy(raw).into();
        self.labels.insert(key, Arc::clone(&label));
        label
    }

    /// The normalized, shared form of a raw label.
    pub fn intern(&mut self, raw: &[u8]) -> Arc<str> {
        if let Some(label) = self.labels.get(raw) {
            return Arc::clone(label);
        }
        let label: Arc<str> = normalize(&String::from_utf8_lossy(raw)).into();
        self.labels.insert(raw.to_vec(), Arc::clone(&label));
        label
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::path::Path;

    #[test]
    fn normalizes_dorado_labels() {
        assert_eq!(normalize("SQK-NBD114-24_barcode01"), "barcode01");
        assert_eq!(normalize("unclassified"), "unclassified");
        assert_eq!(normalize("bc1001--bc1001"), "bc1001--bc1001");
        assert_eq!(normalize("my_sample"), "my_sample");
    }

    #[test]
    fn finds_barcodes_in_paths() {
        let minknow = Path::new("/runs/x/fastq_pass/barcode07/PAY_pass_barcode07_abc_0.fastq.gz");
        assert_eq!(from_path(minknow).as_deref(), Some("barcode07"));
        let lima = Path::new("/data/m84011.hifi_reads.bc2001--bc2001.bam");
        assert_eq!(from_path(lima).as_deref(), Some("bc2001--bc2001"));
        let revio = Path::new("/data/m84011.hifi_reads.bc2045.bam");
        assert_eq!(from_path(revio).as_deref(), Some("bc2045"));
        assert_eq!(
            from_path(Path::new("/data/reads.unassigned.bam")).as_deref(),
            Some(UNCLASSIFIED)
        );
        assert_eq!(from_path(Path::new("/data/abc123/reads.bam")), None);
    }

    #[test]
    fn reads_barcodes_from_dorado_read_groups() {
        let id =
            "7a4eec12-4b56-4de8-984d-756515f320fe_dna_r10.4.1_e8.2_400bps_sup@v5.2.0_barcode03";
        assert_eq!(from_read_group(id).as_deref(), Some("barcode03"));
        assert_eq!(
            from_read_group("abc_dna_sup@v5.2.0_unclassified").as_deref(),
            Some(UNCLASSIFIED)
        );
        assert_eq!(
            from_read_group("70eb6e24_dna_r10.4.1_e8.2_400bps_sup@v5.0.0"),
            None
        );
    }

    #[test]
    fn names_index_pairs() {
        let names = vec!["bc1001".to_owned(), "bc1002".to_owned()];
        assert_eq!(from_indices(1, 1, &names), "bc1002--bc1002");
        assert_eq!(from_indices(0, 1, &[]), "index0--index1");
        assert_eq!(from_indices(5, 5, &names), UNCLASSIFIED);
        assert_eq!(from_indices(u16::MAX, u16::MAX, &[]), UNCLASSIFIED);
    }
}
