//! Warn and fail levels of every check, overridable from a TOML file.

use std::fs;
use std::path::Path;

use anyhow::{Context, Result};
use serde::{Deserialize, Serialize};

/// Keys of checks lrqc no longer has, accepted in older files and ignored.
const RETIRED_KEYS: &[&str] = &["pass_rate_warn", "pass_rate_fail"];

/// Warn and fail levels. Percentages are of reads unless named otherwise.
#[derive(Debug, Clone, Serialize, Deserialize)]
#[serde(default, deny_unknown_fields)]
pub struct Thresholds {
    /// Read Q warn level, by platform: MinKNOW's super-accuracy pass threshold, HiFi's Q20.
    pub read_q_warn_ont: f64,
    pub read_q_warn_pacbio: f64,
    /// The read Q fail level as a share of the warn level.
    pub read_q_fail_ratio: f64,
    /// Largest A/T or G/C difference at a position, in percentage points.
    pub base_content_warn: f64,
    pub base_content_fail: f64,
    /// Reads under a second GC peak.
    pub gc_second_peak_warn: f64,
    pub gc_second_peak_fail: f64,
    /// N bases at a position.
    pub n_content_warn: f64,
    pub n_content_fail: f64,
    /// Sampled reads that duplicate another.
    pub duplication_warn: f64,
    pub duplication_fail: f64,
    /// Sampled reads sharing one 50-base end.
    pub overrepresented_warn: f64,
    pub overrepresented_fail: f64,
    /// Reads with an adapter within 150 bases of an end, and inside reads.
    pub adapter_end_warn: f64,
    pub adapter_end_fail: f64,
    pub adapter_internal_warn: f64,
    pub adapter_internal_fail: f64,
    /// Reads repeating an earlier read ID; any is a warning.
    pub duplicate_ids_fail: f64,
    /// A barcode's reads against the median barcode's, as a ratio.
    pub barcode_ratio_warn: f64,
    pub barcode_ratio_fail: f64,
    /// Unclassified reads, as a share of the run.
    pub unclassified_warn: f64,
    pub unclassified_fail: f64,
    /// Reads matching a contaminant reference.
    pub contaminant_warn: f64,
    pub contaminant_fail: f64,
    /// Reads still holding their own barcode at an end.
    pub untrimmed_barcode_warn: f64,
    pub untrimmed_barcode_fail: f64,
    /// Reads holding another sample's barcode.
    pub foreign_barcode_warn: f64,
    pub foreign_barcode_fail: f64,
    /// Primary reads mapped, for aligned input.
    pub mapped_warn: f64,
    pub mapped_fail: f64,
}

impl Default for Thresholds {
    fn default() -> Self {
        Self {
            read_q_warn_ont: 10.0,
            read_q_warn_pacbio: 20.0,
            read_q_fail_ratio: 0.75,
            base_content_warn: 10.0,
            base_content_fail: 20.0,
            gc_second_peak_warn: 5.0,
            gc_second_peak_fail: 20.0,
            n_content_warn: 5.0,
            n_content_fail: 20.0,
            duplication_warn: 20.0,
            duplication_fail: 50.0,
            overrepresented_warn: 0.1,
            overrepresented_fail: 1.0,
            adapter_end_warn: 5.0,
            adapter_end_fail: 10.0,
            adapter_internal_warn: 2.0,
            adapter_internal_fail: 5.0,
            duplicate_ids_fail: 1.0,
            barcode_ratio_warn: 0.25,
            barcode_ratio_fail: 0.1,
            unclassified_warn: 10.0,
            unclassified_fail: 25.0,
            contaminant_warn: 1.0,
            contaminant_fail: 5.0,
            untrimmed_barcode_warn: 5.0,
            untrimmed_barcode_fail: 20.0,
            foreign_barcode_warn: 1.0,
            foreign_barcode_fail: 5.0,
            mapped_warn: 90.0,
            mapped_fail: 70.0,
        }
    }
}

impl Thresholds {
    /// Defaults overridden by the keys of a TOML file.
    ///
    /// # Errors
    ///
    /// Fails when the file cannot be read, is not TOML, or names an unknown key. Keys of
    /// checks that were removed are ignored with a warning.
    pub fn load(path: &Path) -> Result<Self> {
        let text =
            fs::read_to_string(path).with_context(|| format!("reading {}", path.display()))?;
        Self::parse(&text).with_context(|| format!("parsing thresholds in {}", path.display()))
    }

    fn parse(text: &str) -> Result<Self> {
        let mut table: toml::Table = toml::from_str(text)?;
        for key in RETIRED_KEYS {
            if table.remove(*key).is_some() {
                eprintln!("lrqc: ignoring threshold {key}: its check was removed");
            }
        }
        let thresholds: Self = table.try_into()?;
        thresholds.validate()?;
        Ok(thresholds)
    }

    /// Rejects levels that cannot grade: a negative or non-finite value, a fail level that
    /// is not worse than its warn level, a ratio outside 0 to 1.
    fn validate(&self) -> Result<()> {
        let higher_is_worse = [
            (
                "base_content",
                self.base_content_warn,
                self.base_content_fail,
            ),
            (
                "gc_second_peak",
                self.gc_second_peak_warn,
                self.gc_second_peak_fail,
            ),
            ("n_content", self.n_content_warn, self.n_content_fail),
            ("duplication", self.duplication_warn, self.duplication_fail),
            (
                "overrepresented",
                self.overrepresented_warn,
                self.overrepresented_fail,
            ),
            ("adapter_end", self.adapter_end_warn, self.adapter_end_fail),
            (
                "adapter_internal",
                self.adapter_internal_warn,
                self.adapter_internal_fail,
            ),
            (
                "unclassified",
                self.unclassified_warn,
                self.unclassified_fail,
            ),
            ("contaminant", self.contaminant_warn, self.contaminant_fail),
            (
                "untrimmed_barcode",
                self.untrimmed_barcode_warn,
                self.untrimmed_barcode_fail,
            ),
            (
                "foreign_barcode",
                self.foreign_barcode_warn,
                self.foreign_barcode_fail,
            ),
        ];
        let lower_is_worse = [
            (
                "barcode_ratio",
                self.barcode_ratio_warn,
                self.barcode_ratio_fail,
            ),
            ("mapped", self.mapped_warn, self.mapped_fail),
        ];
        let levels = [
            ("read_q_warn_ont", self.read_q_warn_ont),
            ("read_q_warn_pacbio", self.read_q_warn_pacbio),
            ("duplicate_ids_fail", self.duplicate_ids_fail),
        ];
        let pairs =
            higher_is_worse
                .iter()
                .chain(&lower_is_worse)
                .flat_map(|&(name, warn, fail)| {
                    [
                        (format!("{name}_warn"), warn),
                        (format!("{name}_fail"), fail),
                    ]
                });
        for (key, value) in levels
            .into_iter()
            .map(|(key, value)| (key.to_owned(), value))
            .chain(pairs)
        {
            anyhow::ensure!(
                value.is_finite() && value >= 0.0,
                "threshold {key} must be a number of 0 or more, not {value}"
            );
        }
        for (name, warn, fail) in higher_is_worse {
            anyhow::ensure!(
                warn <= fail,
                "{name}_warn ({warn}) must not be above {name}_fail ({fail}): warn is the lower level"
            );
        }
        for (name, warn, fail) in lower_is_worse {
            anyhow::ensure!(
                warn >= fail,
                "{name}_warn ({warn}) must not be below {name}_fail ({fail}): here lower is worse"
            );
        }
        anyhow::ensure!(
            self.read_q_fail_ratio > 0.0 && self.read_q_fail_ratio <= 1.0,
            "read_q_fail_ratio must be above 0 and at most 1, not {}",
            self.read_q_fail_ratio
        );
        Ok(())
    }

    /// The defaults as TOML, for `--print-thresholds`.
    pub fn to_toml(&self) -> String {
        toml::to_string(self).expect("thresholds serialize")
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn rejects_levels_that_cannot_grade() {
        assert!(Thresholds::parse("contaminant_warn = 10.0\ncontaminant_fail = 5.0\n").is_err());
        assert!(
            Thresholds::parse("mapped_warn = 50.0\n").is_err(),
            "below mapped_fail"
        );
        assert!(Thresholds::parse("read_q_fail_ratio = 1.5\n").is_err());
        assert!(Thresholds::parse("n_content_warn = -1.0\n").is_err());
        assert!(Thresholds::parse("contaminant_warn = 0.1\ncontaminant_fail = 1.0\n").is_ok());
    }

    #[test]
    fn overrides_only_the_keys_given() {
        let t: Thresholds = toml::from_str("read_q_warn_ont = 12.0\n").expect("valid TOML");
        assert_eq!(t.read_q_warn_ont, 12.0);
        assert_eq!(
            t.read_q_warn_pacbio,
            Thresholds::default().read_q_warn_pacbio
        );
        assert!(toml::from_str::<Thresholds>("no_such_key = 1.0\n").is_err());
    }

    #[test]
    fn ignores_keys_of_removed_checks() {
        let t = Thresholds::parse("pass_rate_warn = 90.0\nread_q_warn_ont = 12.0\n")
            .expect("retired keys are ignored");
        assert_eq!(t.read_q_warn_ont, 12.0);
        assert!(Thresholds::parse("no_such_key = 1.0\n").is_err());
    }

    #[test]
    fn round_trips_through_toml() {
        let text = Thresholds::default().to_toml();
        let back: Thresholds = toml::from_str(&text).expect("valid TOML");
        assert_eq!(back.contaminant_fail, 5.0);
    }
}
