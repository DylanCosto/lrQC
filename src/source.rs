//! What the decoders know about the file a record came from.

use std::collections::HashMap;
use std::sync::Arc;

use crate::read::{Filter, Platform};
use crate::summary::Columns;

/// Per-file context shared by every record of the file.
#[derive(Debug, Clone, Default)]
pub struct Source {
    /// Platform from the command line or the header; [`Platform::Unknown`] defers to each read.
    pub platform: Platform,
    /// Barcode named by the file's path, used when the records name none.
    pub barcode: Option<Arc<str>>,
    /// Pass/fail named by the file's path.
    pub filter: Option<Filter>,
    /// The BAM has reference sequences, so unmapped records are unmapped reads.
    pub aligned: bool,
    pub read_groups: HashMap<Vec<u8>, GroupInfo>,
    /// PacBio barcode names by index, for `bc` tags.
    pub barcode_names: Arc<Vec<String>>,
    /// Column positions of a sequencing summary.
    pub columns: Columns,
    /// Run named by the file (MinKNOW's `<flow cell>_pass_...` file names).
    pub run: Option<Arc<str>>,
}

/// What a read group says about its reads.
#[derive(Debug, Clone, Default)]
pub struct GroupInfo {
    pub barcode: Option<String>,
    pub platform: Platform,
    /// Flow cell or movie (`PU`), else the run ID at the front of the group ID.
    pub run: Option<Arc<str>>,
    pub molecule: crate::read::Molecule,
}

impl Source {
    #[cfg(test)]
    pub fn for_tests(platform: Platform) -> Self {
        Self {
            platform,
            ..Self::default()
        }
    }
}
