//! Run information from BAM headers: flow cells and movies, basecallers, models and kits.

use serde::{Deserialize, Serialize};

use crate::pipeline::FileInfo;

/// What the BAM headers say about how the reads were made.
#[derive(Debug, Clone, Default, Serialize, Deserialize)]
pub struct RunInfo {
    pub read_groups: Vec<ReadGroupInfo>,
    /// Programs that touched the reads: name and version.
    pub programs: Vec<(String, String)>,
}

/// One `@RG` line, with the fields people look for named.
#[derive(Debug, Clone, Default, Serialize, Deserialize)]
pub struct ReadGroupInfo {
    pub id: String,
    pub platform: Option<String>,
    pub sample: Option<String>,
    /// ONT flow cell or PacBio movie (`PU`).
    pub flow_cell: Option<String>,
    /// ONT device position (`PM`).
    pub device: Option<String>,
    pub started: Option<String>,
    pub basecall_model: Option<String>,
    pub modbase_models: Option<String>,
    /// PacBio read type, binding and sequencing kits and basecaller version.
    pub pacbio: Option<String>,
}

impl RunInfo {
    /// Collects the distinct read groups and programs of the inputs' headers.
    pub fn from_inputs(files: &[FileInfo]) -> Self {
        let mut info = RunInfo::default();
        for line in files.iter().flat_map(|f| &f.header_lines) {
            let fields: Vec<(&str, &str)> = line
                .split('\t')
                .skip(1)
                .filter_map(|f| f.split_once(':'))
                .collect();
            let field = |key: &str| {
                fields
                    .iter()
                    .find(|(k, _)| *k == key)
                    .map(|(_, v)| (*v).to_owned())
            };
            if line.starts_with("@RG") {
                let description = field("DS").unwrap_or_default();
                let ds = |key: &str| {
                    description
                        .split([' ', ';'])
                        .find_map(|part| part.strip_prefix(key)?.strip_prefix('='))
                        .map(str::to_owned)
                };
                let pacbio: Vec<String> = [
                    "READTYPE",
                    "BINDINGKIT",
                    "SEQUENCINGKIT",
                    "BASECALLERVERSION",
                ]
                .iter()
                .filter_map(|key| ds(key).map(|v| format!("{key}={v}")))
                .collect();
                let group = ReadGroupInfo {
                    id: field("ID").unwrap_or_default(),
                    platform: field("PL"),
                    sample: field("SM"),
                    flow_cell: field("PU"),
                    device: field("PM"),
                    started: field("DT"),
                    basecall_model: ds("basecall_model"),
                    modbase_models: ds("modbase_models"),
                    pacbio: (!pacbio.is_empty()).then(|| pacbio.join(", ")),
                };
                if !info.read_groups.iter().any(|g| g.id == group.id) {
                    info.read_groups.push(group);
                }
            } else if line.starts_with("@PG") {
                let name = field("PN").or_else(|| field("ID")).unwrap_or_default();
                let version = field("VN").unwrap_or_default();
                let program = (name, version);
                if !program.0.is_empty() && !info.programs.contains(&program) {
                    info.programs.push(program);
                }
            }
        }
        info
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::input::Format;

    #[test]
    fn names_dorado_and_pacbio_fields() {
        let file = FileInfo {
            path: "x.bam".into(),
            format: Format::Bam,
            records: 0,
            path_barcode: None,
            header_lines: vec![
                "@RG\tID:abc_dna_sup\tPU:PAY87954\tPM:PC48B260\tDT:2024-08-19T15:21:25\tPL:ONT\tDS:basecall_model=dna_r10.4.1_e8.2_400bps_sup@v5.0.0 modbase_models=x runid=abc".into(),
                "@RG\tID:a57306fa\tPL:PACBIO\tDS:READTYPE=CCS;BINDINGKIT=102-739-100;SEQUENCINGKIT=102-118-800;BASECALLERVERSION=5.0\tPU:m84011_220902_175841_s1".into(),
                "@PG\tID:basecaller\tPN:dorado\tVN:0.8.0".into(),
                "@PG\tID:basecaller\tPN:dorado\tVN:0.8.0".into(),
            ],
        };
        let info = RunInfo::from_inputs(&[file]);
        assert_eq!(info.read_groups.len(), 2);
        let ont = &info.read_groups[0];
        assert_eq!(ont.flow_cell.as_deref(), Some("PAY87954"));
        assert_eq!(
            ont.basecall_model.as_deref(),
            Some("dna_r10.4.1_e8.2_400bps_sup@v5.0.0")
        );
        assert!(info.read_groups[1]
            .pacbio
            .as_deref()
            .is_some_and(|p| p.contains("BINDINGKIT=102-739-100")));
        assert_eq!(
            info.programs,
            vec![("dorado".to_owned(), "0.8.0".to_owned())]
        );
    }
}
