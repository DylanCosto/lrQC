//! Sample sheets: MinKNOW's CSV (`barcode`, `alias`) and SMRT Link's biosample CSV
//! (`Barcode`, `Bio Sample`), mapping barcodes to sample names.

use std::collections::HashMap;
use std::fs;
use std::path::Path;

use anyhow::{bail, Context, Result};

use crate::barcode;

/// Column names holding the barcode, in order of preference.
const BARCODE_COLUMNS: [&str; 5] = ["barcode", "barcode name", "barcodes", "bc", "barcode_name"];

/// Column names holding the sample name, in order of preference.
const NAME_COLUMNS: [&str; 8] = [
    "alias",
    "bio sample",
    "bio sample name",
    "biosample",
    "biosample name",
    "sample_name",
    "sample",
    "sample_id",
];

/// Reads a sample sheet into sample names by normalized barcode.
///
/// # Errors
///
/// Fails when the file cannot be read or has no barcode or sample name column.
pub fn load(path: &Path) -> Result<HashMap<String, String>> {
    let text = fs::read_to_string(path).with_context(|| format!("reading {}", path.display()))?;
    parse(&text).with_context(|| format!("reading the sample sheet {}", path.display()))
}

fn parse(text: &str) -> Result<HashMap<String, String>> {
    // Excel's "CSV UTF-8" starts with a byte order mark.
    let text = text.strip_prefix('\u{feff}').unwrap_or(text);
    let mut lines = text.lines().filter(|l| !l.trim().is_empty());
    let Some(header) = lines.next() else {
        bail!("the sample sheet is empty");
    };
    let delimiter = if header.contains('\t') { '\t' } else { ',' };
    let columns: Vec<String> = split_fields(header, delimiter)
        .iter()
        .map(|c| c.to_ascii_lowercase())
        .collect();
    let find = |names: &[&str]| {
        names
            .iter()
            .find_map(|n| columns.iter().position(|c| c == n))
    };
    let Some(barcode_column) = find(&BARCODE_COLUMNS) else {
        bail!("no barcode column (one of {})", BARCODE_COLUMNS.join(", "));
    };
    let Some(name_column) = find(&NAME_COLUMNS) else {
        bail!("no sample name column (one of {})", NAME_COLUMNS.join(", "));
    };
    let mut samples = HashMap::new();
    for line in lines {
        let fields = split_fields(line, delimiter);
        let (Some(code), Some(name)) = (fields.get(barcode_column), fields.get(name_column)) else {
            continue;
        };
        if !code.is_empty() && !name.is_empty() {
            samples.insert(barcode::normalize(code), name.clone());
        }
    }
    Ok(samples)
}

/// Splits one CSV (or TSV) line into trimmed fields: a field in double quotes may hold the
/// delimiter, and `""` inside quotes is one quote, as spreadsheets write them.
fn split_fields(line: &str, delimiter: char) -> Vec<String> {
    let mut fields = Vec::new();
    let mut field = String::new();
    let mut quoted = false;
    let mut chars = line.chars().peekable();
    while let Some(c) = chars.next() {
        match c {
            '"' if quoted && chars.peek() == Some(&'"') => {
                field.push('"');
                chars.next();
            }
            '"' => quoted = !quoted,
            c if c == delimiter && !quoted => fields.push(std::mem::take(&mut field)),
            c => field.push(c),
        }
    }
    fields.push(field);
    fields.into_iter().map(|f| f.trim().to_owned()).collect()
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn reads_minknow_sheets() {
        let sheet = "flow_cell_id,kit,sample_id,experiment_id,barcode,alias\nPAY1,SQK-NBD114-24,s,e,barcode01,patient_A\nPAY1,SQK-NBD114-24,s,e,barcode02,patient_B\n";
        let samples = parse(sheet).expect("valid sheet");
        assert_eq!(samples["barcode01"], "patient_A");
        assert_eq!(samples.len(), 2);
    }

    #[test]
    fn reads_smrt_link_sheets() {
        let sheet = "Barcode,Bio Sample\nbc2001--bc2001,HG002\n\"bc2002--bc2002\",\"HG003\"\n";
        let samples = parse(sheet).expect("valid sheet");
        assert_eq!(samples["bc2002--bc2002"], "HG003");
    }

    #[test]
    fn reads_spreadsheet_csv() {
        let sheet = "\u{feff}Barcode Name,Bio Sample Name\r\nbc2001--bc2001,\"Smith, J\"\r\nbc2002--bc2002,\"the \"\"B\"\" line\"\r\n";
        let samples = parse(sheet).expect("valid sheet");
        assert_eq!(samples["bc2001--bc2001"], "Smith, J");
        assert_eq!(samples["bc2002--bc2002"], "the \"B\" line");
        // A quoted comma before the barcode column must not shift the columns.
        let sheet = "alias,barcode\n\"a, b\",barcode03\n";
        assert_eq!(parse(sheet).expect("valid sheet")["barcode03"], "a, b");
    }

    #[test]
    fn rejects_sheets_without_barcodes() {
        assert!(parse("name,value\na,b\n").is_err());
    }
}
