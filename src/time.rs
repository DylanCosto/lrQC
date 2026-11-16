//! Parsing the ISO 8601 read start times dorado writes (`st:Z:2024-08-20T10:39:34.457+00:00`).

/// Seconds since the Unix epoch for an ISO 8601 date-time with optional fraction and offset.
pub fn parse_iso8601(text: &[u8]) -> Option<f64> {
    let text = std::str::from_utf8(text).ok()?;
    let (date, rest) = text.split_once(['T', ' '])?;
    let mut date_parts = date.splitn(3, '-');
    let year: i64 = date_parts.next()?.parse().ok()?;
    let month: u32 = date_parts.next()?.parse().ok()?;
    let day: u32 = date_parts.next()?.parse().ok()?;

    let offset_start = rest
        .char_indices()
        .skip(1)
        .find(|&(_, c)| c == 'Z' || c == '+' || c == '-')
        .map_or(rest.len(), |(i, _)| i);
    let (clock, zone) = rest.split_at(offset_start);
    let mut clock_parts = clock.splitn(3, ':');
    let hour: f64 = clock_parts.next()?.parse().ok()?;
    let minute: f64 = clock_parts.next()?.parse().ok()?;
    let second: f64 = clock_parts.next().unwrap_or("0").parse().ok()?;

    let offset_secs = match zone {
        "" | "Z" => 0.0,
        _ => {
            let sign = if zone.starts_with('-') { -1.0 } else { 1.0 };
            let digits = &zone[1..];
            let (h, m) = digits
                .split_once(':')
                .unwrap_or_else(|| digits.split_at(digits.len().min(2)));
            let h: f64 = h.parse().ok()?;
            let m: f64 = if m.is_empty() { 0.0 } else { m.parse().ok()? };
            sign * (h * 3600.0 + m * 60.0)
        }
    };
    // Day counts fit an f64 exactly.
    let days = days_from_civil(year, month, day)? as f64;
    Some(days * 86_400.0 + hour * 3600.0 + minute * 60.0 + second - offset_secs)
}

/// Unix seconds as an ISO 8601 UTC date-time with milliseconds.
pub fn format_iso8601(secs: f64) -> String {
    let millis = (secs * 1000.0).round() as i64;
    let (days, ms_of_day) = (millis.div_euclid(86_400_000), millis.rem_euclid(86_400_000));
    let (year, month, day) = civil_from_days(days);
    let (h, rest) = (ms_of_day / 3_600_000, ms_of_day % 3_600_000);
    let (m, rest) = (rest / 60_000, rest % 60_000);
    let (s, ms) = (rest / 1000, rest % 1000);
    format!("{year:04}-{month:02}-{day:02}T{h:02}:{m:02}:{s:02}.{ms:03}Z")
}

/// The date of a day count since 1970-01-01 (Howard Hinnant's algorithm).
fn civil_from_days(days: i64) -> (i64, u32, u32) {
    let z = days + 719_468;
    let era = z.div_euclid(146_097);
    let day_of_era = z - era * 146_097;
    let year_of_era =
        (day_of_era - day_of_era / 1460 + day_of_era / 36_524 - day_of_era / 146_096) / 365;
    let day_of_year = day_of_era - (365 * year_of_era + year_of_era / 4 - year_of_era / 100);
    let shifted_month = (5 * day_of_year + 2) / 153;
    let day = (day_of_year - (153 * shifted_month + 2) / 5 + 1) as u32;
    let month = if shifted_month < 10 {
        shifted_month + 3
    } else {
        shifted_month - 9
    } as u32;
    let year = year_of_era + era * 400 + i64::from(month <= 2);
    (year, month, day)
}

/// Days since 1970-01-01 (Howard Hinnant's algorithm).
fn days_from_civil(year: i64, month: u32, day: u32) -> Option<i64> {
    if !(1..=12).contains(&month) || !(1..=31).contains(&day) {
        return None;
    }
    let year = if month <= 2 { year - 1 } else { year };
    let era = year.div_euclid(400);
    let year_of_era = year - era * 400;
    let month = i64::from(month);
    let shifted_month = if month > 2 { month - 3 } else { month + 9 };
    let day_of_year = (153 * shifted_month + 2) / 5 + i64::from(day) - 1;
    let day_of_era = year_of_era * 365 + year_of_era / 4 - year_of_era / 100 + day_of_year;
    Some(era * 146_097 + day_of_era - 719_468)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn parses_dorado_start_times() {
        let t = parse_iso8601(b"2024-08-20T10:39:34.457+00:00").expect("valid time");
        assert!((t - 1_724_150_374.457).abs() < 1e-3, "got {t}");
        let shifted = parse_iso8601(b"2024-08-20T12:39:34.457+02:00").expect("valid time");
        assert!((shifted - t).abs() < 1e-3);
        let zulu = parse_iso8601(b"1970-01-02T00:00:00Z").expect("valid time");
        assert!((zulu - 86_400.0).abs() < 1e-9);
        assert_eq!(parse_iso8601(b"not a time"), None);
    }

    #[test]
    fn formats_round_trip() {
        for text in [
            "2024-08-20T10:39:34.457Z",
            "1999-12-31T23:59:59.999Z",
            "2000-02-29T00:00:00.000Z",
        ] {
            let secs = parse_iso8601(text.as_bytes()).expect("valid time");
            assert_eq!(format_iso8601(secs), text);
        }
    }
}
