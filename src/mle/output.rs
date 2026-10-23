//! Deterministic serialization of MLE results to legacy TSV format.

use super::gene_result::GeneResult;
use super::serial::SerialGeneInput;
use crate::error::MageckError;
use std::fmt::Write as _;
use std::io::Write;

fn math(message: impl Into<String>) -> MageckError {
    MageckError::Math(message.into())
}

fn check_field(field: &str) -> Result<(), MageckError> {
    if field.contains(['\t', '\r', '\n']) {
        return Err(math(
            "serial MLE output identifier contains invalid delimiter",
        ));
    }
    Ok(())
}

fn format_general(value: f64, precision: usize) -> Result<String, MageckError> {
    if !value.is_finite() {
        return Err(math(format!(
            "serial MLE output requires finite numerical value, got {value}"
        )));
    }
    if value == 0.0 {
        return Ok(if value.is_sign_negative() { "-0" } else { "0" }.to_owned());
    }

    let scientific = format!("{value:.digits$e}", digits = precision - 1);
    let (mantissa, exponent) = scientific
        .split_once('e')
        .ok_or_else(|| math("standard scientific formatter omitted exponent"))?;
    let exponent: i32 = exponent
        .parse()
        .map_err(|_| math("standard scientific formatter produced invalid exponent"))?;
    let negative = mantissa.starts_with('-');
    let digits = mantissa.trim_start_matches('-').replace('.', "");
    if digits.len() != precision || !digits.bytes().all(|byte| byte.is_ascii_digit()) {
        return Err(math(
            "standard scientific formatter produced invalid significand",
        ));
    }
    let mut rendered = if (-4..(precision as i32)).contains(&exponent) {
        let decimal = exponent + 1;
        if decimal <= 0 {
            format!("0.{}{}", "0".repeat((-decimal) as usize), digits)
        } else if decimal as usize >= digits.len() {
            format!("{}{}", digits, "0".repeat(decimal as usize - digits.len()))
        } else {
            format!(
                "{}.{}",
                &digits[..decimal as usize],
                &digits[decimal as usize..]
            )
        }
    } else {
        format!("{}.{}e{exponent:+03}", &digits[..1], &digits[1..])
    };
    if let Some(point) = rendered.find('.') {
        let exponent_start = rendered.find('e').unwrap_or(rendered.len());
        let fractional = &rendered[point + 1..exponent_start];
        let trimmed = fractional.trim_end_matches('0').to_owned();
        if trimmed.len() != fractional.len() {
            let suffix = rendered[exponent_start..].to_owned();
            rendered.truncate(point);
            if !trimmed.is_empty() {
                rendered.push('.');
                rendered.push_str(&trimmed);
            }
            rendered.push_str(&suffix);
        }
    }
    Ok(if negative {
        format!("-{rendered}")
    } else {
        rendered
    })
}

pub(crate) fn format_5g(value: f64) -> Result<String, MageckError> {
    format_general(value, 5)
}

pub(crate) fn format_3g(value: f64) -> Result<String, MageckError> {
    format_general(value, 3)
}

fn validate(genes: &[SerialGeneInput], results: &[GeneResult]) -> Result<(), MageckError> {
    if genes.len() != results.len() {
        return Err(math(
            "serial MLE output gene input and result length mismatch",
        ));
    }
    let reference = results.first().map(|result| &result.condition_order);
    for (index, (gene, result)) in genes.iter().zip(results).enumerate() {
        if !result.status.converged() {
            return Err(math(
                "MLE primary serialization requires eligible converged fit",
            ));
        }
        check_field(&gene.gene_id)?;
        check_field(&result.gene_id)?;
        for guide in &gene.sgrna_ids {
            check_field(guide)?;
        }
        if gene.gene_id != result.gene_id {
            return Err(math(format!(
                "serial MLE output gene ID mismatch at index {index}: {} vs {}",
                gene.gene_id, result.gene_id
            )));
        }
        let conditions = result.condition_order.len();
        if reference != Some(&result.condition_order)
            || result.condition_beta.len() != conditions
            || result.z_score.len() != conditions
            || result.two_sided_p_value.len() != conditions
            || result.condition_fdr.len() != conditions
        {
            return Err(math("serial MLE output condition alignment mismatch"));
        }
        for label in &result.condition_order {
            check_field(label)?;
        }
        if let Some(efficiencies) = &result.efficiencies {
            if efficiencies.len() != gene.sgrna_ids.len() {
                return Err(math("serial MLE output efficiency length mismatch"));
            }
            for &value in efficiencies {
                format_3g(value)?;
            }
        }
        for &value in result.condition_beta.iter().chain(&result.z_score) {
            format_5g(value)?;
        }
        for &value in result.two_sided_p_value.iter().chain(&result.condition_fdr) {
            if !(0.0..=1.0).contains(&value) {
                return Err(math(format!(
                    "probability value out of [0, 1] domain: {value}"
                )));
            }
            format_5g(value)?;
        }
    }
    Ok(())
}

fn gene_string(genes: &[SerialGeneInput], results: &[GeneResult]) -> Result<String, MageckError> {
    validate(genes, results)?;
    let mut output = String::from("Gene\tsgRNA");
    if let Some(first) = results.first() {
        for condition in &first.condition_order {
            write!(output, "\t{condition}|beta\t{condition}|z\t{condition}|p-value\t{condition}|fdr\t{condition}|wald-p-value\t{condition}|wald-fdr").unwrap();
        }
    }
    output.push('\n');
    for (gene, result) in genes.iter().zip(results) {
        write!(output, "{}\t{}", result.gene_id, gene.sgrna_ids.len()).unwrap();
        for index in 0..result.condition_order.len() {
            write!(
                output,
                "\t{}\t{}\t1.0\t1.0\t{}\t{}",
                format_5g(result.condition_beta[index])?,
                format_5g(result.z_score[index])?,
                format_5g(result.two_sided_p_value[index])?,
                format_5g(result.condition_fdr[index])?
            )
            .unwrap();
        }
        output.push('\n');
    }
    Ok(output)
}

fn sgrna_string(genes: &[SerialGeneInput], results: &[GeneResult]) -> Result<String, MageckError> {
    validate(genes, results)?;
    let mut output = String::from("Gene\tsgRNA\teff\n");
    for (gene, result) in genes.iter().zip(results) {
        for (index, guide) in gene.sgrna_ids.iter().enumerate() {
            let efficiency = match &result.efficiencies {
                Some(values) => format_3g(values[index])?,
                None => "1".to_owned(),
            };
            writeln!(output, "{}\t{guide}\t{efficiency}", gene.gene_id).unwrap();
        }
    }
    Ok(output)
}

pub(crate) fn serialize_gene_summary_to_string(
    genes: &[SerialGeneInput],
    results: &[GeneResult],
) -> Result<String, MageckError> {
    gene_string(genes, results)
}
pub(crate) fn serialize_sgrna_summary_to_string(
    genes: &[SerialGeneInput],
    results: &[GeneResult],
) -> Result<String, MageckError> {
    sgrna_string(genes, results)
}
pub(crate) fn serialize_gene_summary<W: Write>(
    genes: &[SerialGeneInput],
    results: &[GeneResult],
    writer: &mut W,
) -> Result<(), MageckError> {
    writer.write_all(gene_string(genes, results)?.as_bytes())?;
    Ok(())
}
pub(crate) fn serialize_sgrna_summary<W: Write>(
    genes: &[SerialGeneInput],
    results: &[GeneResult],
    writer: &mut W,
) -> Result<(), MageckError> {
    writer.write_all(sgrna_string(genes, results)?.as_bytes())?;
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::{
        format_3g, format_5g, serialize_gene_summary, serialize_gene_summary_to_string,
        serialize_sgrna_summary, serialize_sgrna_summary_to_string,
    };
    use crate::error::MageckError;
    use crate::mle::{gene_result::GeneResult, irls::IrlsStatus, serial::SerialGeneInput};
    use std::io::{self, Write};

    fn input(id: &str, guides: &[&str]) -> SerialGeneInput {
        SerialGeneInput {
            gene_id: id.into(),
            sgrna_ids: guides.iter().map(|v| (*v).into()).collect(),
            counts: vec![],
        }
    }
    fn result(
        id: &str,
        labels: &[&str],
        beta: Vec<f64>,
        z: Vec<f64>,
        p: Vec<f64>,
        fdr: Vec<f64>,
        efficiencies: Option<Vec<f64>>,
    ) -> GeneResult {
        GeneResult {
            gene_id: id.into(),
            condition_order: labels.iter().map(|v| (*v).into()).collect(),
            condition_beta: beta,
            condition_covariance: vec![],
            standard_error: vec![],
            z_score: z,
            two_sided_p_value: p,
            condition_fdr: fdr,
            status: IrlsStatus::Converged,
            attempts: 1,
            efficiencies,
            means: vec![],
            residuals: vec![],
        }
    }

    #[test]
    fn formatter_oracle_transition_and_extreme_vectors() {
        let vectors = [
            (0x0000000000000000, "0", "0"),
            (0x8000000000000000, "-0", "-0"),
            (0x3f1a36e1334e458f, "0.0001", "0.0001"),
            (0x408f3c0000000000, "999.5", "1e+03"),
            (0x40f869f800000000, "1e+05", "1e+05"),
            (0xc08f3c0000000000, "-999.5", "-1e+03"),
            (0xc0f869f800000000, "-1e+05", "-1e+05"),
            (0x0000000000000001, "4.9407e-324", "4.94e-324"),
            (0x0010000000000000, "2.2251e-308", "2.23e-308"),
            (0x7fefffffffffffff, "1.7977e+308", "1.8e+308"),
            (0x3fd2b73f5786af52, "0.29243", "0.292"),
            (0x3fe51d6d642cd43e, "0.65984", "0.66"),
        ];
        for (bits, five, three) in vectors {
            let value = f64::from_bits(bits);
            assert_eq!(format_5g(value).unwrap(), five, "{bits:#018x}");
            assert_eq!(format_3g(value).unwrap(), three, "{bits:#018x}");
        }
    }

    #[test]
    fn complete_authoritative_formatter_oracle_executes_production_helpers() {
        // Hermetic direct oracle: stable IDs, exact f64 bits, and expected output.
        const ORACLE: &str = r#"
OUT1-FMT-5G-001|OUT1-FMT-3G-001|0x0000000000000000|0|0
OUT1-FMT-5G-002|OUT1-FMT-3G-002|0x8000000000000000|-0|-0
OUT1-FMT-5G-003|OUT1-FMT-3G-003|0x3f1a36d1bd105b06|9.9999e-05|0.0001
OUT1-FMT-5G-004|OUT1-FMT-3G-004|0x3f1a36da54164f19|9.9999e-05|0.0001
OUT1-FMT-5G-005|OUT1-FMT-3G-005|0x3f1a36e1334e458f|0.0001|0.0001
OUT1-FMT-5G-006|OUT1-FMT-3G-006|0x3f1a36e2eb1c432d|0.0001|0.0001
OUT1-FMT-5G-007|OUT1-FMT-3G-007|0x3f1a36f419282b53|0.0001|0.0001
OUT1-FMT-5G-008|OUT1-FMT-3G-008|0x3f1a378eb79354b1|0.00010001|0.0001
OUT1-FMT-5G-009|OUT1-FMT-3G-009|0x3f1a302cee759401|9.99e-05|9.99e-05
OUT1-FMT-5G-010|OUT1-FMT-3G-010|0x3f1a3387ecc8eb97|9.995e-05|0.0001
OUT1-FMT-5G-011|OUT1-FMT-3G-011|0x408f380000000000|999|999
OUT1-FMT-5G-012|OUT1-FMT-3G-012|0x408f3b3333333333|999.4|999
OUT1-FMT-5G-013|OUT1-FMT-3G-013|0x408f3c0000000000|999.5|1e+03
OUT1-FMT-5G-014|OUT1-FMT-3G-014|0x408f3f3333333333|999.9|1e+03
OUT1-FMT-5G-015|OUT1-FMT-3G-015|0x408f400000000000|1000|1e+03
OUT1-FMT-5G-016|OUT1-FMT-3G-016|0x408f40cccccccccd|1000.1|1e+03
OUT1-FMT-5G-017|OUT1-FMT-3G-017|0x408f440000000000|1000.5|1e+03
OUT1-FMT-5G-018|OUT1-FMT-3G-018|0x40f869f000000000|99999|1e+05
OUT1-FMT-5G-019|OUT1-FMT-3G-019|0x40f869f666666666|99999|1e+05
OUT1-FMT-5G-020|OUT1-FMT-3G-020|0x40f869f800000000|1e+05|1e+05
OUT1-FMT-5G-021|OUT1-FMT-3G-021|0x40f869fe66666666|1e+05|1e+05
OUT1-FMT-5G-022|OUT1-FMT-3G-022|0x40f86a0000000000|1e+05|1e+05
OUT1-FMT-5G-023|OUT1-FMT-3G-023|0x40f86a019999999a|1e+05|1e+05
OUT1-FMT-5G-024|OUT1-FMT-3G-024|0x40f86a0800000000|1e+05|1e+05
OUT1-FMT-5G-025|OUT1-FMT-3G-025|0xc08f3c0000000000|-999.5|-1e+03
OUT1-FMT-5G-026|OUT1-FMT-3G-026|0xc0f869f800000000|-1e+05|-1e+05
OUT1-FMT-5G-027|OUT1-FMT-3G-027|0xbf1a36e1334e458f|-0.0001|-0.0001
OUT1-FMT-5G-028|OUT1-FMT-3G-028|0x3ee4f8b588e368f1|1e-05|1e-05
OUT1-FMT-5G-029|OUT1-FMT-3G-029|0x3f1a36e2eb1c432d|0.0001|0.0001
OUT1-FMT-5G-030|OUT1-FMT-3G-030|0x3f50624dd2f1a9fc|0.001|0.001
OUT1-FMT-5G-031|OUT1-FMT-3G-031|0x3f847ae147ae147b|0.01|0.01
OUT1-FMT-5G-032|OUT1-FMT-3G-032|0x3fb999999999999a|0.1|0.1
OUT1-FMT-5G-033|OUT1-FMT-3G-033|0x3ff0000000000000|1|1
OUT1-FMT-5G-034|OUT1-FMT-3G-034|0x4024000000000000|10|10
OUT1-FMT-5G-035|OUT1-FMT-3G-035|0x4059000000000000|100|100
OUT1-FMT-5G-036|OUT1-FMT-3G-036|0x408f400000000000|1000|1e+03
OUT1-FMT-5G-037|OUT1-FMT-3G-037|0x40c3880000000000|10000|1e+04
OUT1-FMT-5G-038|OUT1-FMT-3G-038|0x40f86a0000000000|1e+05|1e+05
OUT1-FMT-5G-039|OUT1-FMT-3G-039|0x412e848000000000|1e+06|1e+06
OUT1-FMT-5G-040|OUT1-FMT-3G-040|0xbee4f8b588e368f1|-1e-05|-1e-05
OUT1-FMT-5G-041|OUT1-FMT-3G-041|0xbf1a36e2eb1c432d|-0.0001|-0.0001
OUT1-FMT-5G-042|OUT1-FMT-3G-042|0xbff0000000000000|-1|-1
OUT1-FMT-5G-043|OUT1-FMT-3G-043|0xc059000000000000|-100|-100
OUT1-FMT-5G-044|OUT1-FMT-3G-044|0xc0f86a0000000000|-1e+05|-1e+05
OUT1-FMT-5G-045|OUT1-FMT-3G-045|0x3fe0000000000000|0.5|0.5
OUT1-FMT-5G-046|OUT1-FMT-3G-046|0xbfe0000000000000|-0.5|-0.5
OUT1-FMT-5G-047|OUT1-FMT-3G-047|0x3fbf9aae6c8f7553|0.12345|0.123
OUT1-FMT-5G-048|OUT1-FMT-3G-048|0x3fbf9abf33871609|0.12345|0.123
OUT1-FMT-5G-049|OUT1-FMT-3G-049|0x3fbf9acffa7eb6bf|0.12346|0.123
OUT1-FMT-5G-050|OUT1-FMT-3G-050|0x3fbf9db22d0e5604|0.1235|0.123
OUT1-FMT-5G-051|OUT1-FMT-3G-051|0x3fbfdf3b645a1cac|0.1245|0.124
OUT1-FMT-5G-052|OUT1-FMT-3G-052|0x3ff8000000000000|1.5|1.5
OUT1-FMT-5G-053|OUT1-FMT-3G-053|0x3ff3333333333333|1.2|1.2
OUT1-FMT-5G-054|OUT1-FMT-3G-054|0x3f50624dd2f1a9fc|0.001|0.001
OUT1-FMT-5G-055|OUT1-FMT-3G-055|0x0010000000000000|2.2251e-308|2.23e-308
OUT1-FMT-5G-056|OUT1-FMT-3G-056|0x0000000000000001|4.9407e-324|4.94e-324
OUT1-FMT-5G-057|OUT1-FMT-3G-057|0x0000000000000003|1.4822e-323|1.48e-323
OUT1-FMT-5G-058|OUT1-FMT-3G-058|0x7fe1ccf385ebc8a0|1e+308|1e+308
OUT1-FMT-5G-059|OUT1-FMT-3G-059|0x7fefffffffffffff|1.7977e+308|1.8e+308
OUT1-FMT-5G-060|OUT1-FMT-3G-060|0x3fd2b73f5786af52|0.29243|0.292
OUT1-FMT-5G-061|OUT1-FMT-3G-061|0xbfc596eca2b057f4|-0.16867|-0.169
OUT1-FMT-5G-062|OUT1-FMT-3G-062|0x3ff9753a6e0a3e16|1.5911|1.59
OUT1-FMT-5G-063|OUT1-FMT-3G-063|0xbfed5d60a04e4a52|-0.91765|-0.918
OUT1-FMT-5G-064|OUT1-FMT-3G-064|0x3fbc90ac40afca0b|0.11158|0.112
OUT1-FMT-5G-065|OUT1-FMT-3G-065|0x3fd6f6a0433bf384|0.3588|0.359
OUT1-FMT-5G-066|OUT1-FMT-3G-066|0x3fc56c813083d788|0.16737|0.167
OUT1-FMT-5G-067|OUT1-FMT-3G-067|0x3fe51d6d642cd43e|0.65984|0.66
OUT1-FMT-5G-068|OUT1-FMT-3G-068|0x3fe5d698e072be9c|0.68245|0.682
OUT1-FMT-5G-069|OUT1-FMT-3G-069|0xbfd1c5719498a868|-0.27768|-0.278
OUT1-FMT-5G-070|OUT1-FMT-3G-070|0x3fe9a4d32b7cd67e|0.80137|0.801
OUT1-FMT-5G-071|OUT1-FMT-3G-071|0xbfe0086b018afb57|-0.50103|-0.501
OUT1-FMT-5G-072|OUT1-FMT-3G-072|0x3ff7241e9bf2175e|1.4463|1.45
OUT1-FMT-5G-073|OUT1-FMT-3G-073|0x3fe3b9273790255c|0.61635|0.616
OUT1-FMT-5G-074|OUT1-FMT-3G-074|0x3fc2f48bd4cef6b1|0.14809|0.148
OUT1-FMT-5G-075|OUT1-FMT-3G-075|0x3fcc6ed1bf36720a|0.22213|0.222
OUT1-FMT-5G-076|OUT1-FMT-3G-076|0x3fefe6b4ade9f440|0.99691|0.997
OUT1-FMT-5G-077|OUT1-FMT-3G-077|0x3fefe78bae33e0de|0.99701|0.997
OUT1-FMT-5G-078|OUT1-FMT-3G-078|0x3fefe81c9654e8d3|0.99708|0.997
OUT1-FMT-5G-079|OUT1-FMT-3G-079|0xbfcd3e5eea06aa8b|-0.22847|-0.228
OUT1-FMT-5G-080|OUT1-FMT-3G-080|0x3fd162347eb881ed|0.27162|0.272
OUT1-FMT-5G-081|OUT1-FMT-3G-081|0xbffe925acd7c97f4|-1.9107|-1.91
OUT1-FMT-5G-082|OUT1-FMT-3G-082|0x40022c1dc5bc81c6|2.2715|2.27
OUT1-FMT-5G-083|OUT1-FMT-3G-083|0x3facb129986e9821|0.056039|0.056
OUT1-FMT-5G-084|OUT1-FMT-3G-084|0x3f97ab4243d6c9b4|0.023114|0.0231
OUT1-FMT-5G-085|OUT1-FMT-3G-085|0x3fb1c071b2e11747|0.069343|0.0693
OUT1-FMT-5G-086|OUT1-FMT-3G-086|0x3fe5fdd9b033e652|0.68724|0.687
OUT1-FMT-5G-087|OUT1-FMT-3G-087|0x3fe5365ff6571207|0.66289|0.663
OUT1-FMT-5G-088|OUT1-FMT-3G-088|0x3fe4a2f03cf115b9|0.64489|0.645
OUT1-FMT-5G-089|OUT1-FMT-3G-089|0x3fe3e38c85206446|0.62153|0.622
OUT1-FMT-5G-090|OUT1-FMT-3G-090|0x3fce6761237949bc|0.23753|0.238
OUT1-FMT-5G-091|OUT1-FMT-3G-091|0xbfcc800e1dcd02e0|-0.22266|-0.223
OUT1-FMT-5G-092|OUT1-FMT-3G-092|0x3ff781a604e4fc04|1.4692|1.47
OUT1-FMT-5G-093|OUT1-FMT-3G-093|0xbff608d82f6cdecf|-1.3772|-1.38
OUT1-FMT-5G-094|OUT1-FMT-3G-094|0x3fc22638ccd4d489|0.14179|0.142
OUT1-FMT-5G-095|OUT1-FMT-3G-095|0x3fc590321951a958|0.16846|0.168
OUT1-FMT-5G-096|OUT1-FMT-3G-096|0xbfe1c6db7aadf23f|-0.55552|-0.556
OUT1-FMT-5G-097|OUT1-FMT-3G-097|0x3fe0c14e368373bc|0.5236|0.524
OUT1-FMT-5G-098|OUT1-FMT-3G-098|0xbffccb048c5d61d3|-1.7996|-1.8
OUT1-FMT-5G-099|OUT1-FMT-3G-099|0x3ffb254e9225a544|1.6966|1.7
OUT1-FMT-5G-100|OUT1-FMT-3G-100|0x3fb269f693d69d98|0.071929|0.0719
OUT1-FMT-5G-101|OUT1-FMT-3G-101|0x3fb6fb379d622fce|0.089771|0.0898
OUT1-FMT-5G-102|OUT1-FMT-3G-102|0x3fb3e6563fedd3d0|0.077733|0.0777
OUT1-FMT-5G-103|OUT1-FMT-3G-103|0x3ffb5e9669a5ebd5|1.7106|1.71
OUT1-FMT-5G-104|OUT1-FMT-3G-104|0x3fb64fe0a5530be1|0.087156|0.0872
OUT1-FMT-5G-105|OUT1-FMT-3G-105|0x3fc64fe0a5530be1|0.17431|0.174
OUT1-FMT-5G-106|OUT1-FMT-3G-106|0x3f3318b08bd8921e|0.00029139|0.000291
OUT1-FMT-5G-107|OUT1-FMT-3G-107|0x3f21dcf5db27d5ca|0.00013628|0.000136
OUT1-FMT-5G-108|OUT1-FMT-3G-108|0x3fefff1bf4e37424|0.99989|1
"#;
        let mut count = 0_usize;
        for line in ORACLE.lines().filter(|line| !line.is_empty()) {
            let mut fields = line.split('|');
            let id_5 = fields.next().unwrap();
            let id_3 = fields.next().unwrap();
            let bits =
                u64::from_str_radix(fields.next().unwrap().trim_start_matches("0x"), 16).unwrap();
            let expected_5 = fields.next().unwrap();
            let expected_3 = fields.next().unwrap();
            assert!(
                fields.next().is_none(),
                "malformed embedded oracle row: {line}"
            );
            let value = f64::from_bits(bits);
            let observed_5 = format_5g(value).unwrap();
            assert_eq!(observed_5, expected_5, "{id_5}, bits {bits:#018x}, precision .5g, expected {expected_5:?}, observed {observed_5:?}");
            let observed_3 = format_3g(value).unwrap();
            assert_eq!(observed_3, expected_3, "{id_3}, bits {bits:#018x}, precision .3g, expected {expected_3:?}, observed {observed_3:?}");
            count += 1;
        }
        assert_eq!(count, 108, "all 216 stable-ID formatter cases must execute");
    }

    #[test]
    fn fixture_one_outputs_are_byte_exact() {
        let genes = [
            input("GENE_ALPHA", &["a1", "a2"]),
            input("GENE_BETA", &["b1", "b2", "b3"]),
            input("GENE_GAMMA", &["c1", "c2", "c3", "c4"]),
        ];
        let results = [
            result(
                "GENE_ALPHA",
                &["beta_1", "beta_2"],
                vec![0.292434535479278, -0.1686683458689796],
                vec![1.5911201761663327, -0.9176486140225484],
                vec![0.11158253268494193, 0.35880285802687006],
                vec![0.1673737990274129, 0.35880285802687006],
                Some(vec![0.6598422009034335, 0.6824459441904369]),
            ),
            result(
                "GENE_BETA",
                &["beta_1", "beta_2"],
                vec![-0.2776760054717955, 0.8013702248924373],
                vec![-0.5010275869497985, 1.4463182536723802],
                vec![0.6163517079559528, 0.14808795825040513],
                vec![0.6163517079559528, 0.2221319373756077],
                Some(vec![
                    0.9969123265033275,
                    0.9970148470302111,
                    0.9970839439549998,
                ]),
            ),
            result(
                "GENE_GAMMA",
                &["beta_1", "beta_2"],
                vec![-0.22846590449383472, 0.27161896109649214],
                vec![-1.9107311274360628, 2.271541161361003],
                vec![0.05603914247267051, 0.023114238157811065],
                vec![0.1673737990274129, 0.0693427144734332],
                Some(vec![
                    0.6872375909289963,
                    0.6628875552492567,
                    0.6448899450558009,
                    0.6215269661071339,
                ]),
            ),
        ];
        assert_eq!(serialize_gene_summary_to_string(&genes, &results).unwrap(), "Gene\tsgRNA\tbeta_1|beta\tbeta_1|z\tbeta_1|p-value\tbeta_1|fdr\tbeta_1|wald-p-value\tbeta_1|wald-fdr\tbeta_2|beta\tbeta_2|z\tbeta_2|p-value\tbeta_2|fdr\tbeta_2|wald-p-value\tbeta_2|wald-fdr\nGENE_ALPHA\t2\t0.29243\t1.5911\t1.0\t1.0\t0.11158\t0.16737\t-0.16867\t-0.91765\t1.0\t1.0\t0.3588\t0.3588\nGENE_BETA\t3\t-0.27768\t-0.50103\t1.0\t1.0\t0.61635\t0.61635\t0.80137\t1.4463\t1.0\t1.0\t0.14809\t0.22213\nGENE_GAMMA\t4\t-0.22847\t-1.9107\t1.0\t1.0\t0.056039\t0.16737\t0.27162\t2.2715\t1.0\t1.0\t0.023114\t0.069343\n");
        assert_eq!(serialize_sgrna_summary_to_string(&genes, &results).unwrap(), "Gene\tsgRNA\teff\nGENE_ALPHA\ta1\t0.66\nGENE_ALPHA\ta2\t0.682\nGENE_BETA\tb1\t0.997\nGENE_BETA\tb2\t0.997\nGENE_BETA\tb3\t0.997\nGENE_GAMMA\tc1\t0.687\nGENE_GAMMA\tc2\t0.663\nGENE_GAMMA\tc3\t0.645\nGENE_GAMMA\tc4\t0.622\n");
    }

    #[test]
    fn fixture_two_outputs_are_byte_exact() {
        let genes = [
            input("GENE_ALPHA", &["a1", "a2"]),
            input("GENE_BETA", &["b1", "b2", "b3"]),
        ];
        let results = [
            result(
                "GENE_ALPHA",
                &["beta_1", "beta_2"],
                vec![0.2375298903667281, -0.2226579328071443],
                vec![1.4691524687948467, -1.3771592953252598],
                vec![0.14179143907100336, 0.16846300350455112],
                vec![0.14179143907100336, 0.16846300350455112],
                None,
            ),
            result(
                "GENE_BETA",
                &["beta_1", "beta_2"],
                vec![-0.5555245777209804, 0.5235968651501177],
                vec![-1.7995648844717962, 1.6966081341618784],
                vec![0.07192937001547361, 0.08977077095386152],
                vec![0.14179143907100336, 0.16846300350455112],
                None,
            ),
        ];
        assert_eq!(serialize_gene_summary_to_string(&genes, &results).unwrap(), "Gene\tsgRNA\tbeta_1|beta\tbeta_1|z\tbeta_1|p-value\tbeta_1|fdr\tbeta_1|wald-p-value\tbeta_1|wald-fdr\tbeta_2|beta\tbeta_2|z\tbeta_2|p-value\tbeta_2|fdr\tbeta_2|wald-p-value\tbeta_2|wald-fdr\nGENE_ALPHA\t2\t0.23753\t1.4692\t1.0\t1.0\t0.14179\t0.14179\t-0.22266\t-1.3772\t1.0\t1.0\t0.16846\t0.16846\nGENE_BETA\t3\t-0.55552\t-1.7996\t1.0\t1.0\t0.071929\t0.14179\t0.5236\t1.6966\t1.0\t1.0\t0.089771\t0.16846\n");
        assert_eq!(serialize_sgrna_summary_to_string(&genes, &results).unwrap(), "Gene\tsgRNA\teff\nGENE_ALPHA\ta1\t1\nGENE_ALPHA\ta2\t1\nGENE_BETA\tb1\t1\nGENE_BETA\tb2\t1\nGENE_BETA\tb3\t1\n");
    }

    #[test]
    fn fixture_three_retains_nonconverged_gene_byte_exactly() {
        let genes = [
            input("GENE_NORMAL", &["n1", "n2"]),
            input("GENE_EXHAUSTED", &["e1", "e2"]),
        ];
        let mut normal = result(
            "GENE_NORMAL",
            &["beta_1"],
            vec![0.07773341236845899],
            vec![1.7105926634674138],
            vec![0.08715633426895365],
            vec![0.1743126685379073],
            None,
        );
        let mut exhausted = result(
            "GENE_EXHAUSTED",
            &["beta_1"],
            vec![0.00029138860904456994],
            vec![0.0001362848393708218],
            vec![0.9998912604311312],
            vec![0.9998912604311312],
            None,
        );
        exhausted.status = IrlsStatus::IterationLimitReached;
        normal.attempts = 1;
        let results = [normal, exhausted];
        assert!(serialize_gene_summary_to_string(&genes, &results).is_err());
        assert!(serialize_sgrna_summary_to_string(&genes, &results).is_err());
    }

    #[test]
    fn fixture_four_empty_collections_have_exact_headers() {
        assert_eq!(
            serialize_gene_summary_to_string(&[], &[]).unwrap(),
            "Gene\tsgRNA\n"
        );
        assert_eq!(
            serialize_sgrna_summary_to_string(&[], &[]).unwrap(),
            "Gene\tsgRNA\teff\n"
        );
    }

    struct FailWriter {
        calls: usize,
    }
    impl Write for FailWriter {
        fn write(&mut self, _: &[u8]) -> io::Result<usize> {
            self.calls += 1;
            Err(io::Error::other("fail"))
        }
        fn flush(&mut self) -> io::Result<()> {
            Ok(())
        }
    }
    #[test]
    fn validation_precedes_writer_and_writer_errors_propagate() {
        let genes = [input("A", &["a"])];
        let mut results = [result(
            "A",
            &["c"],
            vec![f64::NAN],
            vec![0.0],
            vec![0.5],
            vec![0.5],
            None,
        )];
        let mut writer = FailWriter { calls: 0 };
        assert!(matches!(
            serialize_gene_summary(&genes, &results, &mut writer),
            Err(MageckError::Math(_))
        ));
        assert_eq!(writer.calls, 0);
        results[0].condition_beta[0] = 0.0;
        assert!(matches!(
            serialize_gene_summary(&genes, &results, &mut writer),
            Err(MageckError::Io(_))
        ));
        assert_eq!(writer.calls, 1);
    }

    #[test]
    fn fixture_five_and_six_rejections_are_prewrite_errors() {
        let genes = [input("A", &["a", "b"])];
        let base = result(
            "A",
            &["c"],
            vec![0.0],
            vec![0.0],
            vec![0.5],
            vec![0.5],
            Some(vec![0.8, 0.9]),
        );
        let assert_prewrite = |case: &str, genes: &[SerialGeneInput], results: &[GeneResult]| {
            let mut gene_writer = FailWriter { calls: 0 };
            assert!(
                matches!(
                    serialize_gene_summary(genes, results, &mut gene_writer),
                    Err(MageckError::Math(_))
                ),
                "{case}"
            );
            assert_eq!(gene_writer.calls, 0, "{case}: gene writer was called");
            let mut sgrna_writer = FailWriter { calls: 0 };
            assert!(
                matches!(
                    serialize_sgrna_summary(genes, results, &mut sgrna_writer),
                    Err(MageckError::Math(_))
                ),
                "{case}"
            );
            assert_eq!(sgrna_writer.calls, 0, "{case}: sgRNA writer was called");
        };
        for (case, field, value) in [
            ("5.1", "beta", f64::NAN),
            ("5.2", "beta", f64::INFINITY),
            ("5.3", "beta", f64::NEG_INFINITY),
            ("5.4", "z", f64::NAN),
            ("5.5", "p", f64::NAN),
            ("5.6", "p", 1.00001),
            ("5.7", "p", -0.0001),
            ("5.8", "fdr", f64::NAN),
            ("5.9", "fdr", 1.00001),
            ("5.10", "fdr", -0.0001),
            ("5.11", "eff", f64::NAN),
            ("5.12", "eff", f64::INFINITY),
        ] {
            let mut invalid = base.clone();
            match field {
                "beta" => invalid.condition_beta[0] = value,
                "z" => invalid.z_score[0] = value,
                "p" => invalid.two_sided_p_value[0] = value,
                "fdr" => invalid.condition_fdr[0] = value,
                "eff" => invalid.efficiencies.as_mut().unwrap()[0] = value,
                _ => unreachable!(),
            }
            assert_prewrite(case, &genes, std::slice::from_ref(&invalid));
        }
        let two_genes = [input("A", &["a", "b"]), input("B", &["c", "d"])];
        assert_prewrite("6.1", &two_genes, std::slice::from_ref(&base));
        let mut mismatch_id = base.clone();
        mismatch_id.gene_id = "B".into();
        assert_prewrite("6.2", &genes, std::slice::from_ref(&mismatch_id));
        let mut two_conditions = base.clone();
        two_conditions.condition_order.push("d".into());
        two_conditions.condition_beta.push(0.0);
        two_conditions.z_score.push(0.0);
        two_conditions.two_sided_p_value.push(0.5);
        two_conditions.condition_fdr.push(0.5);
        let mut one_condition = base.clone();
        one_condition.gene_id = "B".into();
        assert_prewrite(
            "6.3",
            &two_genes,
            &[two_conditions.clone(), one_condition.clone()],
        );
        let mut label_mismatch = one_condition.clone();
        label_mismatch.condition_order[0] = "d".into();
        assert_prewrite("6.4", &two_genes, &[base.clone(), label_mismatch]);
        let mut labels_reordered = two_conditions.clone();
        labels_reordered.gene_id = "B".into();
        labels_reordered.condition_order.swap(0, 1);
        assert_prewrite(
            "6.4-reordered",
            &two_genes,
            &[two_conditions, labels_reordered],
        );
        let mut wrong_efficiency_length = base.clone();
        wrong_efficiency_length.efficiencies = Some(vec![0.8]);
        assert_prewrite(
            "6.5",
            &genes,
            std::slice::from_ref(&wrong_efficiency_length),
        );
        for (case, field) in [("6.6", "beta"), ("6.7", "z"), ("6.8", "p"), ("6.9", "fdr")] {
            let mut invalid = base.clone();
            match field {
                "beta" => invalid.condition_beta.clear(),
                "z" => invalid.z_score.clear(),
                "p" => invalid.two_sided_p_value.clear(),
                "fdr" => invalid.condition_fdr.clear(),
                _ => unreachable!(),
            }
            assert_prewrite(case, &genes, std::slice::from_ref(&invalid));
        }
        for (case, id) in [("6.10", "A\tB"), ("6.11", "A\nB"), ("6.12", "A\rB")] {
            let invalid_genes = [input(id, &["a", "b"])];
            assert_prewrite(case, &invalid_genes, std::slice::from_ref(&base));
        }
        let invalid_guides = [input("A", &["a\tb", "b"])];
        assert_prewrite("6.13", &invalid_guides, std::slice::from_ref(&base));
        let mut invalid_label = base;
        invalid_label.condition_order[0] = "c\td".into();
        assert_prewrite("6.14", &genes, std::slice::from_ref(&invalid_label));
    }
}
