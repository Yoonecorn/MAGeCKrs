//! Resolve opaque control identifiers into aligned directional rank pools.

use std::collections::{HashMap, HashSet};

use super::ControlResolutionDiagnostics;
use crate::error::MageckError;

/// Finalized P4F6 inputs plus raw unresolved values for upstream reporting.
#[derive(Debug, Clone, PartialEq)]
#[allow(dead_code)] // Consumed by the later in-crate pipeline.
pub(crate) struct PreparedControlPools {
    pub(crate) neg_percentiles: Vec<f64>,
    pub(crate) pos_percentiles: Vec<f64>,
    pub(crate) diagnostics: ControlResolutionDiagnostics,
    /// Every unresolved raw occurrence, in supplied order (including repeats).
    pub(crate) unresolved_identifiers: Vec<String>,
}

/// Prepare controls from source-row-aligned P4F4 ranks and validity entries.
///
/// Matching is exact and case-sensitive, without trimming or splitting. Matched
/// identifiers are deduplicated before mask filtering, preserving first supplied
/// occurrence order. Unresolved occurrences count individually; duplicated counts
/// only repeated resolved identifiers. The guide universe must have unique IDs.
/// Only retained controls' ranks are extracted and validated. Missing Some ranks,
/// malformed dimensions, ambiguous universe IDs, and an empty final pool are Math
/// errors, as are non-finite or out-of-domain retained percentiles.
#[allow(dead_code)] // CLI and pipeline wiring are deferred.
pub(crate) fn prepare_control_pools(
    raw_identifiers: &[&str],
    guide_identifiers: &[String],
    valid_mask: &[bool],
    neg_ranks: &[Option<f64>],
    pos_ranks: &[Option<f64>],
) -> Result<PreparedControlPools, MageckError> {
    let n = guide_identifiers.len();
    if valid_mask.len() != n || neg_ranks.len() != n || pos_ranks.len() != n {
        return Err(MageckError::Math(format!(
            "Control preparation length mismatch: universe {n}, validity mask {}, negative ranks {}, positive ranks {}",
            valid_mask.len(), neg_ranks.len(), pos_ranks.len()
        )));
    }
    let mut universe = HashMap::with_capacity(n);
    for (index, identifier) in guide_identifiers.iter().enumerate() {
        if let Some(previous) = universe.insert(identifier.as_str(), index) {
            return Err(MageckError::Math(format!(
                "Control preparation ambiguous universe identifier {identifier:?} at rows {previous} and {index}"
            )));
        }
    }
    let mut prepared = PreparedControlPools {
        neg_percentiles: Vec::new(),
        pos_percentiles: Vec::new(),
        diagnostics: ControlResolutionDiagnostics {
            total_supplied: raw_identifiers.len(),
            uniquely_resolved: 0,
            unresolved: 0,
            duplicated: 0,
            validity_excluded: 0,
            final_valid: 0,
        },
        unresolved_identifiers: Vec::new(),
    };
    let mut seen = HashSet::new();
    for &identifier in raw_identifiers {
        let Some(&index) = universe.get(identifier) else {
            prepared.diagnostics.unresolved += 1;
            prepared.unresolved_identifiers.push(identifier.to_owned());
            continue;
        };
        if !seen.insert(identifier) {
            prepared.diagnostics.duplicated += 1;
            continue;
        }
        prepared.diagnostics.uniquely_resolved += 1;
        if !valid_mask[index] {
            prepared.diagnostics.validity_excluded += 1;
            continue;
        }
        for (direction, rank, pool) in [
            ("negative", neg_ranks[index], &mut prepared.neg_percentiles),
            ("positive", pos_ranks[index], &mut prepared.pos_percentiles),
        ] {
            let value = rank.ok_or_else(|| MageckError::Math(format!(
                "Control preparation {direction} percentile missing for valid guide {identifier:?} at row {index}"
            )))?;
            if !value.is_finite() || !(0.0 < value && value <= 1.0) {
                return Err(MageckError::Math(format!(
                    "Control preparation {direction} percentile {value} for guide {identifier:?} at row {index} is non-finite or outside (0, 1]"
                )));
            }
            pool.push(value);
        }
        prepared.diagnostics.final_valid += 1;
    }
    if prepared.diagnostics.final_valid == 0 {
        return Err(MageckError::Math(format!(
            "Control preparation finalized valid pool is empty: supplied {}, unresolved {} {:?}, resolved-but-invalid {}",
            prepared.diagnostics.total_supplied,
            prepared.diagnostics.unresolved,
            prepared.unresolved_identifiers,
            prepared.diagnostics.validity_excluded
        )));
    }
    Ok(prepared)
}

#[cfg(test)]
mod tests {
    use super::*;

    fn ids(values: &[&str]) -> Vec<String> {
        values.iter().map(|value| (*value).into()).collect()
    }

    fn math_error(result: Result<PreparedControlPools, MageckError>) -> String {
        match result {
            Err(MageckError::Math(message)) => message,
            other => panic!("Expected Math error, got {other:?}"),
        }
    }

    #[test]
    fn empty_raw_list() {
        let message = math_error(prepare_control_pools(
            &[],
            &ids(&["a"]),
            &[true],
            &[Some(1.0)],
            &[Some(1.0)],
        ));
        assert!(message.contains("pool is empty"));
        assert!(message.contains("supplied 0"));
    }

    #[test]
    fn all_unresolved() {
        let message = math_error(prepare_control_pools(
            &["unknown", "unknown"],
            &ids(&["a"]),
            &[true],
            &[Some(1.0)],
            &[Some(1.0)],
        ));
        assert!(message.contains("unresolved 2"));
        assert!(message.contains("[\"unknown\", \"unknown\"]"));
    }

    #[test]
    fn all_invalid_deduplicated_before_filtering() {
        let message = math_error(prepare_control_pools(
            &["a", "a", "b"],
            &ids(&["a", "b"]),
            &[false, false],
            &[None, None],
            &[None, None],
        ));
        assert!(message.contains("resolved-but-invalid 2"));
    }

    #[test]
    fn mixed_partial_resolution() {
        let prepared = prepare_control_pools(
            &["missing", "c", "b", "c", "missing", "a", "b"],
            &ids(&["a", "b", "c"]),
            &[true, false, true],
            &[Some(0.25), None, Some(1.0)],
            &[Some(1.0), None, Some(0.5)],
        )
        .unwrap();
        assert_eq!(prepared.neg_percentiles, [1.0, 0.25]);
        assert_eq!(prepared.pos_percentiles, [0.5, 1.0]);
        assert_eq!(prepared.unresolved_identifiers, ["missing", "missing"]);
        assert_eq!(
            prepared.diagnostics,
            ControlResolutionDiagnostics {
                total_supplied: 7,
                uniquely_resolved: 3,
                unresolved: 2,
                duplicated: 2,
                validity_excluded: 1,
                final_valid: 2,
            }
        );
    }

    #[test]
    fn duplicate_deduplication_preserves_first_supplied_order() {
        let prepared = prepare_control_pools(
            &["b", "a", "b", "a"],
            &ids(&["a", "b"]),
            &[true, true],
            &[Some(0.5), Some(1.0)],
            &[Some(1.0), Some(0.5)],
        )
        .unwrap();
        assert_eq!(prepared.neg_percentiles, [1.0, 0.5]);
        assert_eq!(prepared.pos_percentiles, [0.5, 1.0]);
        assert_eq!(prepared.diagnostics.duplicated, 2);
        assert_eq!(prepared.diagnostics.final_valid, 2);
    }

    #[test]
    fn empty_final_pool_after_mixed_filtering() {
        let message = math_error(prepare_control_pools(
            &["unknown", "a"],
            &ids(&["a"]),
            &[false],
            &[None],
            &[None],
        ));
        assert!(message.contains("pool is empty"));
        assert!(message.contains("unresolved 1"));
        assert!(message.contains("resolved-but-invalid 1"));
    }

    #[test]
    fn non_finite_percentiles_in_either_direction() {
        for value in [f64::NAN, f64::INFINITY, f64::NEG_INFINITY] {
            for (neg, pos, direction) in [(value, 1.0, "negative"), (1.0, value, "positive")] {
                let message = math_error(prepare_control_pools(
                    &["a"],
                    &ids(&["a"]),
                    &[true],
                    &[Some(neg)],
                    &[Some(pos)],
                ));
                assert!(message.contains(direction));
                assert!(message.contains("non-finite"));
                assert!(message.contains("row 0"));
            }
        }
    }

    #[test]
    fn out_of_domain_percentiles_in_either_direction() {
        for value in [0.0, -0.0, -0.1, f64::from_bits(1.0_f64.to_bits() + 1)] {
            for (neg, pos) in [(value, 1.0), (1.0, value)] {
                let message = math_error(prepare_control_pools(
                    &["a"],
                    &ids(&["a"]),
                    &[true],
                    &[Some(neg)],
                    &[Some(pos)],
                ));
                assert!(message.contains("outside (0, 1]"));
            }
        }
    }

    #[test]
    fn successful_directional_output_accepts_endpoints_and_equal_values() {
        let tiny = f64::from_bits(1);
        let prepared = prepare_control_pools(
            &["a", "b", "c"],
            &ids(&["a", "b", "c"]),
            &[true; 3],
            &[Some(tiny), Some(1.0), Some(1.0)],
            &[Some(1.0), Some(tiny), Some(tiny)],
        )
        .unwrap();
        assert_eq!(prepared.neg_percentiles, [tiny, 1.0, 1.0]);
        assert_eq!(prepared.pos_percentiles, [1.0, tiny, tiny]);
        assert_eq!(prepared.diagnostics.final_valid, 3);
        assert!(prepared.unresolved_identifiers.is_empty());
    }

    #[test]
    fn opaque_exact_case_sensitive_identifiers() {
        let prepared = prepare_control_pools(
            &["a|b", "a;b", "a,b", "A", " a", "a"],
            &ids(&["a", "a|b", "a;b", "a,b"]),
            &[true; 4],
            &[Some(0.25), Some(0.5), Some(0.75), Some(1.0)],
            &[Some(1.0); 4],
        )
        .unwrap();
        assert_eq!(prepared.neg_percentiles, [0.5, 0.75, 1.0, 0.25]);
        assert_eq!(prepared.unresolved_identifiers, ["A", " a"]);
    }

    #[test]
    fn malformed_dimensions() {
        let universe = ids(&["a"]);
        for (mask, neg, pos) in [
            (&[][..], &[Some(1.0)][..], &[Some(1.0)][..]),
            (&[true][..], &[][..], &[Some(1.0)][..]),
            (&[true][..], &[Some(1.0)][..], &[][..]),
        ] {
            assert!(
                math_error(prepare_control_pools(&["a"], &universe, mask, neg, pos))
                    .contains("length mismatch")
            );
        }
    }

    #[test]
    fn missing_retained_rank_in_either_direction() {
        for (neg, pos) in [(None, Some(1.0)), (Some(1.0), None)] {
            assert!(math_error(prepare_control_pools(
                &["a"],
                &ids(&["a"]),
                &[true],
                &[neg],
                &[pos],
            ))
            .contains("percentile missing"));
        }
    }

    #[test]
    fn malformed_unselected_ranks_are_not_extracted() {
        let prepared = prepare_control_pools(
            &["a", "b"],
            &ids(&["a", "b", "c"]),
            &[true, false, true],
            &[Some(1.0), Some(f64::NAN), None],
            &[Some(0.5), None, Some(-1.0)],
        )
        .unwrap();
        assert_eq!(prepared.neg_percentiles, [1.0]);
        assert_eq!(prepared.diagnostics.validity_excluded, 1);
    }

    #[test]
    fn ambiguous_universe_ids_fail_without_arbitrary_selection() {
        assert!(math_error(prepare_control_pools(
            &["a"],
            &ids(&["a", "a"]),
            &[true, false],
            &[Some(1.0), None],
            &[Some(1.0), None],
        ))
        .contains("ambiguous universe identifier"));
    }

    #[test]
    fn empty_universe_is_empty_pool_error() {
        assert!(
            math_error(prepare_control_pools(&["missing"], &[], &[], &[], &[]))
                .contains("pool is empty")
        );
    }

    #[test]
    fn finalized_output_is_accepted_by_p4f6_control_mode() {
        use super::super::{
            compute_empirical_significance_default, GeneLoValueResult, NullModelConfig,
        };
        let prepared =
            prepare_control_pools(&["a"], &ids(&["a"]), &[true], &[Some(1.0)], &[Some(0.5)])
                .unwrap();
        let observed = [GeneLoValueResult {
            gene_id: "gene".into(),
            num_sgrnas: 1,
            num_valid_sgrnas: 1,
            neg_score: 1.0,
            pos_score: 0.5,
            neg_goodsgrna: 0,
            pos_goodsgrna: 0,
        }];
        let output = compute_empirical_significance_default(
            &observed,
            &NullModelConfig::Control {
                neg_control_percentiles: &prepared.neg_percentiles,
                pos_control_percentiles: &prepared.pos_percentiles,
                diagnostics: &prepared.diagnostics,
            },
            8,
            1,
            crate::stats::DirectionalAlpha {
                negative: 0.25,
                positive: 0.25,
            },
        )
        .unwrap();
        assert_eq!(output[0].neg_pvalue, 0.5);
        assert_eq!(output[0].pos_pvalue, 0.5);
    }
}
