use mageckrs::data::{CountMatrix, CountMatrixError};

#[test]
fn test_new() {
    let mat = CountMatrix::new(3, 4);
    assert_eq!(mat.n_samples(), 3);
    assert_eq!(mat.n_guides(), 4);
    assert_eq!(mat.len(), 12);
    assert!(!mat.is_empty());
    assert_eq!(mat.as_slice(), &[0.0; 12]);
}

#[test]
fn test_from_vec_valid() {
    let data = vec![1.0, 2.0, 3.0, 4.0, 5.0, 6.0];
    let mat = CountMatrix::from_vec(data.clone(), 2, 3).expect("Valid dimensions");
    assert_eq!(mat.n_samples(), 2);
    assert_eq!(mat.n_guides(), 3);
    assert_eq!(mat.len(), 6);
    assert_eq!(mat.as_slice(), &data[..]);
}

#[test]
fn test_from_slice_valid() {
    let raw = [10.0, 20.0, 30.0, 40.0];
    let mat = CountMatrix::from_slice(&raw, 2, 2).expect("Valid slice dimensions");
    assert_eq!(mat.n_samples(), 2);
    assert_eq!(mat.n_guides(), 2);
    assert_eq!(mat.len(), 4);
    assert_eq!(mat.as_slice(), &raw);
}

#[test]
fn test_rejection_incorrect_data_length() {
    let data = vec![1.0, 2.0, 3.0, 4.0, 5.0]; // 5 elements
    let err = CountMatrix::from_vec(data, 2, 3).unwrap_err();
    assert_eq!(
        err,
        CountMatrixError::DimensionMismatch {
            expected: 6,
            actual: 5,
        }
    );
}

#[test]
fn test_sample_major_index_layout() {
    // 2 samples, 3 guides:
    // Sample 0: guide 0 = 10, guide 1 = 11, guide 2 = 12
    // Sample 1: guide 0 = 20, guide 1 = 21, guide 2 = 22
    let data = vec![10.0, 11.0, 12.0, 20.0, 21.0, 22.0];
    let mat = CountMatrix::from_vec(data, 2, 3).unwrap();

    // Verify index layout: s * n_guides + g
    assert_eq!(mat.get(0, 0), Some(10.0));
    assert_eq!(mat.get(0, 1), Some(11.0));
    assert_eq!(mat.get(0, 2), Some(12.0));
    assert_eq!(mat.get(1, 0), Some(20.0));
    assert_eq!(mat.get(1, 1), Some(21.0));
    assert_eq!(mat.get(1, 2), Some(22.0));

    // Entry & Index operator
    assert_eq!(mat.entry(0, 0), 10.0);
    assert_eq!(mat.entry(1, 2), 22.0);
    assert_eq!(mat[(0, 1)], 11.0);
    assert_eq!(mat[(1, 0)], 20.0);
}

#[test]
fn test_first_and_last_valid_elements() {
    let n_samples = 4;
    let n_guides = 5;
    let mut data = vec![0.0; n_samples * n_guides];
    data[0] = 99.0; // first element
    data[n_samples * n_guides - 1] = 88.0; // last element

    let mat = CountMatrix::from_vec(data, n_samples, n_guides).unwrap();
    assert_eq!(mat.get(0, 0), Some(99.0));
    assert_eq!(mat.get(n_samples - 1, n_guides - 1), Some(88.0));
    assert_eq!(mat[(0, 0)], 99.0);
    assert_eq!(mat[(n_samples - 1, n_guides - 1)], 88.0);
}

#[test]
fn test_sample_slice_boundaries() {
    let data = vec![
        1.0, 2.0, 3.0, // sample 0
        4.0, 5.0, 6.0, // sample 1
        7.0, 8.0, 9.0, // sample 2
    ];
    let mat = CountMatrix::from_vec(data, 3, 3).unwrap();

    assert_eq!(mat.sample_slice(0), &[1.0, 2.0, 3.0]);
    assert_eq!(mat.sample_slice(1), &[4.0, 5.0, 6.0]);
    assert_eq!(mat.sample_slice(2), &[7.0, 8.0, 9.0]);

    assert_eq!(mat.get_sample_slice(0), Some(&[1.0, 2.0, 3.0][..]));
    assert_eq!(mat.get_sample_slice(1), Some(&[4.0, 5.0, 6.0][..]));
    assert_eq!(mat.get_sample_slice(2), Some(&[7.0, 8.0, 9.0][..]));
    assert_eq!(mat.get_sample_slice(3), None);
}

#[test]
fn test_sample_slice_mut_and_updates() {
    let mut mat = CountMatrix::new(2, 3);
    {
        let s0 = mat.sample_slice_mut(0);
        s0[0] = 100.0;
        s0[1] = 101.0;
        s0[2] = 102.0;
    }
    {
        let s1 = mat.sample_slice_mut(1);
        s1[0] = 200.0;
        s1[1] = 201.0;
        s1[2] = 202.0;
    }

    assert_eq!(mat.sample_slice(0), &[100.0, 101.0, 102.0]);
    assert_eq!(mat.sample_slice(1), &[200.0, 201.0, 202.0]);

    // Test set()
    mat.set(0, 1, 999.0);
    assert_eq!(mat.get(0, 1), Some(999.0));

    // Test IndexMut
    mat[(1, 2)] = 777.0;
    assert_eq!(mat[(1, 2)], 777.0);

    // Test get_mut()
    if let Some(val) = mat.get_mut(1, 0) {
        *val = 333.0;
    }
    assert_eq!(mat.get(1, 0), Some(333.0));
}

#[test]
fn test_independence_between_sample_slices() {
    let mut mat = CountMatrix::new(3, 3);
    mat.sample_slice_mut(0)[1] = 42.0;

    // Sample 1 and 2 must remain unmodified
    assert_eq!(mat.sample_slice(0), &[0.0, 42.0, 0.0]);
    assert_eq!(mat.sample_slice(1), &[0.0, 0.0, 0.0]);
    assert_eq!(mat.sample_slice(2), &[0.0, 0.0, 0.0]);
}

#[test]
fn test_guide_iter() {
    // 3 samples, 2 guides:
    // sample 0: [1.0, 2.0]
    // sample 1: [3.0, 4.0]
    // sample 2: [5.0, 6.0]
    let data = vec![1.0, 2.0, 3.0, 4.0, 5.0, 6.0];
    let mat = CountMatrix::from_vec(data, 3, 2).unwrap();

    let guide0: Vec<f64> = mat.guide_iter(0).copied().collect();
    assert_eq!(guide0, vec![1.0, 3.0, 5.0]);

    let guide1: Vec<f64> = mat.guide_iter(1).copied().collect();
    assert_eq!(guide1, vec![2.0, 4.0, 6.0]);

    // Test DoubleEndedIterator and ExactSizeIterator
    let mut iter0 = mat.guide_iter(0);
    assert_eq!(iter0.len(), 3);
    assert_eq!(iter0.next(), Some(&1.0));
    assert_eq!(iter0.len(), 2);
    assert_eq!(iter0.next_back(), Some(&5.0));
    assert_eq!(iter0.len(), 1);
    assert_eq!(iter0.next(), Some(&3.0));
    assert_eq!(iter0.len(), 0);
    assert_eq!(iter0.next(), None);
    assert_eq!(iter0.next_back(), None);
}

#[test]
fn test_as_ptr_and_as_slice() {
    let data = vec![1.5, 2.5, 3.5, 4.5];
    let mut mat = CountMatrix::from_vec(data.clone(), 2, 2).unwrap();

    assert_eq!(mat.as_slice(), &data[..]);
    assert_eq!(mat.as_ptr(), mat.as_slice().as_ptr());

    let mut_slice = mat.as_mut_slice();
    mut_slice[0] = 9.5;
    assert_eq!(mat.sample_slice(0)[0], 9.5);

    let consumed = mat.into_vec();
    assert_eq!(consumed, vec![9.5, 2.5, 3.5, 4.5]);
}

#[test]
fn test_zero_dimensions() {
    // Zero samples, positive guides
    let m0 = CountMatrix::new(0, 5);
    assert_eq!(m0.n_samples(), 0);
    assert_eq!(m0.n_guides(), 5);
    assert_eq!(m0.len(), 0);
    assert!(m0.is_empty());
    assert_eq!(m0.get_sample_slice(0), None);
    assert_eq!(m0.get(0, 0), None);
    let mut iter = m0.guide_iter(0);
    assert_eq!(iter.len(), 0);
    assert_eq!(iter.next(), None);

    // Positive samples, zero guides
    let m1 = CountMatrix::new(5, 0);
    assert_eq!(m1.n_samples(), 5);
    assert_eq!(m1.n_guides(), 0);
    assert_eq!(m1.len(), 0);
    assert!(m1.is_empty());
    assert_eq!(m1.sample_slice(0), &[]);
    assert_eq!(m1.sample_slice(4), &[]);
    assert_eq!(m1.get_sample_slice(5), None);
    assert_eq!(m1.get(0, 0), None);
    assert_eq!(m1.get_guide_iter(0), None);

    // Zero samples, zero guides
    let m2 = CountMatrix::new(0, 0);
    assert_eq!(m2.n_samples(), 0);
    assert_eq!(m2.n_guides(), 0);
    assert_eq!(m2.len(), 0);
    assert!(m2.is_empty());

    // from_vec with empty data and zero dimensions
    let m3 = CountMatrix::from_vec(vec![], 0, 0).unwrap();
    assert_eq!(m3, m2);
}

#[test]
fn test_overflow_safe_dimension_validation() {
    // Attempt dimensions that overflow usize
    let huge_dim = usize::MAX;
    let err = CountMatrix::from_vec(vec![], huge_dim, 2).unwrap_err();
    assert_eq!(
        err,
        CountMatrixError::DimensionOverflow {
            n_samples: huge_dim,
            n_guides: 2,
        }
    );
}

#[test]
#[should_panic(expected = "sample index 2 out of bounds (n_samples = 2)")]
fn test_sample_slice_out_of_bounds_panics() {
    let mat = CountMatrix::new(2, 3);
    let _ = mat.sample_slice(2);
}

#[test]
#[should_panic(expected = "sample index 2 out of bounds (n_samples = 2)")]
fn test_sample_slice_mut_out_of_bounds_panics() {
    let mut mat = CountMatrix::new(2, 3);
    let _ = mat.sample_slice_mut(2);
}

#[test]
#[should_panic(expected = "guide index 3 out of bounds (n_guides = 3)")]
fn test_guide_iter_out_of_bounds_panics() {
    let mat = CountMatrix::new(2, 3);
    let _ = mat.guide_iter(3);
}

#[test]
#[should_panic(expected = "guide index 0 out of bounds (n_guides = 0)")]
fn test_guide_iter_zero_guides_panics() {
    let mat = CountMatrix::new(2, 0);
    let _ = mat.guide_iter(0);
}

#[test]
#[should_panic(expected = "sample index 5 out of bounds (n_samples = 2)")]
fn test_index_out_of_bounds_panics() {
    let mat = CountMatrix::new(2, 3);
    let _ = mat[(5, 0)];
}

#[test]
#[should_panic(expected = "guide index 5 out of bounds (n_guides = 3)")]
fn test_index_guide_out_of_bounds_panics() {
    let mat = CountMatrix::new(2, 3);
    let _ = mat[(0, 5)];
}

// ── StridedIter rigorous tests ─────────────────────────────────────────────

/// Pure reverse iteration (next_back only) must yield samples in reverse order.
#[test]
fn test_guide_iter_reverse_only() {
    // 4 samples, 3 guides; guide 1 values = [10, 20, 30, 40] across samples
    let data = vec![
        0.0, 10.0, 0.0, // sample 0
        0.0, 20.0, 0.0, // sample 1
        0.0, 30.0, 0.0, // sample 2
        0.0, 40.0, 0.0, // sample 3
    ];
    let mat = CountMatrix::from_vec(data, 4, 3).unwrap();
    let collected: Vec<f64> = mat.guide_iter(1).copied().rev().collect();
    assert_eq!(collected, vec![40.0, 30.0, 20.0, 10.0]);
}

/// Interleaved next / next_back must correctly track remaining and not OOB.
#[test]
fn test_guide_iter_interleaved_partial_consumption() {
    let data = vec![
        1.0, 0.0, // sample 0: guide0=1
        2.0, 0.0, // sample 1: guide0=2
        3.0, 0.0, // sample 2: guide0=3
        4.0, 0.0, // sample 3: guide0=4
        5.0, 0.0, // sample 4: guide0=5
    ];
    let mat = CountMatrix::from_vec(data, 5, 2).unwrap();
    let mut it = mat.guide_iter(0);

    assert_eq!(it.size_hint(), (5, Some(5)));
    assert_eq!(it.len(), 5);

    assert_eq!(it.next(), Some(&1.0));
    assert_eq!(it.len(), 4);
    assert_eq!(it.size_hint(), (4, Some(4)));

    assert_eq!(it.next_back(), Some(&5.0));
    assert_eq!(it.len(), 3);

    assert_eq!(it.next_back(), Some(&4.0));
    assert_eq!(it.len(), 2);

    assert_eq!(it.next(), Some(&2.0));
    assert_eq!(it.len(), 1);

    // Only element 3 remains
    assert_eq!(it.next(), Some(&3.0));
    assert_eq!(it.len(), 0);
    assert_eq!(it.size_hint(), (0, Some(0)));

    // Exhausted — both directions must return None (Fused behavior)
    assert_eq!(it.next(), None);
    assert_eq!(it.next_back(), None);
    assert_eq!(it.next(), None);
}

/// ExactSizeIterator: len() matches remaining count at every step.
#[test]
fn test_guide_iter_exact_size_at_every_step() {
    let n_samples = 6_usize;
    let n_guides = 4_usize;
    let data: Vec<f64> = (0..n_samples * n_guides).map(|x| x as f64).collect();
    let mat = CountMatrix::from_vec(data, n_samples, n_guides).unwrap();

    let mut it = mat.guide_iter(2);
    for expected_remaining in (0..=n_samples).rev() {
        assert_eq!(it.len(), expected_remaining);
        assert_eq!(
            it.size_hint(),
            (expected_remaining, Some(expected_remaining))
        );
        if expected_remaining > 0 {
            let _ = it.next();
        }
    }
}

/// Zero-sample iterator: first call to next() and next_back() must both return None.
#[test]
fn test_guide_iter_zero_samples_fused() {
    let mat = CountMatrix::new(0, 4);
    let mut it = mat.guide_iter(0);
    assert_eq!(it.len(), 0);
    assert_eq!(it.size_hint(), (0, Some(0)));
    assert_eq!(it.next(), None);
    assert_eq!(it.next_back(), None);
    // Call again to confirm FusedIterator contract
    assert_eq!(it.next(), None);
    assert_eq!(it.next_back(), None);
}

/// One-sample matrix: forward and backward yield the same single element.
#[test]
fn test_guide_iter_single_sample() {
    let data = vec![7.0, 8.0, 9.0]; // 1 sample, 3 guides
    let mat = CountMatrix::from_vec(data, 1, 3).unwrap();

    // Forward
    let fwd: Vec<f64> = mat.guide_iter(1).copied().collect();
    assert_eq!(fwd, vec![8.0]);

    // Backward
    let bwd: Vec<f64> = mat.guide_iter(1).copied().rev().collect();
    assert_eq!(bwd, vec![8.0]);

    // Interleaved exhausts immediately
    let mut it = mat.guide_iter(2);
    assert_eq!(it.len(), 1);
    assert_eq!(it.next_back(), Some(&9.0));
    assert_eq!(it.len(), 0);
    assert_eq!(it.next(), None);
}

/// First guide (g=0) and last guide (g=n_guides-1) are both reachable.
#[test]
fn test_guide_iter_first_and_last_guide() {
    let data = vec![
        10.0, 20.0, 30.0, // sample 0
        40.0, 50.0, 60.0, // sample 1
        70.0, 80.0, 90.0, // sample 2
    ];
    let mat = CountMatrix::from_vec(data, 3, 3).unwrap();

    let first: Vec<f64> = mat.guide_iter(0).copied().collect();
    assert_eq!(first, vec![10.0, 40.0, 70.0]);

    let last: Vec<f64> = mat.guide_iter(2).copied().collect();
    assert_eq!(last, vec![30.0, 60.0, 90.0]);

    let first_rev: Vec<f64> = mat.guide_iter(0).copied().rev().collect();
    assert_eq!(first_rev, vec![70.0, 40.0, 10.0]);

    let last_rev: Vec<f64> = mat.guide_iter(2).copied().rev().collect();
    assert_eq!(last_rev, vec![90.0, 60.0, 30.0]);
}

/// Collecting forward == reverse of collecting backward for any guide.
#[test]
fn test_guide_iter_forward_reverse_symmetry() {
    let n_samples = 5;
    let n_guides = 7;
    let data: Vec<f64> = (0..n_samples * n_guides).map(|x| x as f64 * 1.5).collect();
    let mat = CountMatrix::from_vec(data, n_samples, n_guides).unwrap();

    for g in 0..n_guides {
        let fwd: Vec<f64> = mat.guide_iter(g).copied().collect();
        let mut rev: Vec<f64> = mat.guide_iter(g).copied().rev().collect();
        rev.reverse();
        assert_eq!(fwd, rev, "guide {g}: forward/backward mismatch");
    }
}
