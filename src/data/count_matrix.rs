use std::fmt;

/// Errors that can occur when constructing or manipulating a [`CountMatrix`].
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum CountMatrixError {
    /// The provided data length does not match `n_samples * n_guides`.
    DimensionMismatch {
        /// Expected element count (`n_samples * n_guides`).
        expected: usize,
        /// Actual slice or vector length provided.
        actual: usize,
    },
    /// The product `n_samples * n_guides` overflows `usize`.
    DimensionOverflow {
        /// Number of samples requested.
        n_samples: usize,
        /// Number of guides requested.
        n_guides: usize,
    },
}

impl fmt::Display for CountMatrixError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::DimensionMismatch { expected, actual } => {
                write!(
                    f,
                    "CountMatrix dimension mismatch: expected data length {expected} (n_samples * n_guides), got {actual}"
                )
            }
            Self::DimensionOverflow {
                n_samples,
                n_guides,
            } => {
                write!(
                    f,
                    "CountMatrix dimension overflow: n_samples ({n_samples}) * n_guides ({n_guides}) overflows usize"
                )
            }
        }
    }
}

impl std::error::Error for CountMatrixError {}

/// Contiguous, sample-major 2D numerical matrix storing count data.
///
/// Data layout is contiguous row-major where each row represents one sample:
/// `index = sample_idx * n_guides + guide_idx`.
///
/// This contiguous layout enables zero-copy passing to Apple Accelerate and
/// Metal GPU shared memory buffers.
#[derive(Debug, Clone, PartialEq)]
pub struct CountMatrix {
    data: Vec<f64>,
    n_samples: usize,
    n_guides: usize,
}

impl CountMatrix {
    /// Creates a new zero-initialized [`CountMatrix`] with dimensions `n_samples` x `n_guides`.
    ///
    /// # Panics
    /// Panics if `n_samples * n_guides` overflows `usize`.
    pub fn new(n_samples: usize, n_guides: usize) -> Self {
        let total = n_samples
            .checked_mul(n_guides)
            .expect("CountMatrix dimension overflow: n_samples * n_guides exceeds usize::MAX");
        Self {
            data: vec![0.0; total],
            n_samples,
            n_guides,
        }
    }

    /// Creates a [`CountMatrix`] from existing contiguous vector data and dimensions.
    ///
    /// # Errors
    /// Returns [`CountMatrixError::DimensionOverflow`] if `n_samples * n_guides` overflows `usize`.
    /// Returns [`CountMatrixError::DimensionMismatch`] if `data.len() != n_samples * n_guides`.
    pub fn from_vec(
        data: Vec<f64>,
        n_samples: usize,
        n_guides: usize,
    ) -> Result<Self, CountMatrixError> {
        let expected =
            n_samples
                .checked_mul(n_guides)
                .ok_or(CountMatrixError::DimensionOverflow {
                    n_samples,
                    n_guides,
                })?;

        if data.len() != expected {
            return Err(CountMatrixError::DimensionMismatch {
                expected,
                actual: data.len(),
            });
        }

        Ok(Self {
            data,
            n_samples,
            n_guides,
        })
    }

    /// Creates a [`CountMatrix`] by copying data from a slice.
    ///
    /// # Errors
    /// Returns [`CountMatrixError::DimensionOverflow`] if `n_samples * n_guides` overflows `usize`.
    /// Returns [`CountMatrixError::DimensionMismatch`] if `slice.len() != n_samples * n_guides`.
    pub fn from_slice(
        slice: &[f64],
        n_samples: usize,
        n_guides: usize,
    ) -> Result<Self, CountMatrixError> {
        Self::from_vec(slice.to_vec(), n_samples, n_guides)
    }

    /// Returns the number of samples (rows).
    #[inline]
    pub fn n_samples(&self) -> usize {
        self.n_samples
    }

    /// Returns the number of guides (columns).
    #[inline]
    pub fn n_guides(&self) -> usize {
        self.n_guides
    }

    /// Returns the total number of elements (`n_samples * n_guides`).
    #[inline]
    pub fn len(&self) -> usize {
        self.data.len()
    }

    /// Returns `true` if the matrix contains 0 elements.
    #[inline]
    pub fn is_empty(&self) -> bool {
        self.data.is_empty()
    }

    /// Returns a raw pointer to the underlying contiguous data buffer.
    ///
    /// Suitable for passing directly to Apple Accelerate routines or Metal shared buffers.
    #[inline]
    pub fn as_ptr(&self) -> *const f64 {
        self.data.as_ptr()
    }

    /// Returns a reference to the entire underlying data as a flat slice.
    #[inline]
    pub fn as_slice(&self) -> &[f64] {
        &self.data
    }

    /// Returns a mutable reference to the entire underlying data as a flat slice.
    #[inline]
    pub fn as_mut_slice(&mut self) -> &mut [f64] {
        &mut self.data
    }

    /// Consumes the matrix, returning the underlying `Vec<f64>`.
    #[inline]
    pub fn into_vec(self) -> Vec<f64> {
        self.data
    }

    /// Returns a contiguous slice of counts for sample `s`.
    ///
    /// # Panics
    /// Panics if `s >= n_samples`.
    #[inline]
    pub fn sample_slice(&self, s: usize) -> &[f64] {
        assert!(
            s < self.n_samples,
            "sample index {s} out of bounds (n_samples = {})",
            self.n_samples
        );
        let start = s * self.n_guides;
        let end = start + self.n_guides;
        &self.data[start..end]
    }

    /// Returns a mutable contiguous slice of counts for sample `s`.
    ///
    /// # Panics
    /// Panics if `s >= n_samples`.
    #[inline]
    pub fn sample_slice_mut(&mut self, s: usize) -> &mut [f64] {
        assert!(
            s < self.n_samples,
            "sample index {s} out of bounds (n_samples = {})",
            self.n_samples
        );
        let start = s * self.n_guides;
        let end = start + self.n_guides;
        &mut self.data[start..end]
    }

    /// Returns a contiguous slice of counts for sample `s`, or `None` if out of bounds.
    #[inline]
    pub fn get_sample_slice(&self, s: usize) -> Option<&[f64]> {
        if s < self.n_samples {
            let start = s * self.n_guides;
            let end = start + self.n_guides;
            Some(&self.data[start..end])
        } else {
            None
        }
    }

    /// Returns a mutable contiguous slice of counts for sample `s`, or `None` if out of bounds.
    #[inline]
    pub fn get_sample_slice_mut(&mut self, s: usize) -> Option<&mut [f64]> {
        if s < self.n_samples {
            let start = s * self.n_guides;
            let end = start + self.n_guides;
            Some(&mut self.data[start..end])
        } else {
            None
        }
    }

    /// Returns a strided iterator over all samples for guide `g`.
    ///
    /// Stride equals `n_guides`.
    ///
    /// # Panics
    /// Panics if `g >= n_guides`.
    #[inline]
    pub fn guide_iter(&self, g: usize) -> StridedIter<'_, f64> {
        assert!(
            g < self.n_guides,
            "guide index {g} out of bounds (n_guides = {})",
            self.n_guides
        );
        StridedIter {
            slice: &self.data,
            offset: g,
            stride: self.n_guides,
            remaining: self.n_samples,
        }
    }

    /// Returns a strided iterator over all samples for guide `g`, or `None` if `g >= n_guides`.
    #[inline]
    pub fn get_guide_iter(&self, g: usize) -> Option<StridedIter<'_, f64>> {
        if g < self.n_guides {
            Some(StridedIter {
                slice: &self.data,
                offset: g,
                stride: self.n_guides,
                remaining: self.n_samples,
            })
        } else {
            None
        }
    }

    /// Returns the element at sample `s` and guide `g`, or `None` if out of bounds.
    #[inline]
    pub fn get(&self, s: usize, g: usize) -> Option<f64> {
        if s < self.n_samples && g < self.n_guides {
            Some(self.data[s * self.n_guides + g])
        } else {
            None
        }
    }

    /// Returns a mutable reference to the element at sample `s` and guide `g`, or `None` if out of bounds.
    #[inline]
    pub fn get_mut(&mut self, s: usize, g: usize) -> Option<&mut f64> {
        if s < self.n_samples && g < self.n_guides {
            Some(&mut self.data[s * self.n_guides + g])
        } else {
            None
        }
    }

    /// Returns the element at sample `s` and guide `g`.
    ///
    /// # Panics
    /// Panics if `s >= n_samples` or `g >= n_guides`.
    #[inline]
    pub fn entry(&self, s: usize, g: usize) -> f64 {
        assert!(
            s < self.n_samples,
            "sample index {s} out of bounds (n_samples = {})",
            self.n_samples
        );
        assert!(
            g < self.n_guides,
            "guide index {g} out of bounds (n_guides = {})",
            self.n_guides
        );
        self.data[s * self.n_guides + g]
    }

    /// Sets the element at sample `s` and guide `g` to `value`.
    ///
    /// # Panics
    /// Panics if `s >= n_samples` or `g >= n_guides`.
    #[inline]
    pub fn set(&mut self, s: usize, g: usize, value: f64) {
        assert!(
            s < self.n_samples,
            "sample index {s} out of bounds (n_samples = {})",
            self.n_samples
        );
        assert!(
            g < self.n_guides,
            "guide index {g} out of bounds (n_guides = {})",
            self.n_guides
        );
        self.data[s * self.n_guides + g] = value;
    }
}

impl std::ops::Index<(usize, usize)> for CountMatrix {
    type Output = f64;

    #[inline]
    fn index(&self, (s, g): (usize, usize)) -> &Self::Output {
        assert!(
            s < self.n_samples,
            "sample index {s} out of bounds (n_samples = {})",
            self.n_samples
        );
        assert!(
            g < self.n_guides,
            "guide index {g} out of bounds (n_guides = {})",
            self.n_guides
        );
        &self.data[s * self.n_guides + g]
    }
}

impl std::ops::IndexMut<(usize, usize)> for CountMatrix {
    #[inline]
    fn index_mut(&mut self, (s, g): (usize, usize)) -> &mut Self::Output {
        assert!(
            s < self.n_samples,
            "sample index {s} out of bounds (n_samples = {})",
            self.n_samples
        );
        assert!(
            g < self.n_guides,
            "guide index {g} out of bounds (n_guides = {})",
            self.n_guides
        );
        &mut self.data[s * self.n_guides + g]
    }
}

/// A strided iterator over a slice of items.
///
/// Yields items with a fixed stride step between successive indices.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct StridedIter<'a, T = f64> {
    slice: &'a [T],
    offset: usize,
    stride: usize,
    remaining: usize,
}

impl<'a, T> Iterator for StridedIter<'a, T> {
    type Item = &'a T;

    #[inline]
    fn next(&mut self) -> Option<Self::Item> {
        if self.remaining == 0 {
            None
        } else {
            let item = &self.slice[self.offset];
            self.offset += self.stride;
            self.remaining -= 1;
            Some(item)
        }
    }

    #[inline]
    fn size_hint(&self) -> (usize, Option<usize>) {
        (self.remaining, Some(self.remaining))
    }
}

impl<'a, T> ExactSizeIterator for StridedIter<'a, T> {
    #[inline]
    fn len(&self) -> usize {
        self.remaining
    }
}

impl<'a, T> std::iter::FusedIterator for StridedIter<'a, T> {}

impl<'a, T> DoubleEndedIterator for StridedIter<'a, T> {
    #[inline]
    fn next_back(&mut self) -> Option<Self::Item> {
        if self.remaining == 0 {
            None
        } else {
            self.remaining -= 1;
            let item = &self.slice[self.offset + self.remaining * self.stride];
            Some(item)
        }
    }
}
