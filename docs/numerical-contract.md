# MLE numerical contract v2

MAGeCKrs fits a negative-binomial model with guide intercepts and condition coefficients. This document describes the numerical acceptance rules used by the implementation. Its SHA-256 is recorded in MLE diagnostic sidecars. Convergence is a numerical property and does not establish biological accuracy or a global optimum.

All guide and condition coefficients must be finite and have absolute value at most 100. For each update, the maximum coefficient change scaled by `1e-8 + 1e-6 * max(abs(new), abs(old))` must be at most one. At the returned vector, the scaled IRLS fixed-point residual must be at most `1e-9`. Both criteria must hold for three consecutive accepted updates. A failed criterion resets the streak.

The fit uses constant relaxation 1.0, then eligible retries use 0.5 and 0.25. Each attempt allows at most 1,000 candidate updates and restarts from the original initialization. Dispersion and ridge parameters remain unchanged during retries; under default scalar MLE each is 0.01. Retrying changes the numerical trajectory of the same model. Invalid initial states and inference-only failures do not trigger futile repeated solves.

Ordinary Wald inference is rebuilt at the returned coefficient vector using the existing scaled sandwich covariance. Standard errors must be finite and strictly positive. Invalid fit or inference states produce diagnostics rather than verified primary summaries. Primary output requires every otherwise planned targeting hypothesis to be eligible for all requested conditions.

Wald Benjamini-Hochberg correction uses unrounded p-values over the emitted, eligible targeting-gene universe, separately per condition. Control, predeclared skipped genes and invalid hypotheses are excluded. Empty eligible universes fail; a singleton retains q=p. Original gene, guide and condition order is preserved.

Diagnostic states distinguish primary convergence, convergence after a same-model retry and failed attempts. MLE permutation columns remain 1.0 placeholders. A CNV adapter changes matched betas after ordinary MLE inference without recomputing uncertainty or significance.
