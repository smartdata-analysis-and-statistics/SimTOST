# SimTOST 1.1.0

## New features

- The `adjust = "t"` option now supports Mielke's strong k-out-of-m
  adjustment.
- Diagnostic and sample-size plots now support endpoint selection.
- Count-outcome workflows now support Poisson and negative-binomial
  sensitivity analyses, with expanded documentation.

## Minor improvements and fixes

- `update()` now handles comparator, endpoint, mean, standard-deviation, and
  covariance information more consistently.
- Negative-binomial dispersion scaling for aggregate counts in parallel
  designs is now correct.
- Monte Carlo stability and precision diagnostics are improved.

# SimTOST 1.0.2

## Bug fixes

- A negative value (`-1`) is now handled correctly in relevant C++ code,
  avoiding an `unsigned int` conversion error.

# SimTOST 1.0.1
- Fixed CRAN review issues: expanded description, added references, documented function outputs.

# SimTOST 1.0.0
* Initial CRAN submission.

# SimTOST 0.6.0

# SimTOST 0.5.0

# SimTOST 0.4.0

# SimTOST 0.3.0

# SimTOST 0.2.0

# SimTOST 0.1.0
