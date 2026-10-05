## Version 1.1.0

This is a resubmission of version 1.1.0, updating the CRAN release 1.0.2.

## R CMD check results

0 errors | 0 warnings | 0 notes

* `devtools::check(remote = TRUE, manual = TRUE)` passed on macOS arm64 with
  R 4.6.0.
* `urlchecker::url_check()` reported that all URLs are correct.
* `revdepcheck::revdep_check(num_workers = 4)` found no reverse dependencies.
* Win-Builder R-devel Windows check (R Under development 2026-09-30,
  r90605 ucrt) passed with status OK.

This release adds the strong k-out-of-m adjustment, improves update and
plotting methods, and corrects negative-binomial dispersion scaling for
aggregate counts in parallel designs.
