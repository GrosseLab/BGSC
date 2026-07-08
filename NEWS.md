# BGSC 0.0.1.4

This release fixes correctness and robustness issues found during a code and
mathematical review of the package (see `CODE_REVIEW.md` for the full write-up
of what was checked, what was fixed, and what was intentionally left as-is).

## Bug fixes

* `make.plot.data.exp.Ill.qPCR()` was missing `qgenesIDs` as a formal argument
  and instead relied on it existing as a variable in the caller's global
  environment. Because this is a package function, its lexical scope is the
  package namespace, not the calling script's environment, so the call in
  `ReproducibleScript.Rmd` / `scripts/MAIN.R` would fail with
  `object 'qgenesIDs' not found`. `qgenesIDs` is now an explicit parameter,
  and both call sites pass it explicitly.
* `union5()` called `Reduce(intersect, ...)` internally instead of
  `Reduce(union, ...)` (a copy-paste leftover from `intersect5()`), so it
  silently returned the intersection of its five sets instead of their union.
* `comparativeMethod_qPCR.analysis()`: the `ddCTPower2`/`ddCTRelExp` columns
  used `2^(+ddCt)` instead of the standard `2^(-ddCt)` fold-change convention,
  i.e. they returned the reciprocal of the intended relative expression. Fixed
  to match the (correct) `Potenz`/`relative_values` columns computed earlier
  in the same function. These columns are not consumed anywhere else in the
  package, so this did not affect any previously published result.
* `comparativeMethod_qPCR.RNAi.log2FC()`: the `dCT.C1C0` and `relVal.C1C0`
  result data frames assigned the same column name twice (e.g.
  `"dCT.mean.C1.stderr"` used for both the C1 and the C0 standard error),
  so R silently renamed the second occurrence. The C0 columns are now named
  correctly (`dCT.mean.C0.stderr`, `relVal.geomean.C0.stderr`). The actually
  reported `dCT.logFC.C1C0.stderr` / `relVal.geomean.logFC.C1C0.stderr`
  values were computed from the correct underlying variables before this
  collision and are unchanged.
* `InformationCriterion(..., IC = "AIC")` used the `k` argument (number of
  observations) as the AIC penalty multiplier. Standard AIC uses a fixed
  penalty of 2, independent of sample size; `k` was only ever meant for BIC's
  `log(k)` term. Fixed so `AIC = -2*logL + 2*npar`. The `BIC` branch (the only
  one actually used in `ReproducibleScript.Rmd`/`scripts/MAIN.R`) is
  numerically unchanged.
* `InformationCriterion()` now raises an error (`stop()`) for mismatched
  argument lengths or an unknown `IC` value instead of printing a message and
  silently returning `NULL`.

## Numerical robustness

* `get.Posterior()` computed `exp(-B/2)`, then `log(...)` of that, then
  `exp(...)` again before normalizing — a mathematical identity that added an
  unnecessary extra rounding pass and, more importantly, performed no
  numerical stabilization. It now works directly in log-space and applies the
  standard log-sum-exp trick (subtracting the row-wise maximum before
  exponentiating), so posterior probabilities stay finite even for genes with
  extreme BIC values that would overflow a naive `exp(-B/2)`. Results are
  numerically identical to before for realistic input ranges.

## Portability

* `makeExpData()` and `.makeQPCR()` (internal, unexported helpers used once
  to build the shipped `.rda` data files from the raw Illumina/qPCR files)
  had the original author's absolute local path
  (`/Users/weinhol/GitHub/BGSC/`) and `setwd()` calls hardcoded. They now take
  a `pkgDir` parameter (defaulting to `getwd()`) and never call `setwd()`.
  These functions are not part of the reproducible analysis pipeline.

## Housekeeping

* Removed dead/always-`FALSE` debug code (the `FYI` block in
  `get.log2Mean.and.log2FC()`) and stray commented-out `print()` calls.
  The dead block's check (that the pooled standard error equals the
  Satterthwaite approximation) is now a real regression test in
  `tests/testthat/test-log2FC.R`.
* Removed an unused, side-effect-free call to `normexp.fit.detection.p()` in
  `normalizeExpData()`.
* Added a `tests/testthat` suite covering `InformationCriterion`, `get.IC`,
  `get.Posterior`, `logLikelihoodOfnormData`, the qPCR ddCt sign fix, the
  duplicate-column-name fix, and the Venn set operations (`intersect*`,
  `union*`).
