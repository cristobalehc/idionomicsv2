# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

**idionomicsv2** is an R package implementing the I-ARIMAX algorithm — individual-level ARIMAX (AutoRegressive Integrated Moving Average with eXogenous variables) modeling for panel/longitudinal data. It fits separate time-series models per subject, then aggregates results via random-effects meta-analysis.

## Common Commands

All commands assume you are in an R session in the project root:

```r
# Load package for interactive development
devtools::load_all()

# Rebuild documentation from Roxygen comments
devtools::document()

# Run full package check (CRAN-style)
devtools::check()

# Install package locally
devtools::install()

# Run tests
devtools::test()
```

From the shell:
```bash
R CMD build .
R CMD check idionomicsv2_*.tar.gz
```

## Architecture

### Exported Functions

**`iarimax()`** — Core algorithm (`R/iarimax.R`)
- Signature: `iarimax(dataframe, min_n_subject = 20, minvar = 0.01, y_series, x_series, focal_predictor = NULL, id_var, timevar, correlation_method = 'pearson', keep_models = FALSE, verbose = FALSE)`
- Algorithm flow:
  1. **Validation** — checks all named variables exist; timevar must be complete (no NAs)
  2. **Filtering** — removes subjects with too few pairwise-complete observations or insufficient variance in any series
  3. **Per-subject loop** — orders by `timevar`, computes raw correlation, fits ARIMAX via `forecast::auto.arima()`, extracts coefficients via `broom::tidy()`
  4. **Aggregation** — pivots coefficients wide, builds summary dataframe
  5. **Meta-analysis** — runs `metafor::rma()` on the focal predictor coefficients
- Returns S3 class `iarimax_results` (list):
  - `results_df` — per-subject ARIMA orders, coefficients, SE, n_valid, n_params, raw_cor
  - `meta_analysis` — `metafor::rma` object (or NULL if failed)
  - `case_number_detail` — list with `n_original_df`, `n_filtered_out`, `error_arimax_skipped`, `n_used_iarimax`
  - `models` — raw model objects (only if `keep_models = TRUE`, else NULL)
- Attributes: `outcome`, `focal_predictor`, `id_var`, `timevar`

**`i_pval()`** — Per-subject p-values (`R/i_pval.R`)
- Signature: `i_pval(iarimax_object, feature = NULL)`
- Computes two-tailed t-test p-values using ML-based df = n_valid - n_params
- Returns updated `iarimax_results` object with new `pval_<feature>` column in `results_df`
- p-values differ from `lm()` because ARIMA uses ML sigma² (SSR/n), not OLS (SSR/(n-k))

**`sden_test()`** — Sign Divergence / Equisyncratic Null tests (`R/sden_test.R`)
- Signature: `sden_test(iarimax_object, alpha_arimax = 0.05, alpha_binom = NULL, test = "auto", feature = NULL)`
- Modes: `"auto"` (selects ENT/SDT based on pooled REMA p-value at fixed 0.05), `"SDT"`, `"ENT"`
- Auto selection threshold is fixed at 0.05 regardless of `alpha_arimax` — intentional design
- Returns S3 class `sden_results` (list: `sden_parameters`, `binomial_test`)
- Attributes: `outcome`, `focal_predictor`, `id_var`, `timevar` (inherited from iarimax object)

### S3 Methods (`R/methods.R`)

- **`summary.iarimax_results(object, alpha = 0.05, ...)`** — prints subject counts, per-subject direction/significance counts, REMA estimates, heterogeneity; handles NULL meta_analysis gracefully
- **`plot.iarimax_results(x, feature = NULL, y_series_name = NULL, x_series_name = NULL, alpha_crit_t = 0.05, lims = c(-1, 1), ...)`** — caterpillar plot with per-subject CIs (green/red/black), REMA band overlay; returns ggplot object
- **`summary.sden_results(object, ...)`** — prints test type, selection mechanism, hypothesis, and binomial test results

### case_number_detail accounting

The four quantities are cleanly separated and must always satisfy:
`n_original_df == n_filtered_out + length(error_arimax_skipped) + n_used_iarimax`

- `n_filtered_out` = lost to var/n filter (computed as `n_original - length(subjects)`)
- `error_arimax_skipped` = lost to ARIMAX failure during the loop
- `n_used_iarimax` = successfully modeled (`nrow(results_df) - length(exclude)`)

## Key Dependencies

- `forecast` — `auto.arima()` for model selection
- `metafor` — `rma()` for random-effects meta-analysis
- `broom` — model tidying
- `ggplot2`, `forcats` — caterpillar plot
- `dplyr`, `tidyr`, `tibble`, `rlang` — data manipulation and NSE
- `stats`, `utils` — base R

## Test Suite

Two-layer approach throughout:
- **Layer 1** — fake objects, no model fitting, no `skip_on_cran()`; fast
- **Layer 2** — real `iarimax()` output, `skip_on_cran()` inside each test

| File | Coverage |
|------|----------|
| `helper-data.R` | `make_panel(n_subjects, n_obs, seed)` — AR(1) synthetic panel |
| `test-iarimax-validation.R` | Upfront errors: misspelled vars, NA timevar, focal_predictor, correlation_method, filter threshold |
| `test-iarimax-output-structure.R` | Class/attributes, list structure, models, results_df columns, verbose equivalence |
| `test-iarimax-core-correctness.R` | Coefficient values, ARIMA orders, raw_cor, n_valid, temporal ordering, signal recovery |
| `test-iarimax-edge-cases.R` | Filtering boundaries, NA handling, id coercion, multi-predictor, non-sequential timevar |
| `test-i-pval.R` | Formula, p-value range/symmetry, NA handling, non-focal features, df≤0 warning |
| `test-sden-test.R` | Validation, output structure, significance counts, auto/manual selection, pnull, binomial test |
| `test-methods.R` | summary.iarimax_results, summary.sden_results, plot.iarimax_results (guard, ggplot return, lims, labels, NULL meta) |

## Development Notes

- Roxygen2 (v7.3.2) with markdown enabled. Run `devtools::document()` after editing `@` tags.
- `man/` files are auto-generated — do not edit them manually.
- **fable/forecast class workaround** at `class(model) <- setdiff(class(model), "ARIMA")` in `iarimax.R` — do not remove; prevents `broom::tidy()` hijacking by fable.
- **Missing data in y/x** handled implicitly by `stats::arima`'s Kalman Filter — do NOT add manual NA syncing.
- **timevar must be complete** (no NAs) — by design, to force researchers to be mindful about data ordering.
- `utils::globalVariables()` declared in `iarimax.R` (`count`, `var_y`) and `methods.R` (`n_valid`, `n_params`, `df_mod`, `crit_val`, `line_color`) for tidyverse NSE.
- `CLAUDE.md` is in `.Rbuildignore`.

## CRAN Readiness

Current status: `0 errors, 0 warnings, 0 notes` (as of last check)

Remaining before submission:
- Version bump: `0.0.0.9000` → `0.1.0`
- Add `@examples` to all 3 exported functions (`iarimax`, `i_pval`, `sden_test`)
- Run rhub cross-platform check
- Consider adding a vignette (strongly recommended for statistical packages)
- Package name `idionomicsv2` has "v2" suffix — CRAN discourages version numbers in package names
