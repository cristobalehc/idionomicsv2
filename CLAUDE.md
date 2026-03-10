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

# Run tests (once tests/ directory exists)
devtools::test()
```

From the shell:
```bash
R CMD build .
R CMD check idionomicsv2_*.tar.gz
```

## Architecture

The package currently has a single exported function `iarimax()` in `R/iarimax.R`.

**Algorithm flow:**
1. **Validation** — checks that all named variables exist in the dataframe
2. **Filtering** — removes subjects with too few observations (`min_n_subject`) or insufficient variance (`minvar`) in either series
3. **Per-subject loop** — for each subject: orders by `timevar`, computes raw correlation, fits ARIMAX via `forecast::auto.arima()`, extracts p/d/q parameters and the `xreg` coefficient using `broom::tidy()`
4. **Aggregation** — pivots coefficients wide, builds a summary dataframe
5. **Meta-analysis** — runs `metafor::rma()` random-effects meta-analysis on the xreg coefficients across subjects

**Return value:** S3 class `iarimax_results` with fields:
- `results_df` — per-subject ARIMA order params and coefficients
- `meta_analysis` — `metafor::rma` object
- `error_arimax_skipped` — list of subject IDs where modeling failed
- `models` — original model objects (only if `keep_models = TRUE`)

## Key Dependencies

- `forecast` — `auto.arima()` for model selection
- `metafor` — `rma()` for random-effects meta-analysis
- `broom` — model tidying
- `dplyr`, `tidyr`, `tibble`, `rlang` — data manipulation and NSE

## Development Notes

- Roxygen2 (v7.3.2) with markdown enabled. Run `devtools::document()` after editing `@` tags.
- `man/` is empty — all `.Rd` files are generated; do not edit them manually.
- `DESCRIPTION` still has placeholder author/title/description fields that need finalizing before CRAN submission.
- No test suite exists yet (`tests/` directory absent).
- Missing data in time series is handled implicitly by `stats::arima`'s Kalman Filter (passed through `forecast::auto.arima`).
- A class dispatch workaround exists in `iarimax.R` for `fable`/`forecast` package compatibility — do not remove it.
