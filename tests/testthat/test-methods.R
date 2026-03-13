# Tests for S3 methods and downstream behavior of NULL meta_analysis.
# All Layer 1 (fake objects, no auto.arima, no skip_on_cran needed).

# ── Helpers ───────────────────────────────────────────────────────────────────

make_fake_sden_result <- function(test_type = "ENT",
                                  selection_mechanism = "auto",
                                  rema_beta = 0.2,
                                  rema_pval = 0.40,
                                  all_sig = 4L, pos = 3L, neg = 1L,
                                  n_effects = 6L, sig_effects = 4L,
                                  pnull = 0.05) {
  btest <- stats::binom.test(x = sig_effects, n = n_effects,
                             p = pnull, alternative = "greater")

  params <- list(
    test_type          = test_type,
    selection_mechanism = selection_mechanism,
    rema_beta          = rema_beta,
    pnull              = pnull,
    rema_pval          = rema_pval,
    all_sig_sum        = all_sig,
    positive_sig_sum   = pos,
    negative_sig_sum   = neg,
    number_of_effects  = n_effects,
    test_pval          = btest$p.value,
    sig_effects        = sig_effects
  )

  result <- list(sden_parameters = params, binomial_test = btest)
  class(result) <- c("sden_results", "list")
  result
}

# ── summary.sden_results: ENT / auto ─────────────────────────────────────────

test_that("summary.sden_results runs without error for ENT (auto)", {
  r <- make_fake_sden_result()
  expect_no_error(capture.output(summary(r)))
})

test_that("summary.sden_results output contains test type label for ENT", {
  r <- make_fake_sden_result()
  expect_output(summary(r), regexp = "ENT")
})

test_that("summary.sden_results output contains 'Automatic' for auto selection", {
  r <- make_fake_sden_result(selection_mechanism = "auto")
  expect_output(summary(r), regexp = "Automatic")
})

test_that("summary.sden_results output contains REMA beta", {
  r <- make_fake_sden_result(rema_beta = 0.2)
  expect_output(summary(r), regexp = "0.2")
})

test_that("summary.sden_results output contains the binomial p-value", {
  r <- make_fake_sden_result()
  expect_output(summary(r), regexp = "p-value")
})

test_that("summary.sden_results returns the object invisibly", {
  r <- make_fake_sden_result()
  captured <- withVisible(summary(r))
  expect_false(captured$visible)
  expect_s3_class(captured$value, "sden_results")
})

# ── summary.sden_results: SDT counter-positive / manual ──────────────────────

test_that("summary.sden_results output contains 'SDT' for SDT counter-positive", {
  r <- make_fake_sden_result(test_type = "SDT counter-positive",
                             selection_mechanism = "SDT",
                             rema_beta = 0.5, rema_pval = 0.02,
                             sig_effects = 1L, pnull = 0.025)
  expect_output(summary(r), regexp = "SDT")
})

test_that("summary.sden_results output contains 'Manual' for manual selection", {
  r <- make_fake_sden_result(test_type = "SDT counter-positive",
                             selection_mechanism = "SDT",
                             rema_beta = 0.5, rema_pval = 0.02,
                             sig_effects = 1L, pnull = 0.025)
  expect_output(summary(r), regexp = "Manual")
})

# ── sden_test errors when meta_analysis is NULL ───────────────────────────────

test_that("sden_test stops with informative error when meta_analysis is NULL", {
  # Simulates the case where too few subjects produced valid models and
  # metafor::rma() failed, leaving meta_analysis = NULL in the iarimax object.
  fake <- list(
    results_df = data.frame(
      id            = c("1", "2"),
      estimate_x    = c(0.5, -0.3),
      "std.error_x" = c(0.1,  0.15),
      n_valid       = c(25L,  25L),
      n_params      = c(2L,   2L),
      stringsAsFactors = FALSE,
      check.names   = FALSE
    ),
    meta_analysis        = NULL,
    error_arimax_skipped = character(0),
    models               = NULL
  )
  class(fake) <- c("iarimax_results", "list")
  attr(fake, "focal_predictor") <- "x"
  attr(fake, "id_var")          <- "id"
  attr(fake, "timevar")         <- "time"

  expect_error(suppressMessages(sden_test(fake)), regexp = "SDEN test stopped")
})
