#' Summary method for iarimax_results objects.
#'
#' @param object An object of class \code{iarimax_results}.
#' @param alpha Significance threshold used for per-subject positive/negative
#'   significant counts. Defaults to 0.05.
#' @param ... Additional arguments (ignored).
#' @returns Invisibly returns \code{object}, called for its side effect of printing.
#' @export
#' @method summary iarimax_results
summary.iarimax_results <- function(object, alpha = 0.05, ...) {

  focal  <- attr(object, "focal_predictor")
  id_var <- attr(object, "id_var")
  outcome <- attr(object, "outcome")

  n_original <- object$case_number_detail$n_original_df
  n_included <- object$case_number_detail$n_used_iarimax
  n_filtered <- object$case_number_detail$n_filtered_out
  n_skipped  <- length(object$case_number_detail$error_arimax_skipped)

  # Compute per-subject p-values (suppress df<=0 warnings — user already gets
  # those from i_pval directly; no need to repeat them in a summary call).
  obj_pval <- suppressWarnings(i_pval(object, feature = focal))
  est_col  <- paste0("estimate_", focal)
  pval_col <- paste0("pval_",     focal)
  df       <- obj_pval$results_df

  # Per-subject direction and significance counts.
  n_positive <- sum(df[[est_col]] > 0, na.rm = TRUE)
  n_negative <- sum(df[[est_col]] < 0, na.rm = TRUE)
  n_pos_sig  <- sum(df[[est_col]] > 0 & df[[pval_col]] < alpha, na.rm = TRUE)
  n_neg_sig  <- sum(df[[est_col]] < 0 & df[[pval_col]] < alpha, na.rm = TRUE)

  # -- Header --
  cat("\n")
  cat("              IARIMAX Results Summary\n")
  cat("           ----------------------------\n")
  cat("   Features:\n")
  cat("     Outcome          :", outcome, "\n")
  cat("     Focal predictor  :", focal, "\n")
  cat("     ID variable      :", id_var, "\n\n")

  # -- Subjects --
  cat("   Subjects:\n")
  cat("     Original Dataframe          :", n_original, "\n")
  cat("     Excl first filter (var or n):", n_filtered, "\n")
  cat("     Skipped (auto.arima failed) :", n_skipped, "\n")
  cat("     Analyzed cases              :", n_included, "\n\n")


  # -- Per-subject effects --
  cat("   Per-subject effects  (alpha =", alpha, "):\n", sep = "")
  cat("     Positive            :", n_positive,"\n")
  cat("     Significant positive:", n_pos_sig, "\n")
  cat("     Negative            :", n_negative,"\n")
  cat("     Significant negative:", n_neg_sig, "\n\n")

  # -- REMA and heterogeneity --
  if (is.null(object$meta_analysis)) {

    cat("   Random-Effects Meta-Analysis : not available (model failed)\n\n")

  } else {

    ma   <- object$meta_analysis
    pred <- predict(ma)

    cat("   Random-Effects Meta-Analysis (REML)\n")
    cat("     beta =", round(as.numeric(ma$beta), 4),
        "  SE =",   round(as.numeric(ma$se),   4),
        "  z =",    round(as.numeric(ma$zval),  4),
        "  p =",    round(as.numeric(ma$pval),  4), "\n")
    cat("     95% Confidence Interval: [", round(ma$ci.lb, 4), ", ", round(ma$ci.ub, 4), "]\n", sep = "")
    cat("     95% Prediction Interval: [", round(pred$pi.lb, 4), ", ", round(pred$pi.ub, 4), "]\n\n", sep = "")

    cat("   Heterogeneity\n")
    cat("     tau2 =", round(ma$tau2, 4),
        "  I2 =",  round(ma$I2, 2), "%",
        "  Q(df =", ma$k - 1, ") =", round(ma$QE, 4),
        "  p =",  round(ma$QEp, 4), "\n\n")

  }

  cat("\n")
  invisible(object)
}


#' Summary method for sden_results objects.
#'
#' @param object An object of class \code{sden_results}.
#' @param ... Additional arguments (ignored).
#' @returns Invisibly returns \code{object}, called for its side effect of printing.
#' @export
#' @method summary sden_results
summary.sden_results <- function(object, ...) {

  #Import attributes.
  focal  <- attr(object, "focal_predictor")
  id_var <- attr(object, "id_var")
  outcome <- attr(object, "outcome")

  p  <- object$sden_parameters
  tt <- p$test_type

  # -- Header --
  cat("\n")
  cat("                  SDEN Test Results\n")
  cat("               -----------------------\n")
  cat("   Features:\n")
  cat("     Outcome          :", outcome, "\n")
  cat("     Focal predictor  :", focal, "\n")
  cat("     ID variable      :", id_var, "\n\n")
  cat("   Test          :", tt,
      ifelse(tt == "ENT", "(Equisyncratic Null Test)", "(Sign Divergence Test)"), "\n")
  cat("   Test selection:",
      ifelse(p$selection_mechanism == "auto", "Automatic", "Manual"), "\n\n")

  # -- Explanation --
  cat("   Explanation:\n")
  cat("  -------------\n")

  if (p$selection_mechanism == "auto") {
    cat("   Automatic selection based on:\n")
    if (tt == "ENT") {
      cat("   The pooled effect (", round(p$rema_beta, 4),
          ") was not statistically significant (p = ", round(p$rema_pval, 4),
          ")\n   or was exactly zero.\n", sep = "")
    } else if (tt == "SDT counter-positive") {
      cat("   The pooled effect (", round(p$rema_beta, 4),
          ") was positive and statistically significant (p = ", round(p$rema_pval, 4), ").\n", sep = "")
    } else {
      cat("   The pooled effect (", round(p$rema_beta, 4),
          ") was negative and statistically significant (p = ", round(p$rema_pval, 4), ").\n", sep = "")
    }
  } else {
    cat("   You manually selected the test.\n")
  }

  # -- Hypothesis --
  cat("\n")
  if (tt == "ENT") {
    cat("   Testing whether the number of all significant effects (both sides)\n",
        "   is greater than", p$pnull, "\n")
  } else if (tt == "SDT counter-positive") {
    cat("   Testing whether the number of negative significant effects\n",
        "   is greater than", p$pnull, "\n")
  } else {
    cat("   Testing whether the number of positive significant effects\n",
        "   is greater than", p$pnull, "\n")
  }

  # -- Results --
  cat("\n   Results:\n")
  cat("  ---------\n")
  cat("   Context - REMA pooled effect (p-val):", round(p$rema_beta, 3),
      "(", round(p$rema_pval, 4), ")\n\n")
  cat("   Total number of significant effects (both sides):", p$all_sig_sum, "\n")
  cat("   Number of significant positive cases            :", p$positive_sig_sum, "\n")
  cat("   Number of significant negative cases            :", p$negative_sig_sum, "\n")
  cat("   Number of valid cases                           :", p$number_of_effects, "\n\n")

  relevant_label <- ifelse(tt == "ENT", "all (both sides)",
                    ifelse(tt == "SDT counter-positive", "just negatives", "just positives"))
  cat("   Given the test, relevant significant effects are", relevant_label, ":", p$sig_effects, "\n")

  test_label <- ifelse(tt == "ENT", "ENT", "SDT")
  cat("  ", test_label, "p-value =", round(object$binomial_test$p.value, 7), "\n")

  invisible(object)
}
