#' Summary method for sden_results objects.
#'
#' @param object An object of class \code{sden_results}.
#' @param ... Additional arguments (ignored).
#' @returns Invisibly returns \code{object}, called for its side effect of printing.
#' @export
#' @method summary sden_results
summary.sden_results <- function(object, ...) {

  p  <- object$sden_parameters
  tt <- p$test_type

  # -- Header --
  cat("\n")
  cat("                  SDEN Test Results\n")
  cat("               -----------------------\n")
  cat("   Test:", tt,
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
  cat("   Number of significant positive cases:", p$positive_sig_sum, "\n")
  cat("   Number of significant negative cases:", p$negative_sig_sum, "\n")
  cat("   Number of valid cases:", p$number_of_effects, "\n\n")

  relevant_label <- ifelse(tt == "ENT", "all (both sides)",
                    ifelse(tt == "SDT counter-positive", "just negatives", "just positives"))
  cat("   Given the test, relevant significant effects are", relevant_label, ":", p$sig_effects, "\n")
  cat("   Number of valid cases:", p$number_of_effects, "\n\n")

  test_label <- ifelse(tt == "ENT", "ENT", "SDT")
  cat("  ", test_label, "p-value =", round(object$binomial_test$p.value, 7), "\n")

  invisible(object)
}
