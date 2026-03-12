#' Case by case p-value calculation based on t-distribution.
#'
#' @export
#'
#' @param iarimax_object An iarimax object.
#' @param feature Feature name to calculate p-value. Defaults to iarimax focal predictor. Use your original name, function will automatically append "estimate_"
#' @returns Returns an updated version of results_df within the iarimax_object with p-values for the specific feature.

i_pval <- function(iarimax_object, feature = NULL) {

  #Only works for iarimax object now.
  if (!inherits(iarimax_object, "iarimax_results")) { #will add other classes if needed.
    stop("iarimax_object must be an iarimax_results object.")
  }

  #Defaults to focal predictor if feature is null.
  if (is.null(feature)) {
    feature <- attr(iarimax_object, "focal_predictor")
  }

  #Construct the column names based on the 'feature' argument
  feature_name <- paste0("estimate_", feature)
  std_feature_name <- paste0("std.error_", feature)
  pval_col_name <- paste0("pval_",feature)

  #Check if the necessary columns exist
  if (!feature_name %in% names(iarimax_object$results_df)) {
    stop(paste0("Coefficient column '", feature, "' not found in the results data frame."))
  }
  if (!std_feature_name %in% names(iarimax_object$results_df)) {
    stop(paste0("Standard error column '", std_feature_name, "' not found."))
  }

  #Calculate p_value
  iarimax_object$results_df[[pval_col_name]] <- 2*stats::pt(-abs(iarimax_object$results_df[[feature_name]] /
                                                                   iarimax_object$results_df[[std_feature_name]]),
                                                            iarimax_object$results_df$n_valid - iarimax_object$results_df$n_params)


  #Return the modified object
  return(invisible(iarimax_object))
}
