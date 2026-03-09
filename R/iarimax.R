#' Run I-ARIMAX algorithm.
#'
#' @export
#' @importFrom utils globalVariables
#'
#' @param dataframe Your dataframe.
#' @param min_n_subject The minimum number of non NA cases to run the analyses. It will filter cases with more NA's than the threshold. Defaults to 20.
#' @param minvar The minimum variance for both series (&) to include a case. Defaults to 0.01.
#' @param y_series A string containing the name of your dependent variable y.
#' @param x_series A string containing the name of your independent variable x.
#' @param id_var A string containing your id variable.
#' @param timevar Requiered to arrange timeseries.
#' @param correlation_method Select method fr raw correlations. Options are: 'spearman', 'pearson' or 'kendall'. Defaults to 'pearson'.
#' @param keep_models If TRUE, will keep original arimax models in a list.
#'
#' @returns A list containing a dataframe with the ARIMA parameters, plus the xreg parameter (the beta value for your x_series) together with their std.errors. If metaanalysis = TRUE, will also output a random effects meta analysis. If hlm_compare = TRUE, will also output a model comparison with HLM.


#######################################
############ I ARIMAX FUNCTION #######
#####################################

iarimax <- function(dataframe, min_n_subject = 20, minvar = 0.01, y_series, x_series, id_var,
                    timevar, correlation_method = 'pearson', keep_models = FALSE) {



  # CHeck wether variables are in the in the dataset.
  required_vars <- c(y_series, x_series, id_var, timevar)

  if (!all(required_vars %in% colnames(dataframe))) {
    missing_vars <- required_vars[!required_vars %in% colnames(dataframe)]
    stop(paste("Cannot find required variables. Check if you spelled the following variables correctly:", paste(missing_vars, collapse = ", ")))
  }

  # Convert strings to rlang::symbols
  y_series_sym <- rlang::sym(y_series)
  x_series_sym <- rlang::sym(x_series)
  id_var_sym <- rlang::sym(id_var)
  timevar_sym <- rlang::sym(timevar)

  #Id var as character, to avoid numerics.
  dataframe[[id_var]] <- as.character(dataframe[[id_var]])


  #Subset only relevant features.
  dataframe <- dataframe[, required_vars, drop = FALSE]

  ##                       ##
  #Add a filtering message ##
  ##                       ##

  # Filter N pairwise complete observations with variance conditions
  subjects <- dataframe %>%
    dplyr::group_by(!!id_var_sym) %>% #Group by id variable.
    dplyr::filter(!is.na(!!y_series_sym) & !is.na(!!x_series_sym) & !is.na(!!timevar_sym)) %>% #Filter all that are complete in all variables.
    dplyr::summarise(
      count = dplyr::n(), #Use counts to filter.
      var_y = stats::var(!!y_series_sym, na.rm = TRUE), #use variance to filter.
      var_x = stats::var(!!x_series_sym, na.rm = TRUE), #Use variance to filter.
      .groups = 'drop'
    ) %>%
    dplyr::filter(count >= min_n_subject, var_y >= minvar, var_x >= minvar)  %>%
    dplyr::pull(!!id_var_sym) #Filter data.

  #ID as character, to avoid giant lists.
  subjects <- as.character(subjects)

  #
  # End of filtering message.
  #

  #
  # Running auto-arima algorithm message.
  #

  #Storage lists.

  #Number of valid cases & parameters.
  n_valid <- list()
  n_params <- list()

  #Exclude cases where arimax don't work.
  exclude <- character(0)

  #Model lists.
  arimax_models <- list()
  tidy_models <- list()
  raw_correlation <- list()
  AR_N <- list()
  I_N <- list()
  MA_N <- list()


  #Start case number counter.
  casen = 0
  for(i in subjects) {

    #Update case number.
    casen <- casen + 1

    #Start text.

    # IARIMAX MESSAGE EACH CASE.
    #cat(paste('  Applying auto ARIMAX to case: ', as.character(i), ' ... '))
    #

    #Extract the current subject & arrange timeseries by timevar.
    subject_n <- dataframe %>%
      dplyr::filter(!!id_var_sym == i, !is.na(!!timevar_sym)) %>% #also filter out missing timevars, to avoid sending them to the bottom.
      dplyr::arrange(!!timevar_sym) #Ensure time-series order.

    ##########################################
    ### A COMMENT OF MISSING DATA HANDLING ###
    #########################################

    # Note: We do NOT need to manually sync NAs.
    # stats::arima's Kalman Filter (C source) automatically treats
    # any row with a missing predictor as a missing observation
    # in the state-space update, preserving temporal structure.

    #Extract y vector.
    y_vector <- subject_n %>%
      dplyr::pull(!!y_series_sym)

    #Extract x vector.
    x_vector <- subject_n %>%
      dplyr::pull(!!x_series_sym)

    #Count number of valid observations.
    n_valid_val <- sum(!is.na(y_vector) & !is.na(x_vector))


    ########################
    #### Calculate cor #####
    ########################

    correlation <- tryCatch(
      {#Supress spearman's warning about p-values with ties.
        suppressWarnings(stats::cor.test(x = x_vector, y = y_vector, method = correlation_method))
      },
      error = function(e) {
        cat("\n","\n",' Error computing correlation for case: ',as.character(i), "\n","  ",e$message,"\n")
        NULL #Set model as NULL
      }
    )

    #Handle when correlations are null.
    if (is.null(correlation)) {
      raw_correlation[[i]] <- NA
    } else {
      raw_correlation[[i]] <- correlation$estimate[1]
    }

    ##############################
    ### CALCULATE AUTO.ARIMA ####
    ############################

    model <- tryCatch(
      {
        forecast::auto.arima(y = y_vector, xreg = x_vector, approximation = FALSE, stepwise = FALSE)
      },
      error = function(e) {
        cat("\n","\n",'  Error running ARIMAX model for case: ',as.character(i), "\n","  ",e$message,"\n")
        NULL #Set model as NULL
      }
    )

    #Fill the list with NA if the model is null

    if (is.null(model)) {
      AR_N[[i]] <- NA
      I_N[[i]] <- NA
      MA_N[[i]] <- NA
      n_valid[[i]] <- NA
      n_params[[i]] <- NA
      tidy_models[[i]] <- NULL
      arimax_models[[i]] <- NULL

      exclude <- c(exclude,i)

      ##### MESSAGE ABOUT SKIPPING CASE OR PROGRESS.

      #cat("\n","   Skipping case due to error: ")
      #cat("        ... ",round((casen/(length(subjects))*100),digits = 1),'% completed',"\n","\n") #Keep printing advance percentage.
      next #Skip the next part of this iteration of the loop, so it doesn't get overriden and throws an error.
    }

    # Tidy the model (may be hijacked by fable::tidy.ARIMA)
    tidymodel <- broom::tidy(model)

    # If we didn't get a tibble/data.frame (e.g., fable's tidy.ARIMA returned NULL),
    # fall back by stripping the "ARIMA" class so broom's tidy.Arima is used.
    if (is.null(tidymodel) || !is.data.frame(tidymodel)) {

      if (inherits(model, "Arima")) {

        # Work on a copy so we don't touch `model` used later
        model2 <- model
        class(model2) <- setdiff(class(model2), "ARIMA")

        # Now S3 dispatch sees classes c("forecast_ARIMA", "Arima"),
        # so it will choose broom's tidy.Arima method instead of fable's tidy.ARIMA
        tidymodel <- broom::tidy(model2)

      } else {
        stop(
          "IARIMAXoid_Pro: could not tidy model of class ",
          paste(class(model), collapse = ", ")
        )
      }
    }

    #Add id to tidy model.
    tidymodel <- tidymodel %>%
      dplyr::mutate(!!id_var_sym := i)

    #Fill the tidy model list.
    tidy_models[[i]] <- tidymodel
    # Conditional storage of raw models
    if (keep_models) {
      arimax_models[[i]] <- model
    } else {
      arimax_models[[i]] <- NULL # Placeholder to keep list lengths consistent
    }


    #Fill the number of ARIMA parameteres lists: Just the number of AR I MA processes involved.
    AR_N[[i]] <- model$arma[1] #Fill AR list.
    I_N[[i]] <- model$arma[6] #Fill I list.
    MA_N[[i]] <- model$arma[2] #Fill MA List.

    #Fill n valid and n params.
    n_valid[[i]] <- n_valid_val
    if (!is.null(model)){
      n_params[[i]] <- length(model$coef)}


    ### FINISH PROGRESS TEXT:

    #Finish the text.
    #cat(round((casen/(length(subjects))*100),digits = 1),'% completed',"\n")
  }

  ###############################################
  ### CREATE DATAFRAME WITH PARAMETERS ##########
  ###############################################
  #Filter null models.
  tidy_list <- Filter(Negate(is.null), tidy_models)

  #Stop if no valid model is present.
  if (length(tidy_list) == 0) {
    stop(
      "iarimax: all ARIMAX fits failed after filtering. ",
      "No subject produced a valid model. ",
      "Check (a) time ordering / timevar, (b) missingness in y/x, ",
      "(c) min_n_subject/minvar thresholds, and (d) whether xreg has NA patterns."
    )
  }

  #Create tidy long.
  tidy_long <- dplyr::bind_rows(tidy_list)


  # Pivot the coefficients to wide format
  # This creates: estimate_xreg, std.error_xreg, estimate_ar1, etc.
  tidy_wide <- tidy_long %>%
    tidyr::pivot_wider(
      id_cols     = !!id_var_sym,
      names_from  = "term",
      values_from = c("estimate", "std.error")
    )


  # Convert the lists to vectors
  AR_vector <- unlist(AR_N)
  I_vector <- unlist(I_N)
  MA_vector <- unlist(MA_N)
  n_valid_vector <- unlist(n_valid)
  n_params_vector <- unlist(n_params)
  raw_correlation_vector <- unlist(raw_correlation)


  #Create individual level summaries.
  summary_df <- tibble::tibble(
    !!id_var_sym := subjects,
    p = AR_vector,
    d = I_vector,
    q = MA_vector,
    raw_cor = raw_correlation_vector,
    n_valid = n_valid_vector,
    n_params = n_params_vector
  )

  #Create dataframe to return.

  results_df <- summary_df %>%
    dplyr::left_join(tidy_wide, by = id_var)


  #Add fixed attribute to id column, to be used in other functions.
  attr(results_df[[id_var]], "is_id_column") <- TRUE




  ##############################################
  ##### RUN RANDOM EFFECTS META ANALYSIS ######
  #############################################

  #Try running the random effects meta analysis.

  meta_analysis <-
    tryCatch(
      {
        metafor::rma(yi = results_df$estimate_xreg, sei = results_df$std.error_xreg, method = "REML")
      },
      error = function(e) {
        cat('Error running RME:', e$message, '\n')
        NULL
      }
    )

  final_obj <- list(results_df = results_df,
                    meta_analysis = meta_analysis,
                    error_arimax_skipped = exclude,
                    models = if (keep_models) arimax_models else NULL)

  class(final_obj) <- c("iarimax_results", "list")
  return(final_obj)

  return()

}


utils::globalVariables(c("count", "var_y", "var_x")) #Declare symbolic global variables.
