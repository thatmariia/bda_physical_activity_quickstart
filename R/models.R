# ==========================================================
# == FUNCTIONS FOR TRAINING MODELS
# ==========================================================

#' Construct a list of specifications for fitting a model
#' based on provided options
get_fit_spec <- function(
  df, pp_utils,
  method, term_mode, term, weighting
) {
  key <- paste(method, term_mode, term, weighting, sep = "_")

  # Add terms (uses term option)
  data <- make_data(df, key, pp_utils, TRUE)

  maxNWts = 25000 # 1000
  maxit = 100

  # Add weights (uses weighting option)
  data <- data |>
    mutate(
      weights = if (weighting == "weighted") {
        df$activity_confidence
      } else {
        rep(1, nrow(df))
      }
    )

  # Create formula (uses term_mode and term options)
  if (term_mode == "inter" && term == "km_cluster") {
    formula <- aggr_activity ~ (. - weights) * cluster
  } else if (term_mode == "inter" && term == "km_dist") {
    feat_names <- names(apply_preprocess(df, pp_utils$preproc))
    dist_names <- grep("^km_dist_", names(data), value = TRUE)
    # ==> START LLM https://chatgpt.com/share/6aaab81f-4d5c-83eb-bfac-021950015cd3
    rhs <- paste0(
      "(", paste(feat_names, collapse = " + "), ") * ",
      "(", paste(dist_names, collapse = " + "), ")"
    )
    formula <- as.formula(paste("aggr_activity ~", rhs))
    # ==> END LLM
  } else {
    formula <- aggr_activity ~ (. - weights)
  }

  return(list(data = data, formula = formula, method = method, MaxNWts = maxNWts, maxit = maxit, key = key))
}

#' Preprocess and fit a list of models based on the provided data and options
fit_models <- function(df, opts, number = 2, repeats = 1, nstart = 2) {
  # Define train control with repeated cross-validation
  trcntr <- caret::trainControl(method = "repeatedcv", number = number, repeats = repeats, verboseIter = FALSE, allowParallel = TRUE)

  # Precompute relevant data
  pp <- fit_preprocess(df)
  df_pp <- apply_preprocess(df, pp)
  k <- length(unique(df$aggr_activity))
  km <- kmeans(df_pp, centers = k, nstart = nstart)

  # Construct preprocessing utils
  pp_utils <- list(preproc = pp, km = km)

  # Get specifications for training models
  specs <- pmap(opts, \(method, term_mode, term, weighting) {
    get_fit_spec(
      df, pp_utils,
      method, term_mode, term, weighting
    )
  })

  # Fit the models
  models <- list()
  i <- 1
  for (spec in specs) {
    cat("Fitting model", i, "/", length(specs), ":", spec$key, "...\n")
    i <- i + 1
    model <- caret::train(
      spec$formula,
      data = spec$data,
      weights = weights,
      method = spec$method,
      trControl = trcntr,
      MaxNWts = spec$MaxNWts,
      maxit = spec$maxit,
      trace = FALSE
    )
    models[[spec$key]] <- model
  }

  return(list(models = models, pp_utils = pp_utils))
}

#' Fit all models and return the results
#' @param df The input data frame
#' @param opts The options for training the models
#' @param number The number of folds for cross-validation
#' @param repeats The number of times to repeat the cross-validation
#' @param nstart The number of random starts for k-means
#' @return A list of fitted models, their results, and preprocessing utils
fit_all <- function(df, opts, number = 2, repeats = 1, nstart = 2) {
  fitted_models <- fit_models(df, opts, number, repeats, nstart)
  models <- fitted_models$models
  pp_utils <- fitted_models$pp_utils

  results <- data.frame(
    accuracy = sapply(models, function(x) max(x$results$Accuracy)),
    kappa = sapply(models, function(x) max(x$results$Kappa))
  )

  return(list(models = models, results = results, pp_utils = pp_utils))
}
