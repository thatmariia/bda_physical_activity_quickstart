# ==========================================================
# == FUNCTIONS FOR TRAINING MODELS
# ==========================================================

#' Construct a list of specifications for fitting a model
#' based on provided options
get_fit_spec <- function(
  df, pp_utils,
  method, term_mode, term, weighting, pca
) {
  key <- paste(method, term_mode, term, weighting, ifelse(pca < 1, paste0("pca", pca), "nopca"), sep = "_")

  # Add terms (uses term option)
  data <- make_data(df, key, pp_utils, TRUE)

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

  # Create extra arguments (uses method)
  extra_args <- switch(method,
    glmnet = list(
      # alpha = 0 makes it a ridge penalty
      tuneGrid = expand.grid(alpha = 0, lambda = 10^seq(-2.5, -0.5, length.out = 9))
    ),
    knn = list(
      tuneGrid = expand.grid(k = seq(1, 51, by = 2))
    ),
    naive_bayes = list(
      tuneGrid = expand.grid(laplace = 0, usekernel = c(FALSE, TRUE), adjust = c(0.5, 1, 2))
    ),
    list()
  )

  return(list(data = data, formula = formula, weights = data$weights, method = method, extra_args = extra_args, pp_utils = pp_utils, key = key))
}

#' Preprocess and fit a list of models based on the provided data and options
fit_models <- function(df, opts, number = 2, repeats = 1, nstart = 2) {
  # Define train control with repeated cross-validation
  trcntr <- caret::trainControl(method = "repeatedcv", number = number, repeats = repeats, verboseIter = FALSE, allowParallel = TRUE)

  # Construct preprocessing utils, one set per pca option used
  k <- length(unique(df$aggr_activity))
  pp_utils_per_pca <- map(set_names(unique(opts$pca)), \(use_pca) {
    pp <- fit_preprocess(df, pca = use_pca)
    df_pp <- apply_preprocess(df, pp)
    km <- kmeans(df_pp, centers = k, nstart = nstart)
    list(preproc = pp, km = km)
  })

  # Get specifications for training models
  specs <- pmap(opts, \(method, term_mode, term, weighting, pca) {
    get_fit_spec(
      df, pp_utils_per_pca[[as.character(pca)]],
      method, term_mode, term, weighting, pca
    )
  })

  # Fit the models
  models <- list()
  i <- 1
  for (spec in specs) {
    cat("Fitting model", i, "/", length(specs), ":", spec$key, "...\n")
    i <- i + 1

    #' Train a model with caret, passing extra arguments through
    train_model <- function(formula, data, weights, method, trControl, ...) {
      caret::train(formula, data = data, weights = weights, method = method, trControl = trControl, ...)
    }

    fit <- tryCatch(
      do.call(train_model, c(
        list(spec$formula, spec$data, spec$weights, spec$method, trcntr),
        spec$extra_args
      )),
      # ==> START LLM
      # Report the model that failed, but keep fitting the other ones
      error = function(e) {
        cat("  could not fit", spec$key, ":", conditionMessage(e), "\n")
        NULL
      }
      # ==> END LLM
    )
    if (!is.null(fit)) {
      models[[spec$key]] <- fit
    }
  }

  # Keep the preprocessing utils that belong to each fitted model
  pp_utils <- set_names(map(specs, \(spec) spec$pp_utils), map_chr(specs, "key"))
  pp_utils <- pp_utils[names(models)]

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
    accuracy = sapply(models, \(x) merge(x$results, x$bestTune)$Accuracy),
    kappa = sapply(models, \(x) merge(x$results, x$bestTune)$Kappa)
  )

  return(list(models = models, results = results, pp_utils = pp_utils))
}
