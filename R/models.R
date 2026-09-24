# ==========================================================
# == FUNCTIONS FOR TRAINING MODELS
# ==========================================================

#' Construct a specification for fitting one model, based on the provided options
#' @param df The input data frame
#' @param method The `caret::train()` method to fit with
#' @param weighting Whether to weight the epochs by the confidence of their label
#' @param pca The share of the variance the principal components should keep
#' @param corr Whether to remove highly correlated features ("corr" or "nocorr")
#' @return A list with everything `caret::train()` needs for this model
get_fit_spec <- function(df, method, weighting, pca, corr) {
  key <- paste(
    method, weighting, ifelse(pca < 1, paste0("pca", pca), "nopca"), corr,
    sep = "_"
  )

  # Weight the epochs by the confidence of their label (uses weighting option).
  # Unlabelled epochs have no confidence by definition, but they are not less
  # certain than the others, so they keep their full weight.
  weights <- if (weighting == "weighted") {
    ifelse(df$aggr_activity == "-", 1, df$activity_confidence)
  } else {
    rep(1, nrow(df))
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

  return(list(
    x = df_feat(df), y = factor(df$aggr_activity), weights = weights,
    method = method, extra_args = extra_args,
    preprocess = preprocess_methods(pca, corr == "corr"), pca = pca, key = key
  ))
}

#' Fit a list of models based on the provided data and options
#' @param df The input data frame (with a user_id column to group the folds by)
#' @param opts The options for training the models
#' @param number The number of folds for cross-validation
#' @param repeats The number of times to repeat the cross-validation
#' @return A named list of the fitted models
fit_models <- function(df, opts, number = 2, repeats = 1) {
  # Define train control with repeated cross-validation, keeping all epochs of
  # a user in the same fold, since the test data comes from users not seen in training
  # ==> START LLM https://chatgpt.com/share/6ab0e864-2034-83eb-be85-7968bad11e46
  folds <- unlist(
    lapply(seq_len(repeats), function(r) {
      f <- groupKFold(df$user_id, k = number)
      names(f) <- paste0("Fold", seq_along(f), ".Rep", r)
      f
    }),
    recursive = FALSE
  )
  # ==> END LLM

  # Get specifications for training models
  specs <- pmap(opts, \(method, weighting, pca, corr) get_fit_spec(df, method, weighting, pca, corr))

  # Fit the models
  models <- list()
  for (i in seq_along(specs)) {
    spec <- specs[[i]]
    cat("Fitting model", i, "/", length(specs), ":", spec$key, "...\n")

    # caret fits the preprocessing within each fold, so the held out users play no
    # part in it, and `predict()` applies it to new data by itself
    trcntr <- caret::trainControl(
      method = "cv", index = folds, verboseIter = FALSE, allowParallel = TRUE,
      preProcOptions = list(thresh = spec$pca)
    )

    #' Train a model with caret, passing extra arguments through
    train_model <- function(x, y, weights, method, preProcess, trControl, ...) {
      caret::train(
        x = x, y = y, weights = weights, method = method,
        preProcess = preProcess, trControl = trControl, ...
      )
    }

    fit <- tryCatch(
      do.call(train_model, c(
        list(spec$x, spec$y, spec$weights, spec$method, spec$preprocess, trcntr),
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

  return(models)
}

#' Fit all models and return them with their results
#' @param df The input data frame (with a user_id column to group the folds by)
#' @param opts The options for training the models
#' @param number The number of folds for cross-validation
#' @param repeats The number of times to repeat the cross-validation
#' @return A list of the fitted models and their results
fit_all <- function(df, opts, number = 2, repeats = 1) {
  models <- fit_models(df, opts, number, repeats)

  results <- data.frame(
    accuracy = sapply(models, \(x) merge(x$results, x$bestTune)$Accuracy),
    kappa = sapply(models, \(x) merge(x$results, x$bestTune)$Kappa)
  )

  return(list(models = models, results = results))
}

#' Tabulate what varied while each model was trained: the values of its tuning
#' parameters, or the folds themselves for a method with nothing to tune
#' @param models A named list of models fitted with `caret::train()`
#' @return A data frame with one row per model and value
training_table <- function(models) {
  imap_dfr(models, \(fit, key) {
    parameters <- fit$modelInfo$parameters$parameter
    tuned <- parameters[map_lgl(parameters, \(p) n_distinct(fit$results[[p]]) > 1)]

    if (length(tuned) == 0) {
      return(tibble(
        key = key, parameter = "fold", setting = "none",
        value = seq_len(nrow(fit$resample)),
        accuracy = fit$resample$Accuracy, accuracy_sd = NA_real_
      ))
    }

    # the parameter with the most values goes on the x axis, the others tell the lines apart
    x_parameter <- tuned[which.max(map_int(tuned, \(p) n_distinct(fit$results[[p]])))]
    others <- setdiff(tuned, x_parameter)
    settings <- if (length(others) == 0) {
      "none"
    } else {
      fit$results |>
        select(all_of(others)) |>
        imap(\(x, name) paste0(name, "=", x)) |>
        pmap_chr(paste)
    }
    value <- as.numeric(fit$results[[x_parameter]])

    tibble(
      key = key,
      # the penalty is tuned on a log scale
      parameter = if (x_parameter == "lambda") "log10(lambda)" else x_parameter,
      setting = settings,
      value = if (x_parameter == "lambda") log10(value) else value,
      accuracy = fit$results$Accuracy, accuracy_sd = fit$results$AccuracySD
    )
  })
}
