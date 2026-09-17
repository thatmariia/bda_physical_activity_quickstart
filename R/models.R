# ==========================================================
# == FUNCTIONS FOR DEFINING AND TRAINING MODELS
# ==========================================================

#' Get feature data excluding outcomes and confidence cols
#' @param df The input data frame
#' @return A data frame with only the feature columns
df_feat <- function(df) {
  df |> select(-any_of(c("aggr_activity", "activity_confidence")))
}

fit_preprocess <- function(df) {
  caret::preProcess(df_feat(df), method = c("nzv", "corr", "center", "scale"))
}

apply_preprocess <- function(df, pp) {
  predict(pp, newdata = df_feat(df))
}

make_data <- function(df, key, pp_utils, inc_outcome = FALSE) {
  x <- apply_preprocess(df, pp_utils$preproc)

  if (grepl("km_cluster", key)) {
    x <- x |> mutate(cluster = get_kmeans_clusters(x, pp_utils$km))
  } else if (grepl("km_dist", key)) {
    x <- x |> bind_cols(get_kmeans_dist_df(x, pp_utils$km))
  }

  x <- x |> mutate(weights = rep(1, nrow(x)))

  if (inc_outcome) {
    x <- x |> mutate(aggr_activity = df$aggr_activity)
  }

  return(x)
}

#' Get a data frame with feature cols and additional
#' cols for distances between data points and k-means centers
#' @param df The input data frame
#' @param km The k-means object
#' @return A data frame with the k-means distances
get_kmeans_dist_df <- function(df, km) {
  dists <- proxy::dist(as.matrix(df), as.matrix(km$centers))
  df_km_dists <- unclass(dists) |> data.frame()
  # Change col names to start with "km_dist_"
  names(df_km_dists) <- paste0("km_dist_", names(df_km_dists))
  return(df_km_dists)
}

#' Get the nearest k-means cluster for each data point
#' @param df The input data frame
#' @param km The k-means object
#' @return A vector with the nearest cluster for each data point
get_kmeans_clusters <- function(df, km) {
  # ==> START LLM https://chatgpt.com/share/6aab3c44-4174-83eb-b07d-5f27ba5ca1bb
  dists <- get_kmeans_dist_df(df, km)
  nearest_cluster <- max.col(-dists, ties.method = "first")
  return(nearest_cluster)
  # ==> END LLM
}

#' Define options for model training, which are used to
#' construct different model configurations
#' @return A data frame with the training options
define_train_options <- function() {
  opts <- expand_grid(
    method = c("multinom"),
    term_mode = c("nomode", "add"),
    term = c("noterm", "km_cluster", "km_dist"),
    weighting = c("unweighted")
    # method = c("multinom", "lda", "naive_bayes"),
    # term_mode = c("nomode", "inter", "add"),
    # term = c("noterm", "km_cluster", "km_dist"),
    # weighting = c("unweighted", "weighted")
  ) |>
    filter(!(term_mode == "nomode" & term != "noterm")) |>
    filter(!(term == "noterm" & term_mode != "nomode")) |>
    filter(!(weighting == "weighted" & method != "multinom")) |>
    # filter(!(method == "knn" & term != "noterm")) |>
    filter(term_mode != "inter") # For now, taking too long
  return(opts)
}

fit_preprocess <- function(df) {
  caret::preProcess(
    df_feat(df),
    method = c("nzv", "corr", "center", "scale")
  )
}

#' Construct a list of specifications for fitting a model
#' based on provided options
#' @param df The input data frame
#' @param pp_utils The list of preprocessing utilities
#' @param method The method for training
#' @param term_mode The term mode option
#' @param term The term option
#' @param weighting The weighting option
#' @return A list of specifications for fitting a model
get_fit_spec <- function(
  df, pp_utils,
  method, term_mode, term, weighting
) {
  key <- paste(method, term_mode, term, weighting, sep = "_")

  # Add terms
  data <- make_data(df, key, pp_utils, TRUE)

  # Add weights
  data <- data |>
    mutate(
      weights = if (weighting == "weighted") {
        df$activity_confidence
      } else {
        rep(1, nrow(df))
      }
    )

  # Create formula
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

  return(list(data = data, formula = formula, method = method, key = key))
}

#' Preprocess and fit a list of models based on the provided data and options
#' @param df The input data frame
#' @param number The number of folds for cross-validation
#' @param repeats The number of times to repeat the cross-validation
#' @param nstart The number of random starts for k-means
#' @return A list of fitted models
fit_models <- function(df, number = 2, repeats = 1, nstart = 2) {
  # Define train control with repeated cross-validation
  trcntr <- caret::trainControl(method = "repeatedcv", number = number, repeats = repeats, verboseIter = FALSE, allowParallel = TRUE)

  # Precompute relevant data
  pp <- fit_preprocess(df)
  df_pp <- apply_preprocess(df, pp)
  k <- length(unique(df$aggr_activity))
  km <- kmeans(df_pp, centers = k, nstart = nstart)

  pp_utils <- list(preproc = pp, km = km)

  # Define options for training
  opts <- define_train_options()

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
      MaxNWts = 25000,
      maxit = 300,
      trace = FALSE
    )
    models[[spec$key]] <- model
  }

  return(list(models = models, pp_utils = pp_utils))
}

#' Fit all models and return the results
#' @param df The input data frame
#' @param weights The weights for the model (NULL for default)
#' @param number The number of folds for cross-validation
#' @param repeats The number of times to repeat the cross-validation
#' @param nstart The number of random starts for k-means
#' @return A list of fitted models and their results
fit_all <- function(df, number = 2, repeats = 1, nstart = 2) {
  fitted_models <- fit_models(df, number, repeats, nstart)
  models <- fitted_models$models
  pp_utils <- fitted_models$pp_utils

  results <- data.frame(
    accuracy = sapply(models, function(x) max(x$results$Accuracy)),
    kappa = sapply(models, function(x) max(x$results$Kappa))
  )

  return(list(models = models, results = results, pp_utils = pp_utils))
}
