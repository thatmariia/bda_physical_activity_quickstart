# ==========================================================
# == FUNCTIONS FOR DATA CONSTRUCTION FOR MODELS
# ==========================================================

#' Get data with features cols only (excluding outcomes and confidence cols)
df_feat <- function(df) {
  df |> select(-any_of(c("aggr_activity", "activity_confidence")))
}

#' Get a data frame with the k-means distances for each data point
get_kmeans_dist_df <- function(df, km) {
  dists <- proxy::dist(as.matrix(df), as.matrix(km$centers))
  df_km_dists <- unclass(dists) |> data.frame()
  # Change col names to start with "km_dist_"
  names(df_km_dists) <- paste0("km_dist_", names(df_km_dists))
  return(df_km_dists)
}

#' Get a vector of the nearest k-means cluster for each data point
get_kmeans_clusters <- function(df, km) {
  # ==> START LLM https://chatgpt.com/share/6aab3c44-4174-83eb-b07d-5f27ba5ca1bb
  dists <- get_kmeans_dist_df(df, km)
  nearest_cluster <- max.col(-dists, ties.method = "first")
  return(nearest_cluster)
  # ==> END LLM
}

#' Fit a preprocessing pipeline to the data
fit_preprocess <- function(df) {
  caret::preProcess(df_feat(df), method = c("nzv", "corr", "center", "scale"))
}

#' Apply a preprocessing pipeline to the data
apply_preprocess <- function(df, pp) {
  predict(pp, newdata = df_feat(df))
}

#' Prepare data for model training
#' @param df The input data frame
#' @param key The key name of the model for the data configuration
#' @param pp_utils The list of preprocessing utilities (preproc, km - k-means object)
#' @param inc_outcome Whether to include the outcome variable (default: FALSE)
#' @return A data frame ready for training a model with key name
make_data <- function(df, key, pp_utils, inc_outcome = FALSE) {
  x <- apply_preprocess(df, pp_utils$preproc)

  # Add columns for additional terms
  if (grepl("km_cluster", key)) {
    x <- x |> mutate(cluster = get_kmeans_clusters(x, pp_utils$km))
  } else if (grepl("km_dist", key)) {
    x <- x |> bind_cols(get_kmeans_dist_df(x, pp_utils$km))
  }

  # Add a placeholder column for weights
  x <- x |> mutate(weights = rep(1, nrow(x)))

  # Add the outcome variable if requested
  if (inc_outcome) {
    x <- x |> mutate(aggr_activity = df$aggr_activity)
  }

  return(x)
}
