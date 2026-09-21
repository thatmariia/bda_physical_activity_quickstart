# ==========================================================
# == FUNCTIONS FOR DATA PREPROCESSING
# ==========================================================

#' Get data with features cols only (excluding outcomes and confidence cols)
df_feat <- function(df) {
  df |> select(-any_of(c("aggr_activity", "activity_confidence")))
}

#' Fit a preprocessing pipeline to the data
#' @param df The input data frame
#' @param pca The share of the variance the principal components should keep
#'   (1 = keep the features themselves, no pca)
#' @return A preprocessing pipeline
fit_preprocess <- function(df, pca = 1) {
  methods <- c("nzv", "corr", "center", "scale")
  if (pca < 1) {
    methods <- c(methods, "pca")
  }
  caret::preProcess(df_feat(df), method = methods, thresh = pca)
}

#' Apply a preprocessing pipeline to the data
apply_preprocess <- function(df, pp) {
  predict(pp, newdata = df_feat(df))
}

#' Find redundant columns in the data
find_redundant_cols <- function(df) {
  df_numeric <- df |> select_if(is.numeric)
  linear_combos <- caret::findLinearCombos(df_numeric)
  redundant <- names(df_numeric)[linear_combos$remove]
  return(redundant)
}
