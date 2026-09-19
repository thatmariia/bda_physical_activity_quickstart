# ==========================================================
# == FUNCTIONS FOR DATA PREPROCESSING
# ==========================================================

#' Fit a preprocessing pipeline to the data
fit_preprocess <- function(df) {
  caret::preProcess(df_feat(df), method = c("nzv", "corr", "center", "scale"))
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