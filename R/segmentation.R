# ==========================================================
# == FUNCTIONS FOR SEGMENTATION OF ACTIVITIES
# ==========================================================

#' Compute the most common value (mode) of a vector
#' @param x A vector of values
#' @return The most common value in the vector
most_common_value <- function(x) {
  non_na_x <- x[!is.na(x)]
  if (length(non_na_x) == 0) {
    return("-")
  }
  counts <- table(non_na_x, useNA = "no")
  most_frequent <- which.max(counts)
  return(names(most_frequent))
}

#' Compute the confidence of the most common value in a vector,
#' defined as the proportion of occurrences of the most common value
#' @param mcval The most common value
#' @param x A vector of values
#' @return The confidence of the most common value
value_confidence <- function(mcval, x) {
  if (mcval == "-") {
    return(0)
  }
  counts <- table(x, useNA = "always")
  most_frequent <- which(names(counts) == mcval)
  confidence <- counts[most_frequent] / sum(counts)
  return(confidence)
}
