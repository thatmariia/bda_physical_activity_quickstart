# Feature functions (time domain)

#' Compute the lagged correlation between two signals
#' @param x A numeric vector representing the first signal
#' @param y A numeric vector representing the second signal (default is x)
#' @param lag An integer specifying the lag (default is 0)
#' @return The lagged correlation coefficient between x and y
lagged_cor <- function(x, y = x, lag = 0) {
  # compute correlation between x and a time shifted y
  r_lagged <- cor(x, dplyr::lag(y, lag), use = "pairwise")
  return(r_lagged)
}

#' Extract all time domain features from an accelerometer signal data frame segmented into epochs
#' @param signal_df A data frame containing the signal data with columns: userid, trial, sampleid, X1, X2, X3, activity
#' @param n_samples_per_epoch The number of samples per epoch (default is 128, corresponding to 2.56 seconds at 50 Hz)
#' @return A data frame containing the extracted time domain features for each epoch
get_time_domain_features_acc <- function(signal_df, n_samples_per_epoch = 128) {
  time_domain_features <- signal_df |>
    # partition into epochs and add an epoch ID variable
    mutate(epoch = sampleid %/% n_samples_per_epoch) |>
    # extract statistical features from each epoch
    group_by(epoch) |>
    summarise(
      # Activity label for the epoch = most common value (mode)
      aggr_activity = most_common_value(activity),
      # Confidence of the activity label = proportion of samples with that label
      activity_confidence = value_confidence(aggr_activity, activity),

      # Keep starting sample ID of the epoch
      sampleid = sampleid[1],

      # Signal features
      mean_X1 = mean(X1),
      mean_X2 = mean(X2),
      sd_X1 = sd(X1),
      q25_X1 = quantile(X1, .25),
      skew_X1 = e1071::skewness(X1),
      AR1_X1_lag1 = lagged_cor(X1, lag = 1),
      AR1_X1_lag2 = lagged_cor(X1, lag = 2),
      AR_X1X2_lag1 = lagged_cor(X1, X2, lag = 1),

      # ...
      # ... add your own features here ...
      # ... (to get inspired, look at the histograms above)
      # ...

      # Keep track of epoch lengths (some epochs are less than 128 samples)
      n_samples = n()
    )
  return(time_domain_features)
}

#' Extract all time domain features from a gyroscope signal data frame segmented into epochs
#' @param signal_df A data frame containing the signal data with columns: userid, trial, sampleid, X1, X2, X3, activity
#' @param n_samples_per_epoch The number of samples per epoch (default is 128, corresponding to 2.56 seconds at 50 Hz)
#' @return A data frame containing the extracted time domain features for each epoch
get_time_domain_features_gyro <- function(signal_df, n_samples_per_epoch = 128) {
  time_domain_features <- signal_df |>
    # partition into epochs and add an epoch ID variable
    mutate(epoch = sampleid %/% n_samples_per_epoch) |>
    # extract statistical features from each epoch
    group_by(epoch) |>
    summarise(
      # Activity label for the epoch = most common value (mode)
      aggr_activity = most_common_value(activity),
      # Confidence of the activity label = proportion of samples with that label
      activity_confidence = value_confidence(aggr_activity, activity),

      # Keep starting sample ID of the epoch
      sampleid = sampleid[1],

      # Signal features
      mean_X1 = mean(X1),
      mean_X2 = mean(X2),
      sd_X1 = sd(X1),
      q25_X1 = quantile(X1, .25),
      skew_X1 = e1071::skewness(X1),
      AR1_X1_lag1 = lagged_cor(X1, lag = 1),
      AR1_X1_lag2 = lagged_cor(X1, lag = 2),
      AR_X1X2_lag1 = lagged_cor(X1, X2, lag = 1),

      # ...
      # ... add your own features here ...
      # ... (to get inspired, look at the histograms above)
      # ...

      # Keep track of epoch lengths (some epochs are less than 128 samples)
      n_samples = n(),
    )
  return(time_domain_features)
}
