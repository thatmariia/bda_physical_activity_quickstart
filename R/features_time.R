# ==========================================================
# == FUNCTIONS FOR TIME DOMAIN FEATURES
# ==========================================================

#' Compute the lagged correlation between two signals
lagged_cor <- function(x, y = x, lag = 0) {
  # compute correlation between x and a time shifted y
  r_lagged <- cor(x, dplyr::lag(y, lag), use = "pairwise")
  return(r_lagged)
}

#' Compute the root mean square of a signal
rms <- function(x) {
  sqrt(mean(x^2))
}

#' Compute the magnitude area of a signal
sma <- function(x, y, z) {
  mean(abs(x) + abs(y) + abs(z))
}

#' Compute the mean absolute first difference of a signal
mean_abs_diff <- function(x) {
  mean(abs(diff(x)))
}

#' Compute the maximum absolute first difference of a signal
max_abs_diff <- function(x) {
  max(abs(diff(x)))
}

#' Compute the difference between the means of first and last k% of values of a signal
diff_start_end <- function(x, k = 0.1) {
  if (length(x) < 2) {
    return(0)
  }
  n <- length(x)
  m <- max(1, floor(n * k))
  return(mean(tail(x, m)) - mean(head(x, m)))
}

#' Compute the difference in orientation between three signals
diff_orientation <- function(x, y, z) {
  diffs <- c(
    diff_start_end(x),
    diff_start_end(y),
    diff_start_end(z)
  )
  return(sqrt(sum(diffs^2)))
}

#' Compute the jerk of a signal
jerk <- function(x, y, z) {
  sqrt(diff(x)^2 + diff(y)^2 + diff(z)^2)
}

#' Compute the slope of a signal
slope <- function(x) {
  if (length(x) < 2) {
    return(0)
  }
  # Compute the slope using linear regression
  fit <- lm(x ~ seq_along(x))
  return(coef(fit)[2])
}

#' Compute the dominance of a signal
dominance <- function(x, y, z) {
  max(c(sd(x), sd(y), sd(z))) / sum(c(sd(x), sd(y), sd(z)))
}

#' Compute the integrated rotation of a signal
rotation <- function(x, fs) {
  sum(x) / fs
}

#' Compute the absolute integrated rotation of a signal
abs_rotation <- function(x, fs) {
  sum(abs(x)) / fs
}

#' Compute the asymmetry of the rotation of a signal
rotation_asym <- function(x, fs) {
  pos <- sum(pmax(0, x)) / fs
  neg <- sum(abs(pmin(0, x))) / fs
  if (pos + neg == 0) {
    return(0)
  }
  return((pos - neg) / (pos + neg))
}

#' Compute rate of zero crossings in a signal
zero_cross_rate <- function(x) {
  if (length(x) < 2) {
    return(0)
  }
  side <- sign(x - mean(x))
  return(mean(side[-1] != side[-length(side)]))
}

#' Compute the angle of the gravity vector
gravity_angle <- function(x, y, z) {
  mag <- sqrt(mean(x)^2 + mean(y)^2 + mean(z)^2)
  return(acos(mean(x) / mag) * 180 / pi)
}

#' Compute the range of the gravity angle
gravity_angle_range <- function(x, mag) {
  angles <- acos(pmin(1, pmax(-1, x / mag))) * 180 / pi
  if (all(is.na(angles))) {
    return(0)
  }
  return(max(angles) - min(angles))
}

#' Compute the change in gravity angle
gravity_angle_change <- function(x, y, z, k = 0.1) {
  if (length(x) < 2) {
    return(0)
  }
  m <- max(1, floor(length(x) * k))
  start <- gravity_angle(head(x, m), head(y, m), head(z, m))
  end <- gravity_angle(tail(x, m), tail(y, m), tail(z, m))
  return(end - start)
}

#' Compute the angle between two signals
sigangle <- function(x, y) {
  atan2(mean(x), mean(y)) * 180 / pi
}

#' Compute the variance ratio of a signal
var_ratio <- function(x, y, z) {
  var(x) / (var(x) + var(y) + var(z))
}

#' Compute the entropy of the distribution of signal values
entropy <- function(x, bins = 10) {
  if (length(unique(x)) < 2) {
    return(0)
  }
  p <- hist(x, breaks = bins, plot = FALSE)$counts
  p <- p / sum(p)
  return(-sum(p[p > 0] * log(p[p > 0])))
}

#' Extract all time domain features from a signal data frame segmented into epochs
#' @param signal_df A data frame containing the signal data with columns: userid, trial, sampleid, X1, X2, X3, activity
#' @param n_samples_per_epoch The number of samples per epoch (default is 128, corresponding to 2.56 seconds at 50 Hz)
#' @param sample_rate The sample rate of the signal (default is 50 Hz)
#' @return A data frame containing the extracted time domain features for each epoch
get_time_domain_features <- function(signal_df, n_samples_per_epoch = 128, sample_rate = 50) {
  time_domain_features <- signal_df |>
    # partition into epochs and add an epoch ID variable
    mutate(
      epoch = sampleid %/% n_samples_per_epoch,
      acc_mag = sqrt(X1^2 + X2^2 + X3^2)
    ) |>
    # extract statistical features from each epoch
    group_by(epoch) |>
    summarise(
      # Activity label for the epoch = most common value (mode)
      aggr_activity = most_common_value(activity),
      # Confidence of the activity label = proportion of samples with that label
      activity_confidence = value_confidence(aggr_activity, activity),

      # Keep starting sample ID of the epoch
      sampleid = sampleid[1],

      # X1
      mean_X1 = mean(X1),
      median_X1 = median(X1),
      sd_X1 = sd(X1),
      min_X1 = min(X1),
      max_X1 = max(X1),
      rms_X1 = rms(X1),
      q05_X1 = quantile(X1, 0.05),
      q25_X1 = quantile(X1, 0.25),
      q75_X1 = quantile(X1, 0.75),
      q95_X1 = quantile(X1, 0.95),
      var_rat_X1 = var_ratio(X1, X2, X3),
      skew_X1 = e1071::skewness(X1),
      kurtosis_X1 = e1071::kurtosis(X1),
      ar1_lag1_X1 = lagged_cor(X1, lag = 1),
      ar1_lag2_X1 = lagged_cor(X1, lag = 2),
      mean_abs_diff_X1 = mean_abs_diff(X1),
      max_abs_diff_X1 = max_abs_diff(X1),
      slope_X1 = slope(X1),
      diff_X1 = diff_start_end(X1),
      rot_X1 = rotation(X1, sample_rate),
      rot_abs_X1 = abs_rotation(X1, sample_rate),
      rot_asym_X1 = rotation_asym(X1, sample_rate),
      zcr_X1 = zero_cross_rate(X1),
      grav_range_X1 = gravity_angle_range(X1, acc_mag),
      grav_X1 = gravity_angle(X1, X2, X3),
      grav_change_X1 = gravity_angle_change(X1, X2, X3),
      entropy_X1 = entropy(X1),

      # X2
      mean_X2 = mean(X2),
      median_X2 = median(X2),
      sd_X2 = sd(X2),
      min_X2 = min(X2),
      max_X2 = max(X2),
      rms_X2 = rms(X2),
      q05_X2 = quantile(X2, 0.05),
      q25_X2 = quantile(X2, 0.25),
      q75_X2 = quantile(X2, 0.75),
      q95_X2 = quantile(X2, 0.95),
      var_rat_X2 = var_ratio(X2, X1, X3),
      skew_X2 = e1071::skewness(X2),
      kurtosis_X2 = e1071::kurtosis(X2),
      ar1_lag1_X2 = lagged_cor(X2, lag = 1),
      ar1_lag2_X2 = lagged_cor(X2, lag = 2),
      mean_abs_diff_X2 = mean_abs_diff(X2),
      max_abs_diff_X2 = max_abs_diff(X2),
      slope_X2 = slope(X2),
      diff_X2 = diff_start_end(X2),
      rot_X2 = rotation(X2, sample_rate),
      rot_abs_X2 = abs_rotation(X2, sample_rate),
      rot_asym_X2 = rotation_asym(X2, sample_rate),
      zcr_X2 = zero_cross_rate(X2),
      grav_range_X2 = gravity_angle_range(X2, acc_mag),
      grav_X2 = gravity_angle(X2, X3, X1),
      grav_change_X2 = gravity_angle_change(X2, X3, X1, k = 0.1),
      entropy_X2 = entropy(X2),

      # X3
      mean_X3 = mean(X3),
      median_X3 = median(X3),
      sd_X3 = sd(X3),
      min_X3 = min(X3),
      max_X3 = max(X3),
      rms_X3 = rms(X3),
      q05_X3 = quantile(X3, 0.05),
      q25_X3 = quantile(X3, 0.25),
      q75_X3 = quantile(X3, 0.75),
      q95_X3 = quantile(X3, 0.95),
      var_rat_X3 = var_ratio(X3, X1, X2),
      skew_X3 = e1071::skewness(X3),
      kurtosis_X3 = e1071::kurtosis(X3),
      ar1_lag1_X3 = lagged_cor(X3, lag = 1),
      ar1_lag2_X3 = lagged_cor(X3, lag = 2),
      mean_abs_diff_X3 = mean_abs_diff(X3),
      max_abs_diff_X3 = max_abs_diff(X3),
      slope_X3 = slope(X3),
      diff_X3 = diff_start_end(X3),
      rot_X3 = rotation(X3, sample_rate),
      rot_abs_X3 = abs_rotation(X3, sample_rate),
      rot_asym_X3 = rotation_asym(X3, sample_rate),
      zcr_X3 = zero_cross_rate(X3),
      grav_range_X3 = gravity_angle_range(X3, acc_mag),
      grav_X3 = gravity_angle(X3, X1, X2),
      grav_change_X3 = gravity_angle_change(X3, X1, X2),
      entropy_X3 = entropy(X3),

      # Magnitude
      mean_mag = mean(acc_mag),
      median_mag = median(acc_mag),
      sd_mag = sd(acc_mag),
      min_mag = min(acc_mag),
      max_mag = max(acc_mag),
      rms_mag = rms(acc_mag),
      q05_mag = quantile(acc_mag, 0.05),
      q25_mag = quantile(acc_mag, 0.25),
      q75_mag = quantile(acc_mag, 0.75),
      q95_mag = quantile(acc_mag, 0.95),
      skew_mag = e1071::skewness(acc_mag),
      kurtosis_mag = e1071::kurtosis(acc_mag),
      ar1_lag1_mag = lagged_cor(acc_mag, lag = 1),
      ar1_lag2_mag = lagged_cor(acc_mag, lag = 2),
      mean_abs_diff_mag = mean_abs_diff(acc_mag),
      max_abs_diff_mag = max_abs_diff(acc_mag),
      slope_mag = slope(acc_mag),
      diff_mag = diff_start_end(acc_mag),
      entropy_mag = entropy(acc_mag),

      # Relationships
      cc_lag1_X1X2 = lagged_cor(X1, X2, lag = 1),
      cc_lag1_X1X3 = lagged_cor(X1, X3, lag = 1),
      cc_lag1_X2X3 = lagged_cor(X2, X3, lag = 1),
      sma = sma(X1, X2, X3),
      mean_jerk = mean(jerk(X1, X2, X3)),
      sd_jerk = sd(jerk(X1, X2, X3)),
      rms_jerk = rms(jerk(X1, X2, X3)),
      dom = dominance(X1, X2, X3),
      diff_orient = diff_orientation(X1, X2, X3),
      angle_X1X2 = sigangle(X1, X2),
      angle_X1X3 = sigangle(X1, X3),
      angle_X2X3 = sigangle(X2, X3),

      # Keep track of epoch lengths (some epochs are less than 128 samples)
      n_samples = n()
    )
  return(time_domain_features)
}
