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

#' Compute the rate of local maxima in a signal
peak_rate <- function(x) {
  if (length(x) < 3) {
    return(0)
  }
  dx <- diff(x)
  n_peaks <- sum(dx[-length(dx)] > 0 & dx[-1] < 0)
  return(n_peaks / (length(x) - 2))
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

#' Compute the autocorrelations of a signal at lags 1 to max_lag
#' (from copilot assignment)
autocorrelations <- function(x, max_lag = 64) {
  max_lag <- min(max_lag, length(x) - 2)
  if (max_lag < 1 || sd(x) == 0) {
    return(numeric(0))
  }
  return(drop(acf(x, lag.max = max_lag, plot = FALSE)$acf)[-1])
}

#' Find the lag of the strongest autocorrelation peak of a signal
#' (from copilot assignment)
autocor_peak_lag <- function(x, min_lag = 5, max_lag = 64) {
  r <- autocorrelations(x, max_lag)
  if (length(r) < min_lag + 1) {
    return(0)
  }
  # Local maxima of the autocorrelation, ignoring very short lags
  lags <- seq(min_lag, length(r) - 1)
  is_peak <- r[lags] > r[lags - 1] & r[lags] >= r[lags + 1]
  if (!any(is_peak)) {
    return(0)
  }
  peak_lags <- lags[is_peak]
  return(peak_lags[which.max(r[peak_lags])])
}

#' Compute the height of the strongest autocorrelation peak of a signal
#' (from copilot assignment)
autocor_peak <- function(x, min_lag = 5, max_lag = 64) {
  lag <- autocor_peak_lag(x, min_lag, max_lag)
  if (lag == 0) {
    return(0)
  }
  return(autocorrelations(x, max_lag)[lag])
}

#' Compute when the largest deviation from the mean occurs, as a share of the epoch
#' (from copilot assignment)
excursion_time <- function(x) {
  if (length(x) < 2) {
    return(0)
  }
  return((which.max(abs(x - mean(x))) - 1) / (length(x) - 1))
}

#' Compute the difference between the means of the second and first half
#' (from copilot assignment)
half_mean_diff <- function(x) {
  half <- floor(length(x) / 2)
  if (half < 1) {
    return(0)
  }
  return(mean(tail(x, half)) - mean(head(x, half)))
}

#' Compute the difference between the movement (RMS around the mean)
#' of the second and first half
#' (from copilot assignment)
half_rms_diff <- function(x) {
  half <- floor(length(x) / 2)
  if (half < 2) {
    return(0)
  }
  rms_around_mean <- function(y) sqrt(mean((y - mean(y))^2))
  return(rms_around_mean(tail(x, half)) - rms_around_mean(head(x, half)))
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
      mag = sqrt(X1^2 + X2^2 + X3^2)
    ) |>
    # remove gravity per epoch
    mutate(
      dyn_mag = sqrt((X1 - mean(X1))^2 + (X2 - mean(X2))^2 + (X3 - mean(X3))^2),
      .by = epoch
    ) |>
    # extract features from each epoch
    group_by(epoch) |>
    summarise(
      # Activity label for the epoch = most common value (mode)
      aggr_activity = most_common_value(activity),
      # Confidence of the activity label = proportion of samples with that label
      activity_confidence = value_confidence(aggr_activity, activity),

      # Keep starting sample ID of the epoch
      sampleid = sampleid[1],

      # Features of each axis magnitude (e.g. mean_X1, mean_mag)
      across(
        c(X1, X2, X3, mag, dyn_mag),
        list(
          mean = mean,
          median = median,
          sd = sd,
          min = min,
          max = max,
          rms = rms,
          q05 = \(x) quantile(x, 0.05),
          q25 = \(x) quantile(x, 0.25),
          q75 = \(x) quantile(x, 0.75),
          q95 = \(x) quantile(x, 0.95),
          skew = e1071::skewness,
          kurtosis = e1071::kurtosis,
          ar1_lag1 = \(x) lagged_cor(x, lag = 1),
          ar1_lag2 = \(x) lagged_cor(x, lag = 2),
          mean_abs_diff = mean_abs_diff,
          max_abs_diff = max_abs_diff,
          slope = slope,
          diff = diff_start_end,
          peaks = peak_rate,
          entropy = entropy,
          autocor_peak = autocor_peak,
          autocor_peak_lag = autocor_peak_lag,
          excursion_time = excursion_time,
          half_mean_diff = half_mean_diff,
          half_rms_diff = half_rms_diff
        ),
        .names = "{.fn}_{.col}"
      ),

      # Features of each axis only
      across(
        c(X1, X2, X3),
        list(
          rot = \(x) rotation(x, sample_rate),
          rot_abs = \(x) abs_rotation(x, sample_rate),
          rot_asym = \(x) rotation_asym(x, sample_rate),
          zcr = zero_cross_rate,
          grav_range = \(x) gravity_angle_range(x, mag)
        ),
        .names = "{.fn}_{.col}"
      ),

      # Features comparing one axis with the other two
      var_rat_X1 = var_ratio(X1, X2, X3),
      var_rat_X2 = var_ratio(X2, X1, X3),
      var_rat_X3 = var_ratio(X3, X1, X2),
      grav_X1 = gravity_angle(X1, X2, X3),
      grav_X2 = gravity_angle(X2, X3, X1),
      grav_X3 = gravity_angle(X3, X1, X2),
      grav_change_X1 = gravity_angle_change(X1, X2, X3),
      grav_change_X2 = gravity_angle_change(X2, X3, X1, k = 0.1),
      grav_change_X3 = gravity_angle_change(X3, X1, X2),

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
