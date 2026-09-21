# ==========================================================
# == FUNCTIONS FOR TIME DOMAIN FEATURES
# ==========================================================

# ----------------------------------------------------------
# -- Epochs and magnitudes
# ----------------------------------------------------------

#' Add the epoch each sample belongs to
add_epoch <- function(df, n_samples_per_epoch = 128) {
  df |> mutate(epoch = sampleid %/% n_samples_per_epoch)
}

#' Compute the magnitude of a three-axial signal
magnitude <- function(x, y, z) {
  sqrt(x^2 + y^2 + z^2)
}

#' Compute the magnitude of a three-axial signal with its mean (gravity) removed
dynamic_magnitude <- function(x, y, z) {
  magnitude(x - mean(x), y - mean(y), z - mean(z))
}

# ----------------------------------------------------------
# -- Basic statistics
# ----------------------------------------------------------

#' Compute the root mean square of a signal
rms <- function(x) {
  sqrt(mean(x^2))
}

#' Compute the Shannon entropy of a distribution given by (unnormalised) weights
shannon_entropy <- function(p) {
  p <- p[p > 0] / sum(p)
  return(-sum(p * log(p)))
}

#' Compute the entropy of the distribution of signal values
entropy <- function(x, bins = 10) {
  if (length(unique(x)) < 2) {
    return(0)
  }
  return(shannon_entropy(hist(x, breaks = bins, plot = FALSE)$counts))
}

# ----------------------------------------------------------
# -- Change and trend
# ----------------------------------------------------------

#' Compute the mean absolute first difference of a signal
mean_abs_diff <- function(x) {
  mean(abs(diff(x)))
}

#' Compute the maximum absolute first difference of a signal
max_abs_diff <- function(x) {
  max(abs(diff(x)))
}

#' Compute the standard deviation of the jerk (rate of change) of a signal
#' (from Reyes-Ortiz et al., 2015)
jerk_sd <- function(x, sample_rate = 50) {
  if (length(x) < 3) {
    return(0)
  }
  return(sd(diff(x)) * sample_rate)
}

#' Compute the jerk of a three-axial signal
jerk <- function(x, y, z) {
  magnitude(diff(x), diff(y), diff(z))
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

#' Compute the difference between the means of the last and first k share of a signal
diff_start_end <- function(x, k = 0.1) {
  if (length(x) < 2) {
    return(0)
  }
  m <- max(1, floor(length(x) * k))
  return(mean(tail(x, m)) - mean(head(x, m)))
}

#' Compute the difference between the movement (RMS around the mean)
#' of the last and first k share of a signal
#' (from copilot assignment)
spread_diff_start_end <- function(x, k = 0.1) {
  m <- max(1, floor(length(x) * k))
  if (m < 2) {
    return(0)
  }
  rms_around_mean <- function(y) rms(y - mean(y))
  return(rms_around_mean(tail(x, m)) - rms_around_mean(head(x, m)))
}

#' Compute when the largest deviation from the mean occurs, as a share of the epoch
#' (from copilot assignment)
excursion_time <- function(x) {
  if (length(x) < 2) {
    return(0)
  }
  return((which.max(abs(x - mean(x))) - 1) / (length(x) - 1))
}

# ----------------------------------------------------------
# -- Rhythm
# ----------------------------------------------------------

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
acf_peak_lag <- function(x, min_lag = 5, max_lag = 64) {
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
acf_peak <- function(x, min_lag = 5, max_lag = 64) {
  lag <- acf_peak_lag(x, min_lag, max_lag)
  if (lag == 0) {
    return(0)
  }
  return(autocorrelations(x, max_lag)[lag])
}

#' Compute one autoregression coefficient of a signal (Burg method)
#' (from Reyes-Ortiz et al., 2015)
ar_coefficient <- function(x, k, order = 4) {
  if (length(x) <= 2 * order || sd(x) == 0) {
    return(0)
  }
  fit <- ar.burg(x, aic = FALSE, order.max = order, demean = TRUE)
  return(fit$ar[k])
}

# ----------------------------------------------------------
# -- Rotation
# ----------------------------------------------------------

#' Compute the integrated rotation of a signal
rotation <- function(x, fs) {
  sum(x) / fs
}

#' Compute the absolute integrated rotation of a signal
rotation_abs <- function(x, fs) {
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

# ----------------------------------------------------------
# -- Orientation and gravity
# ----------------------------------------------------------

#' Compute the mean vector of a three-axial signal
mean_vector <- function(x, y, z) {
  c(mean(x), mean(y), mean(z))
}

#' Compute the angle between two vectors in degrees, 0 if either has no length
vector_angle <- function(a, b) {
  norm_a <- sqrt(sum(a^2))
  norm_b <- sqrt(sum(b^2))
  if (norm_a == 0 || norm_b == 0) {
    return(0)
  }
  return(acos(pmin(1, pmax(-1, sum(a * b) / (norm_a * norm_b)))) * 180 / pi)
}

#' Compute the angle between the first axis and gravity (the mean vector)
gravity_angle <- function(x, y, z) {
  vector_angle(mean_vector(x, y, z), c(1, 0, 0))
}

#' Compute the range of the gravity angle
gravity_angle_range <- function(x, mag) {
  angles <- acos(pmin(1, pmax(-1, x / mag))) * 180 / pi
  if (all(is.na(angles))) {
    return(0)
  }
  return(max(angles) - min(angles))
}

#' Compute the change in gravity angle between the start and end of a three-axial signal
gravity_angle_change <- function(x, y, z, k = 0.1) {
  if (length(x) < 2) {
    return(0)
  }
  m <- max(1, floor(length(x) * k))
  start <- gravity_angle(head(x, m), head(y, m), head(z, m))
  end <- gravity_angle(tail(x, m), tail(y, m), tail(z, m))
  return(end - start)
}

#' Compute the difference in orientation between the start and end of a three-axial signal
diff_orientation <- function(x, y, z) {
  diffs <- c(
    diff_start_end(x),
    diff_start_end(y),
    diff_start_end(z)
  )
  return(sqrt(sum(diffs^2)))
}

# ----------------------------------------------------------
# -- Relationships between axes
# ----------------------------------------------------------

#' Compute the magnitude area of a three-axial signal
sma <- function(x, y, z) {
  mean(abs(x) + abs(y) + abs(z))
}

#' Compute the share of the variance that is along the first axis
var_share <- function(x, y, z) {
  var(x) / (var(x) + var(y) + var(z))
}

#' Compute the share of the strongest axis in the overall movement
axis_dominance <- function(x, y, z) {
  max(c(sd(x), sd(y), sd(z))) / sum(c(sd(x), sd(y), sd(z)))
}

#' Compute the angle of the mean signal in the plane of two axes
plane_angle <- function(x, y) {
  atan2(mean(x), mean(y)) * 180 / pi
}

# ----------------------------------------------------------
# -- Correlations
# ----------------------------------------------------------

#' Compute the correlation between a signal and another one delayed by `lag` samples,
#' over the part of the epoch where they overlap, 0 if either is constant
#' (a negative lag delays the first signal instead)
lagged_cor <- function(x, y = x, lag = 0) {
  if (lag < 0) {
    return(lagged_cor(y, x, -lag))
  }
  n <- length(x)
  if (n - lag < 2) {
    return(0)
  }
  x_now <- x[seq.int(lag + 1, n)]
  y_before <- y[seq_len(n - lag)]
  if (sd(x_now) == 0 || sd(y_before) == 0) {
    return(0)
  }
  return(cor(x_now, y_before))
}

#' Name a lagged correlation, e.g. cc_X1_X2_lag1 (a negative lag becomes lagm1)
lag_name <- function(prefix, signals, lag) {
  paste0(prefix, "_", signals, "_lag", ifelse(lag < 0, paste0("m", -lag), lag))
}

#' Signals to correlate with a lagged copy of themselves (autocorrelations)
#' or of another axis (cross-correlations)
time_lag_pairs <- function() {
  bind_rows(
    expand_grid(from = c("X1", "X2", "X3", "mag", "dyn_mag"), lag = 1:2) |>
      mutate(to = from, name = lag_name("acf", from, lag)),
    # the other order is the same correlation at the opposite lag, so it is left out
    expand_grid(from = c("X1", "X2", "X3"), to = c("X1", "X2", "X3"), lag = -2:2) |>
      filter(from < to) |>
      mutate(name = lag_name("cc", paste0(from, "_", to), lag))
  )
}

#' Compute the lagged correlation of every pair in `pairs`, as a one-row data frame
#' (from copilot assignment)
#' @param signals A data frame with the signals of one epoch
#' @param pairs A data frame with the columns from, to, lag and name
lagged_cors <- function(signals, pairs = time_lag_pairs()) {
  values <- pmap_dbl(
    pairs,
    \(from, to, lag, name) lagged_cor(signals[[from]], signals[[to]], lag)
  )
  names(values) <- pairs$name
  return(as_tibble_row(values))
}

# ----------------------------------------------------------
# -- Feature extraction
# ----------------------------------------------------------

#' Extract all time domain features from a signal data frame segmented into epochs
#' @param signal_df A data frame containing the signal data with columns: userid, trial, sampleid, X1, X2, X3, activity
#' @param n_samples_per_epoch The number of samples per epoch (default is 128, corresponding to 2.56 seconds at 50 Hz)
#' @param sample_rate The sample rate of the signal (default is 50 Hz)
#' @return A data frame containing the extracted time domain features for each epoch
get_time_domain_features <- function(signal_df, n_samples_per_epoch = 128, sample_rate = 50) {
  time_domain_features <- signal_df |>
    # partition into epochs and add the magnitude of the signal
    add_epoch(n_samples_per_epoch) |>
    mutate(mag = magnitude(X1, X2, X3)) |>
    # remove gravity per epoch
    mutate(dyn_mag = dynamic_magnitude(X1, X2, X3), .by = epoch) |>
    # extract features from each epoch
    group_by(epoch) |>
    summarise(
      # Activity label for the epoch = most common value (mode)
      aggr_activity = most_common_value(activity),
      # Confidence of the activity label = proportion of samples with that label
      activity_confidence = value_confidence(aggr_activity, activity),

      # Keep starting sample ID of the epoch
      sampleid = sampleid[1],

      # Features of each axis and magnitude (e.g. mean_X1, mean_mag)
      across(
        c(X1, X2, X3, mag, dyn_mag),
        list(
          # basic statistics
          mean = mean,
          median = median,
          sd = sd,
          mad = mad,
          min = min,
          max = max,
          q05 = \(x) quantile(x, 0.05),
          q25 = \(x) quantile(x, 0.25),
          q75 = \(x) quantile(x, 0.75),
          q95 = \(x) quantile(x, 0.95),
          rms = rms,
          skew = e1071::skewness,
          kurt = e1071::kurtosis,
          entropy = entropy,
          # change and trend
          mean_abs_diff = mean_abs_diff,
          max_abs_diff = max_abs_diff,
          jerk_sd = \(x) jerk_sd(x, sample_rate),
          slope = slope,
          diff = diff_start_end,
          half_mean_diff = \(x) diff_start_end(x, k = 0.5),
          half_rms_diff = \(x) spread_diff_start_end(x, k = 0.5),
          excursion_time = excursion_time,
          # rhythm
          peaks = peak_rate,
          acf_peak = acf_peak,
          acf_peak_lag = acf_peak_lag,
          ar1 = \(x) ar_coefficient(x, 1),
          ar2 = \(x) ar_coefficient(x, 2),
          ar3 = \(x) ar_coefficient(x, 3),
          ar4 = \(x) ar_coefficient(x, 4)
        ),
        .names = "{.fn}_{.col}"
      ),

      # Features of each axis only
      across(
        c(X1, X2, X3),
        list(
          zcr = zero_cross_rate,
          rot = \(x) rotation(x, sample_rate),
          rot_abs = \(x) rotation_abs(x, sample_rate),
          rot_asym = \(x) rotation_asym(x, sample_rate),
          grav_angle_range = \(x) gravity_angle_range(x, mag)
        ),
        .names = "{.fn}_{.col}"
      ),

      # Features comparing one axis with the other two
      var_share_X1 = var_share(X1, X2, X3),
      var_share_X2 = var_share(X2, X1, X3),
      var_share_X3 = var_share(X3, X1, X2),
      grav_angle_X1 = gravity_angle(X1, X2, X3),
      grav_angle_X2 = gravity_angle(X2, X3, X1),
      grav_angle_X3 = gravity_angle(X3, X1, X2),
      grav_angle_change_X1 = gravity_angle_change(X1, X2, X3),
      grav_angle_change_X2 = gravity_angle_change(X2, X3, X1),
      grav_angle_change_X3 = gravity_angle_change(X3, X1, X2),

      # Features of all three axes together
      sma = sma(X1, X2, X3),
      mean_jerk = mean(jerk(X1, X2, X3)),
      sd_jerk = sd(jerk(X1, X2, X3)),
      rms_jerk = rms(jerk(X1, X2, X3)),
      axis_dom = axis_dominance(X1, X2, X3),
      diff_orient = diff_orientation(X1, X2, X3),
      axis_angle_X1_X2 = plane_angle(X1, X2),
      axis_angle_X1_X3 = plane_angle(X1, X3),
      axis_angle_X2_X3 = plane_angle(X2, X3),

      # Correlations (e.g. acf_X1_lag1, cc_X1_X2_lag1)
      cors = lagged_cors(pick(everything())),

      # Keep track of epoch lengths (some epochs are less than 128 samples)
      n_samples = n()
    ) |>
    unpack(cors)
  return(time_domain_features)
}
