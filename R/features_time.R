# ==========================================================
# == FUNCTIONS FOR EXTRACTING TIME DOMAIN FEATURES
# ==========================================================

#' Extract all time domain features from a signal data frame segmented into epochs
#' @param signal_df A data frame containing the signal data with columns: userid, trial, sampleid, X1, X2, X3, activity
#' @param n_samples_per_epoch The number of samples per epoch (default is 128, corresponding to 2.56 seconds at 50 Hz)
#' @param sample_rate The sample rate of the signal (default is 50 Hz)
#' @return A data frame containing the extracted time domain features for each epoch
get_time_domain_features <- function(signal_df, n_samples_per_epoch = 128, sample_rate = 50) {
  lag_pairs <- time_lag_pairs()

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
          peaks = peak_rate
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

      # Autoregression coefficients and autocorrelation peaks (e.g. ar1_X1, acf_peak_X1)
      rhythm = rhythm_features(list(X1 = X1, X2 = X2, X3 = X3, mag = mag, dyn_mag = dyn_mag)),

      # Correlations (e.g. acf_X1_lag1, cc_X1_X2_lag1)
      cors = lagged_cors(pick(everything()), lag_pairs),

      # Keep track of epoch lengths (some epochs are less than 128 samples)
      n_samples = n()
    ) |>
    unpack(c(rhythm, cors))
  return(time_domain_features)
}
