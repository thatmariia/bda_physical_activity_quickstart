# ==========================================================
# == FUNCTIONS FOR EXTRACTING GRAVITY FRAME FEATURES
# ==========================================================

#' Extract the features of the gravity frame channels, which only include the ones
#' that mean something for single signals (so no magnitudes or angles between axes)
#' @param joint_df A data frame from `join_sensors()`
#' @param n_samples_per_epoch The number of samples per epoch
#' @param sample_rate The sample rate of the signal
#' @return A data frame with one row per epoch
get_gravity_features <- function(joint_df, n_samples_per_epoch = 128, sample_rate = 50) {
  lag_pairs <- frame_lag_pairs()
  frame_df <- gravity_frame(joint_df, n_samples_per_epoch)
  channels <- c("vert", "horiz", "yaw", "tilt")

  time_features <- frame_df |>
    group_by(epoch) |>
    summarise(
      # Features of each channel (e.g. mean_vert, sd_yaw)
      across(
        all_of(channels),
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
          zcr = zero_cross_rate,
          # how much movement there is in total, and how one-sided it is
          rot = \(x) rotation(x, sample_rate),
          rot_abs = \(x) rotation_abs(x, sample_rate),
          rot_asym = \(x) rotation_asym(x, sample_rate)
        ),
        .names = "{.fn}_{.col}"
      ),

      # How the movement divides between the channels
      var_share_vert = var(vert) / (var(vert) + var(horiz)),
      var_share_yaw = var(yaw) / (var(yaw) + var(tilt)),
      dom_acc = max(sd(vert), sd(horiz)) / (sd(vert) + sd(horiz)),
      dom_rot = max(sd(yaw), sd(tilt)) / (sd(yaw) + sd(tilt)),
      sma_acc = mean(abs(vert) + abs(horiz)),
      sma_rot = mean(abs(yaw) + abs(tilt)),

      # Rate of change of the acceleration and of the rotation
      mean_jerk_acc = mean(magnitude(diff(vert), diff(horiz), 0)),
      sd_jerk_acc = sd(magnitude(diff(vert), diff(horiz), 0)),
      rms_jerk_acc = rms(magnitude(diff(vert), diff(horiz), 0)),
      mean_jerk_rot = mean(magnitude(diff(yaw), diff(tilt), 0)),
      sd_jerk_rot = sd(magnitude(diff(yaw), diff(tilt), 0)),
      rms_jerk_rot = rms(magnitude(diff(yaw), diff(tilt), 0)),

      # Autoregression coefficients and autocorrelation peaks (e.g. ar1_vert)
      rhythm = rhythm_features(list(vert = vert, horiz = horiz, yaw = yaw, tilt = tilt)),

      # Correlations between the channels (e.g. cc_vert_yaw_lag0)
      cors = lagged_cors(pick(everything()), lag_pairs)
    ) |>
    unpack(c(rhythm, cors))

  freq_features <- frame_df |>
    gravity_spectrum_df(n_samples_per_epoch, sample_rate) |>
    get_frequency_domain_features()

  left_join(time_features, freq_features, by = "epoch")
}
