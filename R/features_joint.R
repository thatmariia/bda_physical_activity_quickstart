# ==========================================================
# == FUNCTIONS FOR JOINT (ACC + GYRO) FEATURES
# ==========================================================

#' Join acc and gyro signals of one experiment by sample
#' (from copilot assignment)
join_sensors <- function(acc_df, gyro_df) {
  inner_join(
    acc_df |> rename(acc_X1 = X1, acc_X2 = X2, acc_X3 = X3),
    gyro_df |> select(userid, trial, sampleid, gyro_X1 = X1, gyro_X2 = X2, gyro_X3 = X3),
    by = c("userid", "trial", "sampleid")
  )
}

# ----------------------------------------------------------
# -- Correlations
# ----------------------------------------------------------

#' Every pair of an acc and a gyro channel, at delays of -2 to 2 samples
cross_lag_pairs <- function() {
  expand_grid(
    from = c("acc_X1", "acc_X2", "acc_X3", "acc_dyn_mag"),
    to = c("gyro_X1", "gyro_X2", "gyro_X3", "gyro_mag"),
    lag = -2:2
  ) |>
    mutate(name = lag_name("cc", paste0(from, "_", to), lag))
}

# ----------------------------------------------------------
# -- Feature extraction
# ----------------------------------------------------------

#' Extract joint features from joined acc and gyro signals segmented into epochs
#' @param joint_df A data frame from `join_sensors()`
#' @param n_samples_per_epoch The number of samples per epoch (default is 128, corresponding to 2.56 seconds at 50 Hz)
#' @return A data frame with one row per epoch
get_joint_features <- function(joint_df, n_samples_per_epoch = 128) {
  lag_pairs <- cross_lag_pairs()

  joint_features <- joint_df |>
    add_epoch(n_samples_per_epoch) |>
    group_by(epoch) |>
    mutate(
      acc_dyn_mag = dynamic_magnitude(acc_X1, acc_X2, acc_X3),
      gyro_mag = magnitude(gyro_X1, gyro_X2, gyro_X3),
      # alignment of the acceleration and rotation directions
      dot = acc_X1 * gyro_X1 + acc_X2 * gyro_X2 + acc_X3 * gyro_X3
    ) |>
    summarise(
      mean_dot = mean(dot),
      sd_dot = sd(dot),
      gyro_acc_rms_ratio = rms(gyro_mag) / (rms(acc_dyn_mag) + 1e-6),
      grav_angle_gyro = vector_angle(
        mean_vector(gyro_X1, gyro_X2, gyro_X3),
        mean_vector(acc_X1, acc_X2, acc_X3)
      ),
      # Correlations (e.g. cc_acc_X1_gyro_X1_lag0)
      cors = lagged_cors(pick(everything()), lag_pairs),
      .groups = "drop"
    ) |>
    unpack(cors)
  return(joint_features)
}
