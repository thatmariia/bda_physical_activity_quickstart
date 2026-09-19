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

#' Compute the correlation between two signals, 0 if either is constant
#' (from copilot assignment)
safe_cor <- function(x, y) {
  if (length(x) < 2 || sd(x) == 0 || sd(y) == 0) {
    return(0)
  }
  return(cor(x, y))
}

#' Compute the angle between the mean vectors of two three-axial signals
#' (from Reyes-Ortiz et al., 2015)
mean_vector_angle <- function(ax, ay, az, bx, by, bz) {
  a <- c(mean(ax), mean(ay), mean(az))
  b <- c(mean(bx), mean(by), mean(bz))
  norm_a <- sqrt(sum(a^2))
  norm_b <- sqrt(sum(b^2))
  if (norm_a == 0 || norm_b == 0) {
    return(0)
  }
  return(acos(pmin(1, pmax(-1, sum(a * b) / (norm_a * norm_b)))) * 180 / pi)
}

#' Extract joint features from joined acc and gyro signals segmented into epochs
#' @param joint_df A data frame from `join_sensors()`
#' @param n_samples_per_epoch The number of samples per epoch (default is 128, corresponding to 2.56 seconds at 50 Hz)
#' @return A data frame with one row per epoch
get_joint_features <- function(joint_df, n_samples_per_epoch = 128) {
  joint_features <- joint_df |>
    # partition into epochs
    mutate(epoch = sampleid %/% n_samples_per_epoch) |>
    group_by(epoch) |>
    # remove gravity from the acceleration
    mutate(
      acc_dyn_mag = sqrt(
        (acc_X1 - mean(acc_X1))^2 + (acc_X2 - mean(acc_X2))^2 + (acc_X3 - mean(acc_X3))^2
      ),
      gyro_mag = sqrt(gyro_X1^2 + gyro_X2^2 + gyro_X3^2),
      # alignment of the acceleration and rotation directions
      dot = acc_X1 * gyro_X1 + acc_X2 * gyro_X2 + acc_X3 * gyro_X3
    ) |>
    summarise(
      cor_acc_gyro_mag = safe_cor(acc_dyn_mag, gyro_mag),
      mean_dot = mean(dot),
      sd_dot = sd(dot),
      gyro_acc_rms_ratio = rms(gyro_mag) / (rms(acc_dyn_mag) + 1e-6),
      cor_acc_gyro_X1 = safe_cor(acc_X1, gyro_X1),
      cor_acc_gyro_X2 = safe_cor(acc_X2, gyro_X2),
      cor_acc_gyro_X3 = safe_cor(acc_X3, gyro_X3),
      gyro_gravity_angle = mean_vector_angle(
        gyro_X1, gyro_X2, gyro_X3, acc_X1, acc_X2, acc_X3
      ),
      .groups = "drop"
    )
  return(joint_features)
}

