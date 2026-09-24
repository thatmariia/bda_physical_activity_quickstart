# ==========================================================
# == FUNCTIONS FOR CREATING GRAVITY FRAME & FEATURES
# ==========================================================

# ----------------------------------------------------------
# -- The gravity frame
# ----------------------------------------------------------

#' Describe the movement of an epoch relative to its own gravity direction, which
#' is the same description whatever way the phone happens to be held
#' @param joint_df A data frame from `join_sensors()`
#' @param n_samples_per_epoch The number of samples per epoch
#' @return A data frame with one row per sample and the four gravity frame channels:
#' `vert` (acceleration along gravity), `horiz` (acceleration in the horizontal plane),
#' `yaw` (rotation around the gravity axis) and `tilt` (rotation around a horizontal axis)
gravity_frame <- function(joint_df, n_samples_per_epoch = 128) {
  joint_df |>
    add_epoch(n_samples_per_epoch) |>
    # the direction of gravity is where the acceleration points on average
    mutate(
      gx = mean(acc_X1), gy = mean(acc_X2), gz = mean(acc_X3),
      .by = epoch
    ) |>
    mutate(
      g_norm = pmax(1e-9, magnitude(gx, gy, gz)),
      gx = gx / g_norm, gy = gy / g_norm, gz = gz / g_norm,
      vert = acc_X1 * gx + acc_X2 * gy + acc_X3 * gz,
      horiz = magnitude(acc_X1 - vert * gx, acc_X2 - vert * gy, acc_X3 - vert * gz),
      yaw = gyro_X1 * gx + gyro_X2 * gy + gyro_X3 * gz,
      tilt = magnitude(gyro_X1 - yaw * gx, gyro_X2 - yaw * gy, gyro_X3 - yaw * gz)
    ) |>
    select(epoch, vert, horiz, yaw, tilt)
}

#' Compute the spectrum of each gravity frame channel of each epoch
#' @param frame_df A data frame from `gravity_frame()`
#' @return A data frame with the spectral densities of the channels per epoch
gravity_spectrum_df <- function(frame_df, n_samples_per_epoch = 128, sample_rate = 50) {
  frame_df |>
    reframe(
      {
        x <- cbind(vert, horiz, yaw, tilt)
        n_missing <- n_samples_per_epoch - n()
        if (n_missing > 0) {
          x <- rbind(x, matrix(0, nrow = n_missing, ncol = ncol(x)))
        }

        sp <- spectrum(x, span = 15, plot = FALSE)
        tibble(
          freq = sp$freq,
          freq_hz = sp$freq * sample_rate,
          spec_vert = sp$spec[, 1],
          spec_horiz = sp$spec[, 2],
          spec_yaw = sp$spec[, 3],
          spec_tilt = sp$spec[, 4]
        )
      },
      .by = epoch
    )
}

# ----------------------------------------------------------
# -- Correlations
# ----------------------------------------------------------

#' Every pair of gravity frame channels, at delays of -2 to 2 samples
frame_lag_pairs <- function() {
  channels <- c("vert", "horiz", "yaw", "tilt")
  expand_grid(from = channels, to = channels, lag = -2:2) |>
    # the other order is the same correlation at the opposite lag, so it is left out
    filter(from < to) |>
    mutate(name = lag_name("cc", paste0(from, "_", to), lag))
}
