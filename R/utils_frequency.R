#' Apply Gaussian smoothing to a numeric vector
#' @param x A numeric vector to be smoothed
#' @param sigma The standard deviation of the Gaussian kernel (default is 2)
#' @return A numeric vector containing the smoothed values
gaussian_smooth <- function(x, sigma = 2) {
  radius <- ceiling(3 * sigma)
  positions <- -radius:radius

  weights <- exp(-(positions^2) / (2 * sigma^2))
  weights <- weights / sum(weights)

  smooth_x <- stats::filter(x, filter = weights, sides = 2)
  return(smooth_x)
}

#   apply(
#     x,
#     2,
#     \(column) stats::filter(
#       column,
#       filter = weights,
#       sides = 2
#     )
#   )
# }

#' Compute the spectrum of a signal and optionally apply Gaussian smoothing
#' @param signal A numeric vector representing the signal
#' @param gaussian_sigma The standard deviation of the Gaussian kernel for smoothing (default is 3; set to 0 to disable smoothing)
#' @return A data frame containing the frequency bins and corresponding spectral densities of the signal
compute_spectrum <- function(signal, gaussian_sigma = 3) {
  # Compute the spectrum of the signal
  spec <- spectrum(signal, plot = FALSE)
  spec_df <- data.frame(freq = spec$freq, spec = spec$spec)
  if (gaussian_sigma <= 0) {
    return(spec_df)
  }
  spec_df$spec <- gaussian_smooth(spec$spec, sigma = gaussian_sigma)
  return(spec_df)
}

#' Convert a signal data frame into a segmented frequency domain representation
#' @param signal_df A data frame containing the signal data
#' @param n_samples_per_epoch The number of samples per epoch (default is 128, corresponding to 2.56 seconds at 50 Hz)
#' @param sample_rate The sample rate of the signal in Hz (default is 50 Hz)
#' @return A data frame containing the frequency domain representation for each epoch,
#' including the frequency bins and corresponding spectral densities for each signal.
convert_signal_to_spectrum_df <- function(signal_df, n_samples_per_epoch = 128, sample_rate = 50) {
  spectrum_df <- signal_df |>
    mutate(
      epoch = sampleid %/% n_samples_per_epoch,
      position = sampleid %% n_samples_per_epoch
    ) |>
    reframe(
      {
        # ==> START LLM https://chatgpt.com/share/6aab570f-9e38-83ed-a64e-fcd83ae13599
        n_missing <- n_samples_per_epoch - n()
        x <- cbind(X1, X2, X3)
        if (n_missing > 0) {
          x <- rbind(x, matrix(0, nrow = n_missing, ncol = 3))
        }
        # ==> END LLM

        sp <- spectrum(x, span = 15, plot = FALSE)
        tibble(
          freq = sp$freq, # cycles per sample
          freq_hz = sp$freq * sample_rate, # cycles per second (Hz)
          spec1 = sp$spec[, 1],
          spec2 = sp$spec[, 2],
          spec3 = sp$spec[, 3]
        )
      },
      .by = epoch
    )
  return(spectrum_df)
}
