# ==========================================================
# == FUNCTIONS FOR FREQUENCY DOMAIN FEATURES
# ==========================================================

# ----------------------------------------------------------
# -- Shape of the spectrum
# ----------------------------------------------------------

#' Compute a moment of the frequency distribution described by the spectrum
#' (k = 1 gives the mean frequency, k > 1 the central moment of that order)
spectral_moment <- function(freq, spec, k) {
  total <- sum(spec)
  if (total == 0) {
    return(0)
  }
  # the spectrum as a probability distribution over the frequencies
  p <- spec / total
  mean_freq <- sum(freq * p)
  if (k == 1) {
    return(mean_freq)
  }
  return(sum((freq - mean_freq)^k * p))
}

#' Compute the mean frequency of a signal
mean_frequency <- function(freq, spec) {
  spectral_moment(freq, spec, 1)
}

#' Compute the standard deviation of the frequency of a signal
sd_frequency <- function(freq, spec) {
  sqrt(spectral_moment(freq, spec, 2))
}

#' Compute the skewness of the frequency of a signal
skew_frequency <- function(freq, spec) {
  sd_freq <- sd_frequency(freq, spec)
  if (sd_freq == 0) {
    return(0)
  }
  return(spectral_moment(freq, spec, 3) / sd_freq^3)
}

#' Compute the kurtosis of the frequency of a signal
kurtosis_frequency <- function(freq, spec) {
  sd_freq <- sd_frequency(freq, spec)
  if (sd_freq == 0) {
    return(0)
  }
  return(spectral_moment(freq, spec, 4) / sd_freq^4)
}

#' Compute the spectral entropy of a signal
spectral_entropy <- function(spec) {
  if (length(spec) == 0 || sum(spec) == 0) {
    return(0)
  }
  # normalised by the entropy of a flat spectrum, so it does not depend on its length
  return(shannon_entropy(spec) / log(sum(spec > 0)))
}

#' Compute the spectral edge frequency of a signal
spectral_edge_frequency <- function(freq, spec, k = 0.9) {
  total <- sum(spec)
  if (total == 0) {
    return(0)
  }
  return(freq[which(cumsum(spec) >= k * total)[1]])
}

#' Compute the share of the spectral power in one of `n_bins` equal frequency bins
#' (from Reyes-Ortiz et al., 2015)
band_bin_share <- function(spec, bin, n_bins = 8) {
  total <- sum(spec)
  if (total == 0) {
    return(0)
  }
  bin_id <- ceiling(seq_along(spec) / length(spec) * n_bins)
  return(sum(spec[bin_id == bin]) / total)
}

# ----------------------------------------------------------
# -- Dominant frequency
# ----------------------------------------------------------

#' Find the position of the strongest frequency of a signal above `min_hz`, 0 if there is none
dominant_index <- function(freq, spec, sample_rate = 50, min_hz = 0.5) {
  keep <- which(freq * sample_rate >= min_hz)
  if (length(keep) == 0 || sum(spec[keep]) == 0) {
    return(0)
  }
  return(keep[which.max(spec[keep])])
}

#' Compute the dominant frequency of a signal in Hz
dominant_frequency <- function(freq, spec, sample_rate = 50, min_hz = 0.5) {
  i <- dominant_index(freq, spec, sample_rate, min_hz)
  if (i == 0) {
    return(0)
  }
  return(freq[i] * sample_rate)
}

#' Compute the share of the total power that is at the dominant frequency
dominant_power_ratio <- function(freq, spec, sample_rate = 50, min_hz = 0.5) {
  i <- dominant_index(freq, spec, sample_rate, min_hz)
  if (i == 0) {
    return(0)
  }
  return(spec[i] / sum(spec))
}

# ----------------------------------------------------------
# -- Feature extraction
# ----------------------------------------------------------

#' Extract frequency domain features from a spectrum data frame,
#' which is assumed to be already segmented into epochs
#' @param spectrum_df A data frame containing the spectrum data with columns: epoch, freq, spec_X1, spec_X2, ...
#' @return A data frame containing the extracted frequency domain features for each epoch
get_frequency_domain_features <- function(spectrum_df) {
  frequency_domain_features <- spectrum_df |>
    group_by(epoch) |>
    summarise(
      across(
        starts_with("spec_"),
        list(
          # shape of the spectrum
          mean_freq = \(s) mean_frequency(freq, s),
          sd_freq = \(s) sd_frequency(freq, s),
          skew_freq = \(s) skew_frequency(freq, s),
          kurt_freq = \(s) kurtosis_frequency(freq, s),
          spec_entropy = spectral_entropy,
          edge_freq = \(s) spectral_edge_frequency(freq, s),
          bin1 = \(s) band_bin_share(s, 1),
          bin2 = \(s) band_bin_share(s, 2),
          bin3 = \(s) band_bin_share(s, 3),
          bin4 = \(s) band_bin_share(s, 4),
          bin5 = \(s) band_bin_share(s, 5),
          bin6 = \(s) band_bin_share(s, 6),
          bin7 = \(s) band_bin_share(s, 7),
          bin8 = \(s) band_bin_share(s, 8),
          # dominant frequency
          dom_freq = \(s) dominant_frequency(freq, s),
          dom_power_ratio = \(s) dominant_power_ratio(freq, s)
        ),
        .names = "{.fn}_{sub('spec_', '', .col)}"
      )
    )
  return(frequency_domain_features)
}
