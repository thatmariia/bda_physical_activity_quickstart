# ==========================================================
# == FUNCTIONS FOR FREQUENCY DOMAIN FEATURES
# ==========================================================

#' Compute the mean frequency of a signal
mean_frequency <- function(freq, spec) {
  delta_f <- freq[2] - freq[1] # gap size between frequencies
  normalizing_constant <- sum(spec * delta_f) # ≈ ∫S(f)df
  mean_freq <- sum(freq * spec * delta_f) / normalizing_constant
  return(mean_freq)
}

#' Compute the standard deviation of the frequency of a signal
sd_frequency <- function(freq, spec) {
  mean_freq <- mean_frequency(freq, spec)
  delta_freq <- freq[2] - freq[1]
  normalizing_constant <- sum(spec * delta_freq)
  if (normalizing_constant == 0) {
    return(0)
  }
  ds <- (freq - mean_freq)^2
  var_freq <- sum(ds * spec * delta_freq) / normalizing_constant
  return(sqrt(var_freq))
}

#' Compute the skewness of the frequency of a signal
skew_frequency <- function(freq, spec) {
  mean_freq <- mean_frequency(freq, spec)
  sd_freq <- sd_frequency(freq, spec)
  delta_freq <- freq[2] - freq[1]
  normalizing_constant <- sum(spec * delta_freq)
  if (normalizing_constant == 0 || sd_freq == 0) {
    return(0)
  }
  dc <- (freq - mean_freq)^3
  skew <- sum(dc * spec * delta_freq) / (normalizing_constant * sd_freq^3)
  return(skew)
}

#' Compute the kurtosis of the frequency of a signal
kurtosis_frequency <- function(freq, spec) {
  mean_freq <- mean_frequency(freq, spec)
  sd_freq <- sd_frequency(freq, spec)
  delta_freq <- freq[2] - freq[1]
  normalizing_constant <- sum(spec * delta_freq)
  if (normalizing_constant == 0 || sd_freq == 0) {
    return(0)
  }
  df <- (freq - mean_freq)^4
  kurt <- sum(df * spec * delta_freq) / (normalizing_constant * sd_freq^4)
  return(kurt)
}

#' Compute the total power of a signal
spectral_power <- function(freq, spec) {
  df <- mean(diff(freq))
  sum(spec * df)
}

#' Compute the peak power of a signal
peak_power <- function(freq, spec) {
  max(spec[freq > 0])
}

#' Compute the spectral entropy of a signal
spectral_entropy <- function(spec) {
  if (length(spec) == 0 || sum(spec) == 0) {
    return(0)
  }
  spec <- spec[spec > 0]
  spec <- spec / sum(spec)
  return(-sum(spec * log2(spec)) / log2(length(spec)))
}

#' Compute the spectral edge frequency of a signal
spectral_edge_frequency <- function(freq, spec, k = 0.9) {
  tot <- sum(spec)
  if (tot == 0) {
    return(0)
  }
  cumsum_spec <- cumsum(spec)
  threshold <- k * tot
  edge_freq <- which(cumsum_spec >= threshold)[1]
  return(freq[edge_freq])
}

#' Compute the dominant frequency of a signal in Hz
dominant_frequency <- function(freq, spec, sample_rate = 50, min_hz = 0.5) {
  freq_hz <- freq * sample_rate
  keep <- freq_hz >= min_hz
  if (!any(keep) || sum(spec[keep]) == 0) {
    return(0)
  }
  freq_hz[keep][which.max(spec[keep])]
}

dominant_power_ratio <- function(freq, spec, sample_rate = 50, min_hz = 0.5) {
  dom_freq <- dominant_frequency(freq, spec, sample_rate, min_hz)
  if (dom_freq == 0) {
    return(0)
  }
  return(spec[which.min(abs(freq * sample_rate - dom_freq))] / sum(spec))
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

#' Extract frequency domain features from a spectrum data frame,
#' which is assumed to be already segmented into epochs
#' @param spectrum_df A data frame containing the spectrum data with columns: epoch, freq, spec1, spec2, ...
#' @return A data frame containing the extracted frequency domain features for each epoch
get_frequency_domain_features <- function(spectrum_df) {
  userfreqdom <- spectrum_df %>%
    group_by(epoch) %>%
    summarise(
      # Features of each signal (e.g. mean_freq1 for spec1)
      across(
        c(spec1, spec2, spec3),
        list(
          dom_power_ratio = \(s) dominant_power_ratio(freq, s),
          dom_freq = \(s) dominant_frequency(freq, s),
          mean_freq = \(s) mean_frequency(freq, s),
          sd_freq = \(s) sd_frequency(freq, s),
          skew_freq = \(s) skew_frequency(freq, s),
          kurt_freq = \(s) kurtosis_frequency(freq, s),
          entropy = spectral_entropy,
          edge_freq = \(s) spectral_edge_frequency(freq, s),
          bin1_ = \(s) band_bin_share(s, 1),
          bin2_ = \(s) band_bin_share(s, 2),
          bin3_ = \(s) band_bin_share(s, 3),
          bin4_ = \(s) band_bin_share(s, 4),
          bin5_ = \(s) band_bin_share(s, 5),
          bin6_ = \(s) band_bin_share(s, 6),
          bin7_ = \(s) band_bin_share(s, 7),
          bin8_ = \(s) band_bin_share(s, 8)
        ),
        .names = "{.fn}{sub('spec', '', .col)}"
      )
    )
}
