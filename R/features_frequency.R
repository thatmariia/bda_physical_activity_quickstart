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

#' Compute the spectral bandwidth of a signal
spectral_bandwidth <- function(freq, spec) {
  mean_freq <- mean_frequency(freq, spec)
  tot <- sum(spec)
  if (tot == 0) {
    return(0)
  }
  ds <- (freq - mean_freq)^2
  return(sqrt(sum(ds * spec) / tot))
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
dominant_frequency_hz <- function(freq, spec, sample_rate = 50, min_hz = 0.5) {
  freq_hz <- freq * sample_rate
  keep <- freq_hz >= min_hz
  if (!any(keep) || sum(spec[keep]) == 0) {
    return(0)
  }
  freq_hz[keep][which.max(spec[keep])]
}

dominant_power_ratio <- function(freq, spec, sample_rate = 50, min_hz = 0.5) {
  dom_freq <- dominant_frequency_hz(freq, spec, sample_rate, min_hz)
  if (dom_freq == 0) {
    return(0)
  }
  return(spec[which.min(abs(freq * sample_rate - dom_freq))] / sum(spec))
}

#' Compute the ratio of power in a specific frequency band to the total power
band_power_ratio <- function(freq, spec, lb, ub, sample_rate = 50) {
  freq_hz <- freq * sample_rate
  keep <- freq_hz >= lb & freq_hz <= ub
  if (!any(keep) || sum(spec[keep]) == 0) {
    return(0)
  }
  return(sum(spec[keep]) / sum(spec))
}

#' Extract frequency domain features from a spectrum data frame,
#' which is assumed to be already segmented into epochs
#' @param spectrum_df A data frame containing the spectrum data with columns: epoch, freq, spec1, spec2, ...
#' @return A data frame containing the extracted frequency domain features for each epoch
get_frequency_domain_features <- function(spectrum_df) {
  userfreqdom <- spectrum_df %>%
    group_by(epoch) %>%
    summarise(
      # Signal 1
      dom_power_ratio1 = dominant_power_ratio(freq, spec1),
      dom_freq1 = freq[which.max(spec1)],
      dom_hz_freq1 = dominant_frequency_hz(freq, spec1),
      mean_freq1 = mean_frequency(freq, spec1),
      sd_freq1 = sd_frequency(freq, spec1),
      skew_freq1 = skew_frequency(freq, spec1),
      kurt_freq1 = kurtosis_frequency(freq, spec1),
      entropy1 = spectral_entropy(spec1),
      bandwidth1 = spectral_bandwidth(freq, spec1),
      edge_freq1 = spectral_edge_frequency(freq, spec1),
      i1_band1 = band_power_ratio(freq, spec1, 0, 0.5),
      i2_band1 = band_power_ratio(freq, spec1, 0.5, 3),
      i3_band1 = band_power_ratio(freq, spec1, 3, 10),
      i4_band1 = band_power_ratio(freq, spec1, 10, 25),

      # Signal 2
      dom_power_ratio2 = dominant_power_ratio(freq, spec2),
      dom_freq2 = freq[which.max(spec2)],
      dom_hz_freq2 = dominant_frequency_hz(freq, spec2),
      mean_freq2 = mean_frequency(freq, spec2),
      sd_freq2 = sd_frequency(freq, spec2),
      skew_freq2 = skew_frequency(freq, spec2),
      kurt_freq2 = kurtosis_frequency(freq, spec2),
      entropy2 = spectral_entropy(spec2),
      bandwidth2 = spectral_bandwidth(freq, spec2),
      edge_freq2 = spectral_edge_frequency(freq, spec2),
      i1_band2 = band_power_ratio(freq, spec2, 0, 0.5),
      i2_band2 = band_power_ratio(freq, spec2, 0.5, 3),
      i3_band2 = band_power_ratio(freq, spec2, 3, 10),
      i4_band2 = band_power_ratio(freq, spec2, 10, 25),

      # Signal 3
      dom_power_ratio3 = dominant_power_ratio(freq, spec3),
      dom_freq3 = freq[which.max(spec3)],
      dom_hz_freq3 = dominant_frequency_hz(freq, spec3),
      mean_freq3 = mean_frequency(freq, spec3),
      sd_freq3 = sd_frequency(freq, spec3),
      skew_freq3 = skew_frequency(freq, spec3),
      kurt_freq3 = kurtosis_frequency(freq, spec3),
      entropy3 = spectral_entropy(spec3),
      bandwidth3 = spectral_bandwidth(freq, spec3),
      edge_freq3 = spectral_edge_frequency(freq, spec3),
      i1_band3 = band_power_ratio(freq, spec3, 0, 0.5),
      i2_band3 = band_power_ratio(freq, spec3, 0.5, 3),
      i3_band3 = band_power_ratio(freq, spec3, 3, 10),
      i4_band3 = band_power_ratio(freq, spec3, 10, 25),
    )
}
