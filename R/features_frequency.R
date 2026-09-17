# Feature functions (frequency domain)

#' Compute the mean frequency of a signal
mean_frequency <- function(freq, spec) {
  delta_f <- freq[2] - freq[1] # gap size between frequencies
  normalizing_constant <- sum(spec * delta_f) # ≈ ∫S(f)df
  mean_freq <- sum(freq * spec * delta_f) / normalizing_constant
  return(mean_freq)
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
  if (length(spec) == 0 || sum(spec) == 0) return(0)
  spec <- spec[spec > 0]
  spec <- spec / sum(spec)
  return(-sum(spec * log2(spec)) / log2(length(spec)))
}

#' Compute the spectral bandwidth of a signal
spectral_bandwidth <- function(freq, spec) {
  mean_freq <- mean_frequency(freq, spec)
  tot <- sum(spec)
  if (tot == 0) return(0)
  ds <- (freq - mean_freq)^2
  return(sqrt(sum(ds * spec) / tot))
}

#' Compute the spectral edge frequency of a signal
spectral_edge_frequency <- function(freq, spec, k = 0.9) {
  tot <- sum(spec)
  if (tot == 0) return(0)
  cumsum_spec <- cumsum(spec)
  threshold <- k * tot
  edge_freq <- which(cumsum_spec >= threshold)[1]
  return(freq[edge_freq])
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
      dom_freq1 = freq[which.max(spec1)],
      mean_freq1 = mean_frequency(freq, spec1),
      entropy1 = spectral_entropy(spec1),
      bandwidth1 = spectral_bandwidth(freq, spec1),
      edge_freq1 = spectral_edge_frequency(freq, spec1),

      # Signal 2
      dom_freq2 = freq[which.max(spec2)],
      mean_freq2 = mean_frequency(freq, spec2),
      entropy2 = spectral_entropy(spec2),
      bandwidth2 = spectral_bandwidth(freq, spec2),
      edge_freq2 = spectral_edge_frequency(freq, spec2),

      # Signal 3
      dom_freq3 = freq[which.max(spec3)],
      mean_freq3 = mean_frequency(freq, spec3),
      entropy3 = spectral_entropy(spec3),
      bandwidth3 = spectral_bandwidth(freq, spec3),
      edge_freq3 = spectral_edge_frequency(freq, spec3),
  )
}
