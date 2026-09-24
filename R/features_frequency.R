# ==========================================================
# == FUNCTION FOR EXTRACTING FREQUENCY DOMAIN FEATURES
# ==========================================================


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
