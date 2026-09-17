#' Parse a file and extract joined time and frequency domain features.
#'
#' @param dir The directory containing the file.
#' @param filename The name of the file to parse.
#' @param sample_labels The labels for each sample.
#' @param n_samples_per_epoch The number of samples per epoch.
#' @param sample_rate The sampling rate of the data.
#' @return A data frame containing the extracted features.
parse_and_features <- function(filename, dir, sample_labels, n_samples_per_epoch, sample_rate) {
  # Parse and load the df
  df <- load_and_parse_file(dir, filename, sample_labels)
  df_freq <- convert_signal_to_spectrum_df(df, n_samples_per_epoch, sample_rate)

  # Extract relevant time and freq features based on the sensor type
  params <- get_file_params(filename)
  if (params$sensor_type == "acc") {
    # Process accelerometer data
    features_time <- get_time_domain_features(df, n_samples_per_epoch, sample_rate)
    features_freq <- get_frequency_domain_features(df_freq)
  } else if (params$sensor_type == "gyro") {
    # Process gyroscope data
    features_time <- get_time_domain_features(df, n_samples_per_epoch, sample_rate)
    features_freq <- get_frequency_domain_features(df_freq)
  }

  # Left join features by epoch
  features <- left_join(features_time, features_freq, by = "epoch")

  # Add cols for user_id and exp_id from params
  features <- features |> mutate(user_id = params$user_id, exp_id = params$exp_id)

  return(features)
}
