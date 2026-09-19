# ==========================================================
# == FUNCTIONS FOR COMBINING FEATURES
# ==========================================================

#' Parse a file and extract joined time and frequency domain features.
#'
#' @param dir The directory containing the file.
#' @param filename The name of the file to parse.
#' @param sample_labels The labels for each sample.
#' @param n_samples_per_epoch The number of samples per epoch.
#' @param sample_rate The sampling rate of the data.
#' @param overlap The share of each epoch that overlaps with the next one
#'   (0 = no overlap, 0.5 = a new epoch starts every half epoch).
#' @return A data frame containing the extracted features.
parse_and_features <- function(filename_acc, dir, sample_labels, n_samples_per_epoch, sample_rate, overlap = 0) {
  # check that filename starts with "acc"
  if (!grepl("^acc", filename_acc)) {
    stop("File must start with 'acc'")
  }
  filename_gyro <- sub("^acc", "gyro", filename_acc)

  # Parse and load the dfs
  df_acc <- load_and_parse_file(dir, filename_acc, sample_labels)
  df_gyro <- load_and_parse_file(dir, filename_gyro, sample_labels)

  # With overlap, epochs start every `step` samples instead of every epoch
  step <- max(1, round(n_samples_per_epoch * (1 - overlap)))
  offsets <- seq(0, n_samples_per_epoch - 1, by = step)

  # Extract features for the epochs starting at each offset
  features <- map(offsets, \(offset) {
    # Shift the signals so that the epochs start at `offset`
    shift <- \(df) df |>
      filter(sampleid >= offset) |>
      mutate(sampleid = sampleid - offset)

    get_epoch_features(shift(df_acc), shift(df_gyro), n_samples_per_epoch, sample_rate) |>
      # Shift the start sample of each epoch back to its original position
      mutate(sampleid = sampleid + offset)
  }) |>
    list_rbind()

  # Add cols for user_id and exp_id from params
  params <- get_file_params(filename_acc)
  features <- features |> mutate(user_id = params$user_id, exp_id = params$exp_id)

  return(features)
}

#' Extract joined time, frequency and joint features for each epoch
#'
#' @param df_acc,df_gyro Signal data frames from `load_and_parse_file()`.
#' @param n_samples_per_epoch The number of samples per epoch.
#' @param sample_rate The sampling rate of the data.
#' @return A data frame with one row per epoch.
get_epoch_features <- function(df_acc, df_gyro, n_samples_per_epoch, sample_rate) {
  df_freq_acc <- convert_signal_to_spectrum_df(df_acc, n_samples_per_epoch, sample_rate)
  df_freq_gyro <- convert_signal_to_spectrum_df(df_gyro, n_samples_per_epoch, sample_rate)
  joint_df <- join_sensors(df_acc, df_gyro)

  # Extract features from each sensor
  features_time_acc <- get_time_domain_features(df_acc, n_samples_per_epoch, sample_rate)
  features_freq_acc <- get_frequency_domain_features(df_freq_acc)
  features_time_gyro <- get_time_domain_features(df_gyro, n_samples_per_epoch, sample_rate)
  features_freq_gyro <- get_frequency_domain_features(df_freq_gyro)
  features_joint <- get_joint_features(joint_df, n_samples_per_epoch)

  # Keep the epoch information from the accelerometer only
  epoch_info <- c("sampleid", "n_samples", "aggr_activity", "activity_confidence")
  features_time_gyro <- features_time_gyro |> select(-all_of(epoch_info))

  # Left join all the features
  time_features <- left_join(features_time_acc, features_time_gyro, by = "epoch", suffix = c("_acc", "_gyro"))
  freq_features <- left_join(features_freq_acc, features_freq_gyro, by = "epoch", suffix = c("_acc", "_gyro"))
  sensor_features <- left_join(time_features, freq_features, by = "epoch", suffix = c("_time", "_freq"))
  features <- left_join(sensor_features, features_joint, by = "epoch", suffix = c("", "_joint"))

  return(features)
}
