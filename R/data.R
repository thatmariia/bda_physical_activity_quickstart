# Helper functions

#' Create a data frame with signal data and corresponding activity labels
#' @param exp_id The experiment ID (integer)
#' @param user_id The user ID (integer)
#' @param signal_data A data frame containing the signal data
#' @param sample_labels A data frame containing the sample labels with columns: userid, trial, sampleid, act_code, activity, segment
#' @return A data frame combining the signal data with the corresponding activity labels
get_signal_df <- function(exp_id, user_id, signal_data, sample_labels) {
  user_df <-
    tibble(userid = user_id, trial = exp_id, sampleid = seq.int(0, nrow(signal_data) - 1)) |>
    bind_cols(signal_data) |>
    left_join(sample_labels, by = c("userid", "trial", "sampleid"))
  return(user_df)
}

#' Get the file path for a specific sensor type, experiment ID, and user ID
#' @param dir The directory where the data files are located
#' @param sensor_type The type of sensor (e.g., "acc", "gyro")
#' @param exp_id The experiment ID (integer)
#' @param user_id The user ID (integer)
#' @return The constructed file path as a string
get_file_path <- function(dir, sensor_type, exp_id, user_id) {
  # Construct the filename based on the provided parameters
  exp_id_str <- sprintf("%02d", exp_id)
  user_id_str <- sprintf("%02d", user_id)
  file_path <- file.path(dir, paste0(sensor_type, "_exp", exp_id_str, "_user", user_id_str, ".txt"))
  return(file_path)
}

#' Load signal data from a file
#' @param file_path The path to the signal file
#' @return A data frame containing the signal data
load_signal <- function(file_path) {
  signal_data <- read_delim(file_path, delim = " ", col_names = FALSE, col_types = "ddd", progress = FALSE)
  return(signal_data)
}

#' Extract parameters from a filename
#' @param filename The name of the file
#' @return A list containing the user ID, experiment ID, and sensor type
get_file_params <- function(filename) {
  # extract parts of the filename <sensor>_exp<XX>_user<YY>.txt
  user_id <- str_extract(filename, "(?<=user)\\d+") |> as.integer()
  exp_id <- str_extract(filename, "(?<=exp)\\d+") |> as.integer()
  sensor_type <- str_extract(filename, "^(acc|gyro)")
  return(list(user_id = user_id, exp_id = exp_id, sensor_type = sensor_type))
}

#' Load and parse a signal file into a data frame
#' @param dir The directory where the file is located
#' @param filename The name of the file
#' @param sample_labels A data frame containing the sample labels
#' @return A data frame containing the signal data with corresponding activity labels
load_and_parse_file <- function(dir, filename, sample_labels) {
  file_path <- file.path(dir, filename)
  params <- get_file_params(filename)
  signal_data <- load_signal(file_path)
  signal_df <- get_signal_df(params$exp_id, params$user_id, signal_data, sample_labels)
  return(signal_df)
}

#' Get the range of user IDs and experiment IDs in a directory
#' @param dir The directory to search
#' @return A list with two elements: user_ids and exp_ids
range_users_exps <- function(dir) {
  params <- params_from_files(dir)
  return(list(
    user_ids = unique(params$user_id),
    exp_ids = unique(params$exp_id)
  ))
}
