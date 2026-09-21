#' Get sample labels
#' @param lbl_dir The directory containing the activity labels
#' @param lbl_filename The name of the file containing activity labels
#' @param lbl_train_dir The directory containing the training labels
#' @param lbl_train_filename The name of the file containing training labels
get_sample_labels <- function(lbl_dir, lbl_filename, lbl_train_dir, lbl_train_filename) {
  labels_list <- load_labels(lbl_dir, lbl_filename)
  labels_train <- load_labels_train(lbl_train_dir, lbl_train_filename)

  # Join the activity labels to the labels data frame
  labels <- labels_train |> left_join(labels_list, by = c(act_code = "code"))

  # Expand the labels data frame to have one row per sample
  sample_labels <- labels |>
    mutate(segment = row_number()) |> # identify each separate activity interval
    reframe(
      sampleid = seq.int(start, end), # expand start:end into all sample IDs
      .by = c(trial, userid, act_code, activity, segment)
    )
  return(sample_labels)
}

#' Load activity labels
load_labels <- function(dir, filename) {
  read.table(
    file.path(data_dir, filename),
    header = FALSE, col.names = c("code", "activity")
  )
}

#' Load training labels
load_labels_train <- function(dir, filename) {
  read_delim(file.path(dir, filename),
    delim = " ",
    col_names = c("trial", "userid", "act_code", "start", "end"),
    col_types = "iiiii"
  )
}
