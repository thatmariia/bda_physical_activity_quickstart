library(tidyverse)
library(caret, warn.conflicts = FALSE)
options(width = 100)

on_kaggle <- dir.exists("/kaggle/input")
final_submission <- FALSE
data_dir <- if (on_kaggle) {
  list.files("/kaggle/input/competitions", full.names = TRUE)[1]
} else {
  here::here("data")
}
output_dir <- if (on_kaggle) "/kaggle/working" else here::here("output")
dir.create(output_dir, showWarnings = FALSE)

train_dir <- file.path(data_dir, "RawData", "Train")
test_dir <- file.path(data_dir, "RawData", "Test")
