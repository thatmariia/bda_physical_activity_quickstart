# Converts notebook.Rmd to output/notebook.ipynb for Kaggle.
# The code is not run here; Kaggle runs it.

# Every file in R/ must be shown in one of the notebooks in analysis/
analysis_text <- unlist(lapply(list.files(here::here("analysis"), full.names = TRUE), readLines))
for (r_file in list.files(here::here("R"))) {
  if (!any(grepl(paste0('"R", "', r_file, '"'), analysis_text, fixed = TRUE))) {
    stop("R/", r_file, " is not included in any notebook in analysis/", call. = FALSE)
  }
}

ipynb <- rmarkdown::output_format(
  knitr = rmarkdown::knitr_options(
    opts_chunk = list(eval = FALSE),
    opts_hooks = list(child = function(options) {
      options$eval <- TRUE # still include the notebooks in analysis/
      options
    }),
    knit_hooks = list(source = function(x, options) {
      code <- paste(x, collapse = "\n")
      cell <- paste0("```r\n", code, "\n```\n:::\n")
      if (is.null(options$fold)) {
        return(paste0("\n::: {.cell .code}\n", cell))
      }
      # Folded chunks: a title, and a cell that is collapsed below
      paste0(
        "\n*", options$fold, "*\n\n",
        "::: {.cell .code tags=\"[\\\"hide-input\\\"]\"}\n", cell
      )
    })
  ),
  pandoc = rmarkdown::pandoc_options(to = "ipynb")
)

notebook_file <- rmarkdown::render(
  here::here("notebook.Rmd"),
  output_format = ipynb,
  output_dir = here::here("output")
)

# Collapse the marked cells on Kaggle and in Jupyter
notebook <- jsonlite::read_json(notebook_file)
notebook$cells <- lapply(notebook$cells, function(cell) {
  if ("hide-input" %in% unlist(cell$metadata$tags)) {
    cell$metadata[["_kg_hide-input"]] <- TRUE
    cell$metadata$jupyter <- list(source_hidden = TRUE)
  }
  cell
})
jsonlite::write_json(notebook, notebook_file, auto_unbox = TRUE, pretty = TRUE, null = "null")
