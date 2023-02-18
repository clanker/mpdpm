library(targets)
library(targetypes)

tar_option_set(
  packages = c("future", "ggplot2"),
  workspace_on_error = TRUE
)

options(
  tidyverse.quiet - TRUE,
  clustermq.scheduler = "multiprocess"
)

future::availableCores()
clustermq::register_dopar_cmq(n_jobs = parallel::detectCores(logical = FALSE) - 1)

pipelines_path <- here::here("script/pipelines")
pipelines_files <- list.files(pipelines_path, pattern = 'R$', full.names = TRUE)
purrr::walk(pipelines_files, source)

pipe_data <- list(
  tar_target()
)

pipe_prep <- list(

)

pipe_report <- list(

)

c(pipe_data, pipe_prep, pipe_report)
