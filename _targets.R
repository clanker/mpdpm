library(targets)
library(tarchetypes)

tar_option_set(
	packages = c("rsample", "mclust", "mvtnorm", "MASS", "future", "ggplot2"),
  workspace_on_error = TRUE
)

options(
  tidyverse.quiet - TRUE,
  clustermq.scheduler = "multiprocess"
)
#options(clustermq.scheduler = "multicore")
# use tar_make_clustermq(workers = #) as the tar_make command

#future::availableCores()
clustermq::register_dopar_cmq(n_jobs = parallel::detectCores(logical = FALSE) - 1)

function_path <- here::here("script/functions")
function_files <- list.files(pipelines_path, pattern = 'R$', full.names = TRUE)
purrr::walk(function_files, source)

list(
  tar_target(iter_sample, seq_len(10)),
  tar_target(data_raw, generate_data(seed = iter_sample),
             pattern = map(iter_sample),
             iteration = "list"),
  tar_target(data_split, split_data(data_raw, seed = NULL),
             pattern = map(data_raw),
             iteration = "list"),
  tar_target(mclust_fit, fit_mclust(data_split),
             pattern = map(data_split),
             iteration = "list"),
  tar_target(method_list, "dpmg"),
  tar_target(numerator_list, seq(2, 20, 2)),
  tar_target(model_fits, fit_model(data_split, mclust_fit, param = numerator_list,
                                   seed = iter_sample, method = method_list),
             pattern = cross(map(data_split, mclust_fit, iter_sample), numerator_list)),
  tar_target(save_metrics_csv, readr::write_csv(model_fits, "results/v3-out.csv")),
  tar_render(report, "script/analysis/Report.Rmd", params = list(data_split, model_fits)),
  tar_render(run_testcases, "script/analysis/Work 7-17 no 3_211122.Rmd"),
  tar_render(run_sigmatest, "script/analysis/sigmatest_v1a.Rmd")
)
