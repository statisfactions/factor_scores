## Load packages and preliminaries ----
library(tidyverse)
library(readr)
library(simpr)
library(psych)
library(furrr)
library(tictoc)
library(testthat)
library(glue)
source("sim_fns.R")

## Specify ----
spec = blueprint(x = ~ mes(get_loadings(condition_lists, factor_condition, num_items), 
                           condition_lists[[factor_condition]]$cor, 
                           marginals = marg_fun(test_candidate, num_items), 
                           n_multiplier = n_multiplier,
                           sd_latent = get_sd_latent(r))) %>% 
  meta(r = c(.1, .3, .5, .7, .9),
       n_multiplier = c(10, 15, 20),
       factor_condition = names(condition_lists),
       num_items = c(5, 10, 15, 20)
       )


## Do simulation -----
set.seed(21333)

plan(multicore)
filenames = glue("sim_results/{Sys.time()} {1:50}.Rdata")
for(i in filenames) {
  print(i)
tic()
gen = spec %>% 
  future_produce2(1, fn = fit_and_tidy, globals = TRUE, packages = "GPArotation")
toc()
saveRDS(gen, file = i)
}
# 441 for 1

gen$sim_cell
## Reproduce --------
regen = gen %>% reproduce(colname = "r1") %>% 
  reproduce(colname = "r2")

test_that("reproduce results are always the same", {
  expect_equal(regen$r1, regen$r2)
})



