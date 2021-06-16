## Load packages and preliminaries ----
library(tidyverse)
library(readr)
library(simpr)
library(psych)
# library(furrr)
library(tictoc)
library(testthat)
library(glue)
library(parallel)
source("sim_fns.R")

all_rs = c(.5, .7, .9)
conditions =  names(condition_lists) %>% str_subset("^[1-3]") 
## Specify ----
spec_1 = blueprint(x = ~ mes(get_loadings(condition_lists, factor_condition, num_items), 
                           condition_lists[[factor_condition]]$cor, 
                           marginals = marg_fun(test_candidate, num_items), 
                           n_multiplier = n_multiplier,
                           sd_latent = get_sd_latent(r)))


  

## Do simulation -----
options(mc.cores = parallel::detectCores())

for(k in 1:50) {
  runtime = Sys.time() %>% str_replace_all(":", "_")
	cat("------\n",k, "\n")
for(i in all_rs) {
  cat(i, "\n")
  for(j in conditions) {
    cat(j, "\n")
    tic()
    gen = spec_1 %>%
      meta(r = i,
           n_multiplier = c(10, 15, 20),
           factor_condition = j,
           num_items = c(5, 10, 15)
      ) %>% 
      future_produce2(1, fn = get_nfact, globals = TRUE, use_future = FALSE)
    toc()
    saveRDS(gen, file = glue("nfactor_results/{runtime} {i} {j} {k}.RDS"))
  }
}
  system("git add nfactor_results/*")
  system(glue("git commit -m '{runtime} repetition {k}'"))
  system("git push -u origin main")
}
# 441 for 1

# gen
## Reproduce --------
# regen = gen %>% reproduce(colname = "r1") %>% 
#   reproduce(colname = "r2")
# 
# test_that("reproduce results are always the same", {
#   expect_equal(regen$r1, regen$r2)
# })
# 


