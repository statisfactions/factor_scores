## Load packages and preliminaries ----
library(tidyverse)
library(readr)
library(simpr)
library(psych)
library(furrr)
library(future)
library(tictoc)
library(testthat)
library(glue)
library(parallel)
library(semTools)
source("sim_fns.R")

## Specify ----
nfactor_data = list.files("nfactor_results", pattern = "RDS$", full.names = TRUE) %>% 
  set_names(.)

nfactor_results = map_dfr(nfactor_data, readRDS, .id = "filename") %>% rename(n_factors = sim_cell) 

#nfactor_results

## Fit efas -----
plan(multicore)
tic()
gen_efa = nfactor_results %>% reproduce(colname = "sim_cell", fn = fit_and_tidy_efa_only, 
                                      globals = TRUE, packages = "GPArotation")
toc()
saveRDS(gen_efa, file = "efas_1.RDS")


# tic()
# regen = nfactor_results %>% reproduce(colname = "sim_cell", fn = genify,
#                                      globals = TRUE)
# toc()
# 
# regen
# 
# x = regen$sim_cell[[3]]

