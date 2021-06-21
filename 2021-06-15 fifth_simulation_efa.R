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
library(semTools)
source("sim_fns.R")

## Specify ----

get_nfactor_data = function(done_files) {

nfactor_data = list.files("nfactor_results", pattern = "RDS$", full.names = TRUE) %>%
setdiff(done_files) %>%
  purrr::set_names(.)

nfactor_results = map_dfr(nfactor_data, readRDS, .id = "filename") %>% rename(n_factors = sim_cell) %>% 
  mutate(rep = filename %>% 
           str_remove(".RDS$") %>% str_extract("[0-9]+$") %>% as.integer)            
}

if(file.exists("nfactor_completed.RDS")) {
nfactor_completed = readRDS("nfactor_completed.RDS")
} else {
	nfactor_completed = c()
}


## Fit efas -----
nfiles = 20 # number of EFA files to do simultaneously
plan(multicore)
runtime = Sys.time() %>% str_replace_all(":", "_")
for(j in 1:10000) {

	system("git pull")
nfactor_results = get_nfactor_data(nfactor_completed)

batch_files = sort(unique(nfactor_results$filename))[1:nfiles] 
if(any(str_detect(batch_files, "all_conditions")))
	batch_files = batch_files[1]
	
batch = nfactor_results %>% filter(filename %in% batch_files) 

batch_name = with(batch[1,], paste(r, factor_condition, rep))
cat("Batch ", batch_files, "\n")


  tic()
  gen_efa = batch %>%  reproduce(colname = "sim_cell", fn = fit_and_tidy_efa_only, 
                                                               globals = TRUE, packages = "GPArotation")
  toc() #753 seconds
  nfactor_completed = c(nfactor_completed, gen_efa$filename)
  saveRDS(nfactor_completed, file = "nfactor_completed.RDS")
  saveRDS(gen_efa, file = glue("efa_results/{runtime} {batch_name} efas {length(batch_files)} files.RDS"))
}


# tic()
# regen = nfactor_results[1,] %>% reproduce(colname = "sim_cell", fn = genify,
#                                      globals = TRUE)
# toc()
# 
# 
# regen %>% fit(fit = ~ fa_wlsmv(., nf = 1))

