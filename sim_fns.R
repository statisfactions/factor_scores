list_w_names = function(...) {
  dots_names <- sapply(substitute(list(...))[-1], deparse)
  
  structure(list(...), .Names = dots_names)
  
}

mes <- function(fmodel,effect, marginals, numberofcases=1000, sd_latent = 1, sd_error = 1) {   # define a general function in terms of a factor model and an effects matrix
  ## Function altered from http://personality-project.org/r/r.datageneration.html
  numberofvariables <- dim(fmodel)[1]        #problem size determined by input to the function
  numberoflatent <- dim(fmodel)[2]
  numberofcategories = ncol(marginals)
  tmodel <- t(fmodel)      #transpose of model
  # fmodel  %*% tmodel        #show the resulting  measurement structure
  
  communality=diag(fmodel%*%tmodel)       #find how much to weight true scores and errors given the measurement model
  uniqueness=1-communality
  errorweight=sqrt(uniqueness)
  errorweight=diag(errorweight)      #how much to weight the errors
  
  latentscores <- matrix(rnorm(numberofcases*(numberoflatent), sd = sd_latent), numberofcases) #create true scores for the latent variables
  # round(cor(latentscores),2)      #if uncommented, this shows the sample true score correlation matrix of the factors
  latentscores <- latentscores%*%effect      #create true scores to reflect structural relations between the factors
  # round(cor(latentscores),2)           #note that the factors are now correlated
  
  truescores <- latentscores %*% tmodel
  # round(cor(truescores),2)           #show the true score correlation matrix (without error)
  error<- matrix(rnorm(numberofcases*(numberofvariables), sd = sd_error),numberofcases)  #create normal error scores
  error=error%*%errorweight
  observedscore=truescores+error
  
  ## Adapted from source of psych::poly.mat()
 Y <- matrix(0, nrow = numberofcases, ncol = numberofvariables)
  for (i in 1:(numberofcategories)) {
    Y[observedscore > marginals[, i]] <- i
  }
  return(Y)
}       #end of function mes



genify = function(variables, ..., reproduce = FALSE) {
  # Generate the simulation data based on the simpr specification
  
  eval_environment = rlang::as_environment(list(...), parent = parent.frame())
  
  ## If reproducing, random seed should be in the
  ## "seed" variable provided in ...
  if(reproduce) {
    if(exists(eval_environment$seed)) {
      starting_seed <- eval_environment$seed
      assign(".Random.seed", starting_seed, envir = .GlobalEnv)
    } else
      stop("Seed should be provided in the `seed` column for reproduction.")
  } else {
    ## Save random seed if not reproducing
    
    ## Create random seed if it does not already exist.
    if (!exists(".Random.seed", envir = .GlobalEnv))
      runif(1)
    
    ## save the random seed 
    starting_seed <- .GlobalEnv$.Random.seed
  }
  
  df = purrr::map_dfc(variables, function(y) {
    
    eval_fun = purrr::as_mapper(y)
    environment(eval_fun) <- eval_environment
    
    gen = eval_fun()
    
    varnames = attr(y, "varnames")
    
    if(is.null(ncol(gen))) {
      gen_df = tibble::as_tibble(gen, .name_repair = "minimal")
      names(gen_df) = varnames
      assign(varnames, gen, envir = eval_environment)
      
    } else if(length(dim(gen)) > 3) {
      stop("More than 2 dimensional output in blueprint() not supported")
    } else if(ncol(gen) == 0) {
      stop("variable function returns 0 columns")
    } else if(ncol(gen) == 1) {
      names(gen) = attr(y, "varnames")
      assign(varnames, gen[[1]], envir = eval_environment)
    } else if(ncol(gen) > 1) {
      gen_df = tibble::as_tibble(gen, .name_repair = "minimal")
      
      # rename gen_df
      ## if multiple varnames given via formula lhs, use those
      if(length(varnames) > 1) {
        names(gen_df) = varnames
      } else {
        # Otherwise, use auto-numbering
        names(gen_df) = sprintf(paste0("%s%s%0", nchar(trunc(ncol(gen))), ".0f"),
                                varnames,
                                attr(variables, "sep"),
                                1:ncol(gen))
        ## assign names to the eval_environmnent
      }
      purrr::iwalk(gen_df, ~ assign(.y, .x, envir = eval_environment))
      
      
    }
    
    gen_df
  })
  ## save seed as attribute
  attr(df, "seed") = starting_seed
  
  df
  
}

reproduce <-
  function(x, colname = "reproduce", fn = genify, globals = TRUE, packages = NULL) {
    ## Based on simpr::produce() to reproduce, given a seed.
    
    
    ## Generate all replications
    x[[colname]] = x %>% select(-sim_cell) %>% 
      furrr::future_pmap(., fn, variables = attr(x, "variables_spec"), 
                         .options = furrr::future_options(globals = globals,
                                                          packages = packages),
                         .progress = TRUE)
    
    x
  }


future_produce2 <-
  function(x, reps, fn = simpr:::genify, globals = TRUE, packages = NULL) {
    ## Based on simpr::produce() but altered to allow parallel processing
    
    ## Check reps argument is a whole number > 0
    if(length(reps) > 1) {
      reps = reps[1]
      warning("reps has length > 1, using only first element")
    }
    if(reps != round(reps)) {
      reps = round(reps)
      warning("reps not a whole number, rounding to nearest whole number")
    }
    if(reps < 1) {
      stop("reps should be at least 1")
    }
    
    # Create data frame representing all possible values of meta parameter indices
    specs = expand.grid(c(list(rep = 1:reps), x$meta$indices), stringsAsFactors = FALSE)
    
    if(!is.null(x$meta$lookup)) {
      ## If there are list elements, join cells representing those list-columns
      ## into specs
      specs = purrr::reduce2(x$meta$lookup,
                             dplyr::inner_join,
                             .init = specs,
                             .y = names(x$meta$lookup)) # the "by" argument to the join
    }
    
    ## define variable "." to satisfy R CMD Check
    . = "Defining ."
    
    sim_results = specs %>% tibble::as_tibble()
    
    ## Generate all replications
    sim_results$sim_cell = specs %>%
      furrr::future_pmap(., fn, variables = x$variables, 
                         .options = furrr::future_options(globals = globals,
                                                          packages = packages),
                         .progress = TRUE)
    
    sim_results$seed = purrr::map(sim_results$sim_cell, ~ attr(., "seed"))
    
    ## Add some attributes to the tibble to track meta and variables
    attr(sim_results, "meta") = names(x$meta$indices)
    attr(sim_results, "variables") = purrr::map(x$variables, ~ attr(., "varnames")) %>% unlist
    attr(sim_results, "variables_spec") = x$variables
    
    ## Add "simpr_gen" class
    class(sim_results) = c("simpr_gen", class(sim_results))
    
    sim_results
  }

# Power simulation ------
fit_and_tidy = function(variables, ...) {
  ## For the power simulation, generate the data,
  ## fit the model (same for ES and MF), and tidy
  ## the results, extracting just the fixed
  ## effects (all we will use is the p.value for z.edyrscap.mean
  
    
  
  catchToList <- function(expr) {
    ## https://stackoverflow.com/questions/4948361/how-do-i-save-warnings-and-errors-as-output-from-a-function
    val <- NULL
    myWarnings <- NULL
    wHandler <- function(w) {
      myWarnings <<- c(myWarnings, w$message)
      invokeRestart("muffleWarning")
    }
    myError <- NULL
    eHandler <- function(e) {
      myError <<- e$message
      NULL
    }
    val <- tryCatch(withCallingHandlers(expr, warning = wHandler), error = eHandler)
    return_list = list(value = val, warnings = myWarnings, error=myError)
    attr(return_list, "seed") = attr(gen, "seed")
    return(return_list)
  } 
  
  fa_wlsmv = function(x, nf) {
    x = dplyr::mutate(x, across(everything(), ordered, levels = 0:3))
    
    unrotated = semTools::efaUnrotate(x, nf = nf, varList = names(x), estimator = "wlsmv")
    rotated = semTools::oblqRotate(unrotated, method = "geomin")
    
    rotated
  }
  
  gen = genify(variables, ...)
  
  ## Get number of factors
  pop_num_fact = 3
  
  catchToList({
    
    capture.output(nfact <- psych::fa.parallel(gen, cor = "poly", plot = FALSE, sim = FALSE)$nfact, file = "NUL")
    
    ## Don't bother fitting model if wrong number of factors
    if(nfact != 3) {
      attr(nfact, "seed") = attr(gen, "seed")
      return(seed)
    }
    
    fa_fit <- gen %>%
    fa_wlsmv(nfact)
  
  fa_out <- tibble::as_tibble(fa_fit@loading, rownames = "variable") %>% 
    dplyr::bind_rows(tibble::as_tibble(fa_fit@phi, rownames = "variable"))
  
  attr(fa_out, "seed") = attr(gen, "seed")
  
  return(fa_out)
  })
}

## 

