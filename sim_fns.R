load("condition_matrices.Rdata")

loading_con_df = tibble(loading = names(loading_conditions)) %>% 
  mutate(nfactors = str_extract(loading, pattern = "[0-9]+$"))


cor_con_df = tibble(cor = names(cor_conditions)) %>% 
  mutate(nfactors = str_extract(cor, pattern = "[0-9]+$"))

all_cond_df = inner_join(loading_con_df, cor_con_df, by = "nfactors") %>% 
  mutate(condition_name = paste(nfactors, loading, cor) %>% str_remove_all("_[0-9]+"))

condition_df = all_cond_df %>% 
  group_by(condition_name) %>% 
  summarize(cond = list(list(loading = loading_conditions[[loading]], cor = cor_conditions[[cor]])),
            .groups = "drop")

condition_lists = condition_df$cond
names(condition_lists) = condition_df$condition_name


get_sd_latent = function(r) sqrt(1 / (1-r))
get_sd_total = function(r) sqrt(get_sd_latent(r)^2 +1)

rs = c(.1, .3, .5, .7, .9)

marg_fun = function(x, num_items) {
  matrix(rep(seq(-abs(x), abs(x), length.out = 3), num_items), ncol = 3, byrow = TRUE) 
}

marg_fun_pop = function(r, x) {
  x = abs(x)
  sd_t = get_sd_total(r)
  c0 = pnorm(-x, mean = 0, sd = sd_t, lower.tail = TRUE)
  
  data.frame(c0 = c0, 
             c1 = 0.5 - c0,
             c2 = 0.5 - c0,
             c3 = c0)
}

get_loadings = function(master, cond, num_items, master_items_per_factor = 5) {
  ## Get the population loadings given a master
  ## condition list of matrices
  ## (loading_conditions), the condition
  ## (character vector), the number of items
  ## (numeric), and the number of items per factor
  ## in the master condition list (numeric, default 5)
  do.call(rbind, rep(list(master[[cond]]$loading), round(num_items/master_items_per_factor)))
}
## We want the max-reliability condition to have
## at least 10% (or another %) in the center
## This means the quantile we want is .495/2 = 0.2479
test_candidate = qnorm(0.2475, mean = 0, get_sd_total(.9))

g = marg_fun_pop(.9, test_candidate)


list_w_names = function(...) {
  dots_names <- sapply(substitute(list(...))[-1], deparse)
  
  structure(list(...), .Names = dots_names)
  
}

seed_factory = function(f) {
  ## Function factory: Return a function that
  ## captures .Random.seed and returns a list with
  ## seed and original object
  function(...) {
  ## Create random seed if it does not already exist.
  if (!exists(".Random.seed", envir = .GlobalEnv))
    runif(1)
  
  ## save the random seed 
  starting_seed <- .GlobalEnv$.Random.seed
  res = f(...)
  list(sim_cell = res, seed = starting_seed)
  }
}

mes <- function(fmodel,effect, marginals, n_multiplier, sd_latent = 1, sd_error = 1) {   # define a general function in terms of a factor model and an effects matrix
  ## Function altered from http://personality-project.org/r/r.datageneration.html
  numberofvariables <- dim(fmodel)[1]        #problem size determined by input to the function
  numberoflatent <- dim(fmodel)[2]
  numberofcategories = ncol(marginals)
  numberofcases = n_multiplier * numberofvariables
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
    if(exists("seed", envir = eval_environment)) {
      starting_seed <- eval_environment$seed
      assign(".Random.seed", starting_seed, envir = .GlobalEnv)
    } else
      stop("Seed should be provided in the `seed` column for reproduction.")
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
  
  df
  
}

reproduce <-
  function(x, colname = "reproduce", fn = genify, globals = TRUE, packages = NULL) {
    ## Based on simpr::produce() to reproduce, given a seed.
    
    
    ## Generate all replications
    x[[colname]] = x %>% select(-matches("^sim_cell$")) %>% 
      furrr::future_pmap(., fn, variables = attr(x, "variables_spec"), reproduce = TRUE,
                         .options = furrr::future_options(globals = globals,
                                                          packages = packages),
                         .progress = TRUE)
    
    x
  }


future_produce2 <-
  function(x, reps, fn = simpr:::genify, globals = TRUE, packages = NULL, use_future = TRUE) {
    ## Based on simpr::produce() but altered to allow parallel processing
    if(use_future) {
      pmap_fn = furrr::future_pmap
    } else
      pmap_fn = purrr::pmap
    
    
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
      pmap_fn(., seed_factory(fn), variables = x$variables, 
                         .options = furrr::future_options(globals = globals,
                                                          packages = packages),
                         .progress = TRUE)
    
    sim_results = unnest_wider(sim_results, col = sim_cell)
    
    ## Add some attributes to the tibble to track meta and variables
    attr(sim_results, "meta") = names(x$meta$indices)
    attr(sim_results, "variables") = purrr::map(x$variables, ~ attr(., "varnames")) %>% unlist
    attr(sim_results, "variables_spec") = x$variables
    
    ## Add "simpr_gen" class
    class(sim_results) = c("simpr_gen", class(sim_results))
    
    return(sim_results)
    
  }


get_nfact = function(variables, ...) {
  gen = genify(variables, ...)
  
  capture.output(nfact <- psych::fa.parallel(gen, cor = "poly", plot = FALSE, sim = FALSE, correct = 0.1)$nfact, file = "/dev/null")
  
  return(nfact)
}

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
    return(return_list)
  } 
  
  fa_wlsmv = function(x, nf) {
    x = dplyr::mutate(x, across(everything(), ordered, levels = 0:3))
    
    unrotated = semTools::efaUnrotate(x, nf = nf, varList = names(x), estimator = "wlsmv")
    if(nf == 1)
      return(unrotated)
    else {
      rotated = semTools::oblqRotate(unrotated, method = "geomin")
      return(rotated)
    }
  }
  
  gen = genify(variables, ...)
  
  ## Get number of factors
  
  catchToList({
    
    capture.output(nfact <- psych::fa.parallel(gen, cor = "poly", plot = FALSE, sim = FALSE, correct = 0.1)$nfact, file = "NUL")
    
    fa_fit <- gen %>%
      fa_wlsmv(nfact)
    
    if(nfact == 1) {
      return(lavaan::parameterEstimates(fa_fit, standardized = TRUE))
    } else {
      fa_out <- tibble::as_tibble(fa_fit@loading, rownames = "variable") %>% 
        dplyr::bind_rows(tibble::as_tibble(fa_fit@phi, rownames = "variable"))
      
      return(fa_out)
    }
  })
}

fit_and_tidy_efa_only = function(variables, ..., reproduce = TRUE) {
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
    return(return_list)
  } 
  
  fa_wlsmv = function(x, nf) {
    x = dplyr::mutate(x, across(everything(), ordered, levels = 0:3))
    
    unrotated = semTools::efaUnrotate(x, nf = nf, varList = names(x), estimator = "wlsmv")
    if(nf == 1)
      return(unrotated)
    else {
      rotated = semTools::oblqRotate(unrotated, method = "geomin")
      return(rotated)
    }
  }
  
  gen = genify(variables, ..., reproduce = reproduce)
  
  ## Get number of factors
  nfact = list(...)$n_factors
 
  print(summary(gen))
  catchToList({
    fa_fit <- gen %>%
      fa_wlsmv(nfact)
    
    if(nfact == 1) {
      return(lavaan::parameterEstimates(fa_fit, standardized = TRUE))
    } else {
      fa_out <- tibble::as_tibble(fa_fit@loading, rownames = "variable") %>% 
        dplyr::bind_rows(tibble::as_tibble(fa_fit@phi, rownames = "variable"))
      
      return(fa_out)
    }
  })
}

## Altered psych functions -----
#1/2/14  switched the n.iter loop to a mclapply loop to allow for multicore parallel processing

fa.parallel_fixed <-
  function(x,n.obs=NULL,fm="minres",fa="both",nfactors=1,main="Parallel Analysis Scree Plots",n.iter=20,error.bars=FALSE,se.bars=FALSE,SMC=FALSE,ylabel=NULL,show.legend=TRUE,sim=TRUE,quant=.95,cor="cor",use="pairwise",plot=TRUE,correct=.5)  { 
    a <- "This is the fixed version"
    
    cl <- match.call()
    # if(!require(parallel)) {message("The parallel package needs to be installed to run mclapply")}
    
    ci <- 1.96
    arrow.len <- .05
    nsub <- dim(x)[1]
    nvariables <- dim(x)[2]
    resample <- TRUE  #this is used as a flag for correlation matrices
    if((isCorrelation(x)) && !sim)  {warning("You specified a correlation matrix, but asked to just resample (sim was set to FALSE).  This is impossible, so sim is set to TRUE")
      sim <- TRUE
      resample <- FALSE}
    if (!is.null(n.obs)) { nsub <- n.obs 
    rx <- x
    resample <- FALSE
    if(dim(x)[1] != dim(x)[2]) {warning("You specified the number of subjects, implying a correlation matrix, but do not have a correlation matrix, correlations found ")
      #	rx <- cor(x,use="pairwise") 
      #add the option to choose the type of correlation, this allows us to do fa.parallel.poly inside fa.parallel
      switch(cor, 
             cor = {rx <- cor(x,use=use)},
             cov = {rx <- cov(x,use=use) 
             covar <- TRUE},
             tet = {rx <- tetrachoric(x,correct=correct)$rho},
             poly = {rx <- polychoric(x,correct=correct)$rho},
             mixed = {rx <- mixedCor(x,use=use,correct=correct)$rho},
             Yuleb = {rx <- YuleCor(x,,bonett=TRUE)$rho},
             YuleQ = {rx <- YuleCor(x,1)$rho},
             YuleY = {rx <- YuleCor(x,.5)$rho } 
      )
      
      if(!sim) {warning("You specified a correlation matrix, but asked to just resample (sim was set to FALSE).  This is impossible, so sim is set to TRUE")
        sim <- TRUE
        resample <- FALSE}
    }   	 } else {
      if (isCorrelation(x)) {warning("It seems as if you are using a correlation matrix, but have not specified the number of cases. The number of subjects is arbitrarily set to be 100  ") 
        rx <- x
        nsub = 100
        n.obs=100
        resample <- FALSE
      }  else {
        
        switch(cor, 
               cor = {rx <- cor(x,use=use)},
               cov = {rx <- cov(x,use=use) 
               covar <- TRUE},
               tet = {rx <- tetrachoric(x,correct=correct)$rho},
               poly = {rx <- polychoric(x,correct=correct)$rho},
               mixed = {rx <- mixedCor(x,use=use,correct=correct)$rho},
               Yuleb = {rx <- YuleCor(x,,bonett=TRUE)$rho},
               YuleQ = {rx <- YuleCor(x,1)$rho},
               YuleY = {rx <- YuleCor(x,.5)$rho } 
        )
      } }
    
    
    
    valuesx  <- eigen(rx)$values #these are the PC values
    if(SMC) {diag(rx) <- smc(rx)
    fa.valuesx <- eigen(rx)$values} else {
      fa.valuesx  <- fa(rx,nfactors=nfactors,rotate="none", fm=fm,warnings=FALSE)$values}  #these are the FA values
    
    temp <- list(samp =vector("list",n.iter),samp.fa = vector("list",n.iter),sim=vector("list",n.iter),sim.fa=vector("list",n.iter))
    
    #parallel processing starts here  - the more cores the better!
    #however, mixedCor seems to break this
    # templist <- lapply(1:n.iter,function(XX) {  
    
    templist <- mclapply(1:n.iter, function(XX) { #at least for now, the errors from mixedCor prevent mclapply 
      
      if(is.null(n.obs)) {
        #Sample the data, column wise (to keep the basic distributional properties, but making the correlations 0 (on average)
        bad <- TRUE
        while(bad) {sampledata <- matrix(apply(x,2,function(y)   sample(y,nsub,replace=TRUE)),ncol=nvariables) #do it column wise
        colnames(sampledata) <- colnames(x)   #this allows mixedCor to work               
        switch(cor,          #we can do a number of different types of correlations
               cor = {C <- cor(sampledata,use=use)},
               cov = {C <- cov(sampledata,use=use) 
               covar <- TRUE},
               tet = {C <- tetrachoric(sampledata,correct=correct)$rho},
               poly = {C <- polychoric(sampledata,correct=correct)$rho},
               mixed = {C <- mixedCor(sampledata,use=use,correct=correct)$rho},
               Yuleb = {C <- YuleCor(sampledata,,bonett=TRUE)$rho},
               YuleQ = {C <- YuleCor(sampledata,1)$rho},
               YuleY = {C <- YuleCor(sampledata,.5)$rho } 
        )  
        
        bad <- any(is.na(C))   #some (not frequently) correlations will be improper, particularly if sampling from sparse matrices                 
        }  #Try resampling until we get a correlation matrix that works                    
        values.samp <- eigen(C)$values
        temp[["samp"]] <- values.samp
        if (fa!= "pc") {
          if(SMC) {sampler <- C 
          diag(sampler) <- smc(sampler)
          temp[["samp.fa"]]<- eigen(sampler)$values} else {
            temp[["samp.fa"]]  <- fa(C,fm=fm,nfactors=nfactors, SMC=FALSE,warnings=FALSE)$values
          }
        }
      } 
      
      if(sim) { simdata=matrix(rnorm(nsub*nvariables),nrow=nsub,ncol=nvariables)    #create simulated data based upon normal theory
      sim.cor <- cor(simdata)   #we must use correlations based upon Pearson here, because we are simulating the data
      
      temp[["sim"]] <- eigen(sim.cor)$values
      
      if (fa!="pc") {
        if(SMC) { diag(sim.cor) <- smc(sim.cor)
        temp[["sim.fa"]]<- eigen(sim.cor)$values} else {fa.values.sim <- fa(sim.cor,fm=fm,nfactors=nfactors,rotate="none",SMC=FALSE,warnings=FALSE)$values
        temp[["sim.fa"]]    <- fa.values.sim
        }}}
      replicates <- list(samp=temp[["samp"]],samp.fa=temp[["samp.fa"]],sim=temp[["sim"]],sim.fa=temp[["sim.fa"]])
    }, mc.set.seed = TRUE)
    #parallelism stops here
    #now combine the results   	
    
    if(is.null(ylabel)) {
      ylabel <- switch(fa,           #switch implementation suggested by Meik Michalke   3/20/17
                       pc = "eigen values of principal components",
                       fa = "eigen values of principal factors",
                       both =  "eigenvalues of principal components and factor analysis")
    }
    
    
    values<- t(matrix(unlist(templist),ncol=n.iter))
    
    
    values.sim.mean=colMeans(values,na.rm=TRUE)
    # if(!missing(quant)) {values.ci = apply(values,2,function(x) quantile(x,quant))} else {values.ci <- values.sim.mean} #fixed Sept 22, 2018
    values.ci = apply(values,2,function(x) quantile(x,quant)) #always apply quant
    if(se.bars) {values.sim.se <- apply(values,2,sd,na.rm=TRUE)/sqrt(n.iter)} else {values.sim.se <- apply(values,2,sd,na.rm=TRUE)}
    
    ymin <- min(valuesx,values.sim.mean)
    ymax <- max(valuesx,values.sim.mean)
    sim.pcr <- sim.far <- NA
    
    switch(fa, 
           pc = {    if (plot) { plot(valuesx,type="b", main = main,ylab=ylabel ,ylim=c(ymin,ymax),xlab="Component Number",pch=4,col="blue")}
             if(resample) { sim.pcr <- values.sim.mean[1:nvariables] 
             sim.pcr.ci <- values.ci[1:nvariables]
             sim.se.pcr <- values.sim.se[1:nvariables]
             if (plot) { points(sim.pcr,type ="l",lty="dashed",pch=4,col="red")}} else {sim.pcr <- NA
             sim.se.pc <- NA}
             if(sim) {    
               if(resample) {sim.pc <- values.sim.mean[(nvariables+1):(2*nvariables)] 
               sim.pc.ci <- values.ci[(nvariables+1):(2*nvariables)] 
               sim.se.pc <- values.sim.se[(nvariables+1):(2*nvariables)] 
               } else {sim.pc <- values.sim.mean[1:nvariables]
               sim.pc.ci <- values.ci[1:nvariables]
               sim.se.pc <- values.sim.se[1:nvariables]}
               if (plot) { points(sim.pc,type ="l",lty="dotted",pch=4,col="red")}
               pc.test <- which(!(valuesx > sim.pc.ci))[1]-1} else { 
                 sim.pc <- NA
                 sim.pc.ci <- NA 
                 sim.se.pc <- NA
                 pc.test <- which(!(valuesx > sim.pcr.ci))[1]-1  }
             fa.test <- NA
             sim.far <- NA
             sim.fa <- NA 
           },
           fa = { #ylabel <-  "eigen values of principal factors"     should not be changed if set (reported by Meik Michalke)
             if (plot) {plot(fa.valuesx,type="b", main = main,ylab=ylabel ,ylim=c(ymin,ymax),xlab="Factor Number",pch=2,col="blue")}
             sim.se.pc <- NA
             if(resample) {sim.far <- values.sim.mean[(nvariables+1):(2*nvariables)]
             sim.far.ci <- values.ci[(nvariables+1):(2*nvariables)]
             sim.se.far <- values.sim.se[(nvariables+1):(2*nvariables)]
             if (plot) {  points(sim.far,type ="l",lty="dashed",pch=2,col="red")}}
             if(sim) { if(resample) {sim.fa <- values.sim.mean[(3*nvariables+1):(4*nvariables)] 
             sim.fa.ci <- values.ci[(3*nvariables+1):(4*nvariables)]
             sim.se.fa <- values.sim.se[(3*nvariables+1):(4*nvariables)] } else {
               sim.fa <- values.sim.mean[(nvariables+1):(2*nvariables)]
               sim.fa.ci <- values.sim.mean[(nvariables+1):(2*nvariables)]
               sim.se.fa <- values.sim.se[(nvariables+1):(2*nvariables)]
               sim.far <- NA  #added May 1, 2016
               sim.far.ci <- NA
               sim.se.far <- NA
             }
               if (plot) {points(sim.fa,type ="l",lty="dotted",pch=2,col="red")}
               fa.test <- which(!(fa.valuesx > sim.fa.ci))[1]-1 
             } else {sim.fa <- NA 
             fa.test <- which(!(fa.valuesx > sim.far.ci))[1]-1 }
             sim.pc <- NA
             sim.pcr <- NA
             sim.se.pc <- NA
             pc.test <- NA },
           both = {     if (plot) {plot(valuesx,type="b", main = main,ylab=ylabel ,ylim=c(ymin,ymax),xlab="Factor/Component Number",pch=4,col="blue")
             points(fa.valuesx,type="b",pch=2,col="blue")}
             if(sim) {  
               if(resample) {
                 sim.pcr <- values.sim.mean[1:nvariables]
                 sim.pcr.ci <- values.ci[1:nvariables] 
                 sim.se.pcr <- values.sim.se[1:nvariables]  
                 sim.far <- values.sim.mean[(nvariables+1):(2*nvariables)]
                 sim.se.far <- values.sim.se[(nvariables+1):(2*nvariables)]
                 sim.far.ci <- values.ci[(nvariables+1):(2*nvariables)]
                 sim.pc <- values.sim.mean[(2*nvariables+1):(3*nvariables)]
                 sim.pc.ci <- values.ci[(2*nvariables+1):(3*nvariables)]
                 sim.se.pc <- values.sim.se[(2*nvariables+1):(3*nvariables)]
                 sim.fa <- values.sim.mean[(3*nvariables+1):(4*nvariables)]
                 sim.fa.ci <- values.ci[(3*nvariables+1):(4*nvariables)]
                 sim.se.fa <- values.sim.se[(3*nvariables+1):(4*nvariables)]
                 pc.test <- which(!(valuesx > sim.pcr.ci))[1]-1  
                 fa.test <- which(!(fa.valuesx > sim.far.ci))[1]-1  } else { #not resampling, just sim
                   sim.pc <- values.sim.mean[1:nvariables] 
                   sim.pc.ci <- values.ci[1:nvariables] 
                   sim.se.pc <- values.sim.se[1:nvariables]  
                   sim.fa <- values.sim.mean[(nvariables+1):(2*nvariables)]
                   sim.fa.ci <- values.ci[(nvariables+1):(2*nvariables)]
                   sim.se.fa <- values.sim.se[(nvariables+1):(2*nvariables)]
                   pc.test <- which(!(valuesx > sim.pc.ci))[1]-1  
                   fa.test <- which(!(fa.valuesx > sim.fa.ci))[1]-1 }
               
               if (plot) { points(sim.pc,type ="l",lty="dotted",pch=4,col="red")
                 points(sim.fa,type ="l",lty="dotted",pch=4,col="red")
                 # sim.pcr <- sim.far <- NA   @#removed Dec 31, 2016
                 points(sim.pcr,type ="l",lty="dashed",pch=2,col="red")
                 points(sim.far,type ="l",lty="dashed",pch=2,col="red")}
               pc.test <- which(!(valuesx > sim.pc.ci))[1]-1 
               fa.test <- which(!(fa.valuesx > sim.fa.ci))[1]-1  
             } else {      #sim is false
               sim.pcr <- values.sim.mean[1:nvariables]
               sim.pcr.ci <- values.ci[1:nvariables]
               sim.se.pcr <- values.sim.se[1:nvariables]
               sim.far <- values.sim.mean[(nvariables+1):(2*nvariables)]
               sim.far.ci <- values.ci[(nvariables+1):(2*nvariables)]
               sim.se.far <- values.sim.se[(nvariables+1):(2*nvariables)]
               sim.fa <- NA
               sim.pc <- NA
               sim.se.fa <- NA
               sim.se.pc <- NA
               pc.test <- which(!(valuesx > sim.pcr.ci))[1]-1  
               fa.test <- which(!(fa.valuesx > sim.far.ci))[1]-1 } 
             if(resample) {     if (plot) {points(sim.pcr,type ="l",lty="dashed",pch=4,col="red")
               points(sim.far,type ="l",lty="dashed",pch=4,col="red") }
             } } 
    )  
    
    
    
    if(error.bars) {  
      if(!any(is.na(sim.pc))) {  
        for (i in 1:length(sim.pc))  {
          ycen <- sim.pc[i]
          yse <-  sim.se.pc[i]
          arrows(i,ycen-ci*yse,i,ycen+ci* yse,length=arrow.len, angle = 90, code=3,col = par("fg"), lty = NULL, lwd = par("lwd"), xpd = NULL)} 
      }
      if(!any(is.na(sim.pcr))) {    
        for (i in 1:length(sim.pcr))   {
          ycen <- sim.pcr[i]
          yse <-  sim.se.pcr[i]
          arrows(i,ycen-ci*yse,i,ycen+ci* yse,length=arrow.len, angle = 90, code=3,col = par("fg"), lty = NULL, lwd = par("lwd"), xpd = NULL)} 
      }
      if(!any(is.na(sim.fa))) {      	       
        for (i in 1:length(sim.fa))  {
          ycen <- sim.fa[i]
          yse <-  sim.se.fa[i]
          arrows(i,ycen-ci*yse,i,ycen+ci* yse,length=arrow.len, angle = 90, code=3,col = par("fg"), lty = NULL, lwd = par("lwd"), xpd = NULL)} 
      }
      if(!any(is.na(sim.far))) {      	
        for (i in 1:length(sim.far))  {
          ycen <- sim.far[i]
          yse <-  sim.se.far[i]
          arrows(i,ycen-ci*yse,i,ycen+ci* yse,length=arrow.len, angle = 90, code=3,col = par("fg"), lty = NULL, lwd = par("lwd"), xpd = NULL)} 
      } 
    }
    
    
    if(show.legend && plot) {
      if(is.null(n.obs)) {   #that is, do we have real data or a correlation matrix
        switch(fa,  
               both = {if(sim) {legend("topright", c("  PC  Actual Data", "  PC  Simulated Data", " PC  Resampled Data","  FA  Actual Data", "  FA  Simulated Data", " FA  Resampled Data"),
                                       col = c("blue","red","red","blue","red","red"),pch=c(4,NA,NA,2,NA,NA),
                                       text.col = "green4", lty = c("solid","dotted", "dashed","solid","dotted", "dashed"),
                                       merge = TRUE, bg = 'gray90')
               } else {legend("topright", c("  PC  Actual Data",  " PC  Resampled Data","  FA  Actual Data",  " FA  Resampled Data"),
                              col = c("blue","red","blue","red"),pch=c(4,NA,2,NA,NA),
                              text.col = "green4", lty = c("solid","dashed", "solid","dashed"),
                              merge = TRUE, bg = 'gray90')}},
               
               pc = {if(sim) {legend("topright", c("  PC  Actual Data", "  PC  Simulated Data", " PC  Resampled Data"), col = c("blue","red","red","blue","red","red"),pch=c(4,NA,NA,2,NA,NA),
                                     text.col = "green4", lty = c("solid","dotted", "dashed","solid","dotted", "dashed"),
                                     merge = TRUE, bg = 'gray90')} else {
                                       legend("topright", c("  PC  Actual Data",  " PC  Resampled Data"), col = c("blue","red","red","blue","red","red"),pch=c(4,NA,NA,2,NA,NA),
                                              text.col = "green4", lty = c("solid", "dashed","solid","dotted", "dashed"),
                                              merge = TRUE, bg = 'gray90')
                                     } } ,  
               
               fa = {if(sim) {legend("topright", c("  FA  Actual Data", "  FA  Simulated Data", " FA  Resampled Data"), col = c("blue","red","red","blue","red","red"),pch=c(4,NA,NA,2,NA,NA),
                                     text.col = "green4", lty = c("solid","dotted", "dashed","solid","dotted", "dashed"),
                                     merge = TRUE, bg = 'gray90')} else {
                                       legend("topright", c("  FA  Actual Data",  " FA  Resampled Data"), col = c("blue","red","red","blue","red","red"),pch=c(4,NA,NA,2,NA,NA),
                                              text.col = "green4", lty = c("solid", "dashed","solid","dotted", "dashed"),
                                              merge = TRUE, bg = 'gray90')
                                     }
               }   
        ) } else {
          switch(fa,
                 
                 both= {      legend("topright", c("PC  Actual Data", " PC  Simulated Data","FA  Actual Data", " FA  Simulated Data"), col = c("blue","red","blue","red"),pch=c(4,NA,2,NA),
                                     text.col = "green4", lty = c("solid","dotted","solid","dotted"),
                                     merge = TRUE, bg = 'gray90')},
                 
                 pc= {   legend("topright", c("PC  Actual Data", " PC  Simulated Data"), col = c("blue","red","blue","red"),pch=c(4,NA,2,NA),
                                text.col = "green4", lty = c("solid","dotted","solid","dotted"),
                                merge = TRUE, bg = 'gray90')},
                 fa =  {legend("topright", c("FA  Actual Data", " FA  Simulated Data"), col = c("blue","red","blue","red"),pch=c(4,NA,2,NA),
                               text.col = "green4", lty = c("solid","dotted","solid","dotted"),
                               merge = TRUE, bg = 'gray90')})}
    }
    
    
    colnames(values) <- paste0("Sim",1:ncol(values))
    if(fa!= "pc" && plot) abline(h=1)
    results <- list(fa.values = fa.valuesx,pc.values=valuesx,pc.sim=sim.pc,pc.simr = sim.pcr,fa.sim=sim.fa,fa.simr = sim.far,nfact=fa.test,ncomp=pc.test, Call=cl) 
    
    
    
    if (fa == "pc" )  {
      colnames(values)[1:nvariables] <- paste0("C",1:nvariables)
    } else {
      
      
      colnames(values)[1:(2*nvariables)] <- c(paste0("C",1:nvariables),paste0("F",1:nvariables))
      if(sim) {
        if(resample) colnames(values)[(2*nvariables +1):ncol(values)] <- c(paste0("CSim",1:nvariables),paste0("Fsim",1:nvariables))
      } 
      
      results$nfact <- fa.test}
    
    results$ncomp <- pc.test
    results$values <- values
    cat("Parallel analysis suggests that ")
    cat("the number of factors = ",fa.test, " and the number of components = ",pc.test,"\n")
    class(results) <- c("psych","parallel")
    return(invisible(results))
  }

assignInNamespace("fa.parallel", fa.parallel_fixed, pos = "package:psych")

