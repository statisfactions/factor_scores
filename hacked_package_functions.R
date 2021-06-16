## psych -----
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
if (("package:psych" %in% search())) {
  assignInNamespace("fa.parallel", fa.parallel_fixed, pos = "package:psych")
}
## semTools -----
### Sunthud Pornprasertmanit & Terrence D. Jorgensen
### Last updated: 27 May 2020
### fit and rotate EFA models in lavaan

efaUnrotate_hacked <- function(data = NULL, nf, varList = NULL, aux = NULL, ...) {
  "Hacked version"
  start = TRUE
  efaArgs <- list(...)
  if (is.null(data)) {
    ## check for summary statistics
    sample.cov  <- efaArgs$sample.cov
    sample.mean <- efaArgs$sample.mean
    sample.nobs <- efaArgs$sample.nobs
    sample.th   <- efaArgs$sample.th
    WLS.V       <- efaArgs$WLS.V
    NACOV       <- efaArgs$NACOV
    
    if (is.null(sample.cov)) stop('User must supply either raw data or ',
                                  'summary statistics to pass to lavaan().')
    if (is.null(varList)) varList <- colnames(sample.cov)
    
    anyOrdered <- !is.null(sample.th)
    ordNames <- efaArgs$ordered
    if (anyOrdered & (is.logical(ordNames) | is.null(ordNames))) {
      if (is.null(ordNames)) {
        message('Thresholds supplied, but not an ordered= argument. Must ',
                'assume all model variables are ordered.')
      }
      ordNames <- varList
    }
    
  } else {
    sample.cov  <- NULL
    sample.mean <- NULL
    sample.nobs <- NULL
    sample.th   <- NULL
    WLS.V       <- NULL
    NACOV       <- NULL
    
    efaArgs$data <- data
    if (is.null(varList)) varList <- colnames(data)
    
    anyOrdered <- semTools:::checkOrdered(data, varList, ...)
    ordNames <- semTools:::checkOrdered(data, varList, ..., return.names = TRUE)
    
    if (!is.null(efaArgs$group)) stop("Multi-group EFA is not currently supported.")
  }
  
  if (!is.null(aux)) {
    if (anyOrdered) {
      stop("The analysis model or the analysis data have ordered categorical",
           " variables. The auxiliary variable feature is not available for",
           " the models for categorical variables with the weighted least",
           " square approach.")
    }
    efaArgs$fixed.x <- FALSE
    efaArgs$missing <- "fiml"
    efaArgs$aux <- aux
    lavaancfa <- function(...) { cfa.auxiliary(...)}
  } else lavaancfa <- function(...) { lavaan::cfa(...)}
  nvar <- length(varList)
  facnames <- paste0("factor", 1:nf)
  loading <- outer(1:nvar, 1:nf, function(x, y) paste0("load", x, "_", y))
  syntax <- ""
  for (i in 1:nf) {
    variablesyntax <- paste(paste0(loading[,i], "*", varList), collapse = " + ")
    factorsyntax <- paste0(facnames[i], " =~ NA*", varList[1], " + ", variablesyntax, "\n")
    syntax <- paste(syntax, factorsyntax)
  }
  syntax <- paste(syntax, paste(paste0(facnames, " ~~ 1*", facnames),
                                collapse = "\n"), "\n")
  
  dataSupportsMeans <- length(setdiff(varList, ordNames)) && !(is.null(data) && is.null(sample.mean))
  meanstructure <- efaArgs$meanstructure
  if (is.null(meanstructure)) meanstructure <- anyOrdered #FIXME: wise default for EFA?
  stopifnot(is.logical(meanstructure))
  if (meanstructure && dataSupportsMeans) {
    syntax <- paste(syntax, paste(paste0(setdiff(varList, ordNames),
                                         " ~ 1"), collapse = "\n"), "\n")
  }
  
  if (nf > 1) {
    covsyntax <- outer(facnames, facnames,
                       function(x, y) paste0(x, " ~~ 0*", y, "\n"))[lower.tri(diag(nf), diag = FALSE)]
    syntax <- paste(syntax, paste(covsyntax, collapse = " "))
    for (i in 2:nf) {
      for (j in 1:(i - 1)) {
        loadconstraint <- paste(paste0(loading[,i], "*", loading[,j]), collapse=" + ")
        syntax <- paste(syntax, paste0("0 == ", loadconstraint), "\n")
      }
    }
  }
  
  List <- c(list(model = syntax, data = data), list(...))
  List$do.fit <- FALSE
  outtemp <- do.call(lavaancfa, List)
  covtemp <- lavInspect(outtemp, "sampstat")$cov
  partemp <- parTable(outtemp)
  err <- try(startload <- factanal(factors = nf, covmat = covtemp)$loadings[],
             silent = TRUE)
  if (is(err, "try-error")) {
    warning("The starting values from the factanal",
            " function cannot be calculated. Defaulting to  ",
            " start = FALSE instead.")
    startval <- sqrt(diag(diag(covtemp))) %*% startload
    partemp$ustart[match(as.vector(loading), partemp$label)] <- as.vector(startval)
    partemp$est <- partemp$se <- NULL
    syntax <- partemp
  } else {
    ## FIXME: optimizer can't escape loadings == 0 without starting values from
    ##        factanal() or by using theta parameterization
    ##        https://groups.google.com/d/msg/lavaan/ujkHmCVirEY/-LGut4ewAwAJ
    parameterization <- efaArgs$parameterization
    if (is.null(parameterization)) parameterization <- lavaan::lavOptions("parameterization")
    if (anyOrdered && parameterization != "theta")
      warning('If the default parameterization = "delta" returns results with ',
              'all factor loadings equal to zero, try either setting start ',
              '= TRUE or setting parameterization = "theta" instead.')
  }
  efaArgs$model <- syntax
  do.call(lavaancfa, efaArgs)
}

oblqRotate_hacked <- function(object, method = "quartimin", ...) {
  "Hacked version"
  requireNamespace("GPArotation")
  if (!("package:GPArotation" %in% search())) attachNamespace("GPArotation")
  mc <- match.call()
  initL <- semTools:::getLoad(object)
  if(ncol(initL) == 1) {
    phi <- matrix(1)
    lv.names <- colnames(initL)
    dimnames(phi) <- list(lv.names, lv.names)
    class(initL) = "matrix"
    return(new("EFA", 
        loading= initL, 
        rotate = matrix(NA), 
        gradRotate = matrix(NA),
        convergence = NA, 
        phi = phi, 
        se = matrix(NA), 
        method = NA_character_, 
        call = mc))
    
  } else {
  rotated <- GPArotation::GPFoblq(initL, method = method, ...)
  # rotateMat <- t(solve(rotated$Th)) # defined but never used
  
  phi <- rotated$Phi
  lv.names <- colnames(rotated$loading)
  dimnames(phi) <- list(lv.names, lv.names)
  return(new("EFA", 
      loading=rotated$loadings, 
      rotate = rotated$Th, 
      gradRotate = rotated$Gq,
      convergence = rotated$convergence, 
      phi = phi, 
      se = matrix(NA), 
      method = rotated$method, 
      call = mc))
  }
}
if (("package:semTools" %in% search())) {
  assignInNamespace("efaUnrotate", efaUnrotate_hacked, pos = "package:semTools")
  assignInNamespace("oblqRotate", oblqRotate_hacked, pos = "package:semTools")
}
