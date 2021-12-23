
mKrigMAGMA.MLE <- function(x, y, weights = rep(1, nrow(x)), 
                                Z = NULL, ..., par.grid = NULL, lambda = NULL, lambda.profile = TRUE, 
                                verbose = FALSE, relative.tolerance = 1e-04, 
                                nGPUs = 0, MAGMArank = 0, singlePrecision=FALSE) {
  
  #Get distance function and arguments if available.  Otherwise use 'dist' function
  #to compute upper triangle of distance matrix
  #
  cov.fun = list(...)$cov.fun
  cov.args= list(...)$cov.args
  onlyUpper = supportsOnlyUpper(cov.fun, cov.args)
  
  noDist = FALSE
  if(!is.null(list(...)$Distance))
    Dist.fun = list(...)$Distance
  else if(!is.null(list(...)$cov.args) && ! is.null(list(...)$cov.args$Distance))
    Dist.fun = list(...)$cov.args$Distance
  else
    Dist.fun = "dist"
  
  if(!is.null(list(...)$Dist.args))
    Dist.args = list(...)$Dist.args
  else if(!is.null(list(...)$cov.args) && !is.null(list(...)$cov.args$Dist.args))
    Dist.args = list(...)$cov.args&Dist.args
  else
    Dist.args = NULL
  
  #if it's optimal to use rdist (in the case that the covariance function supports 
  #evaluation over its upper triangle only), use it
  if(identical(match.fun(Dist.fun), dist) && onlyUpper && is.null(Dist.args))
    Dist.fun = "rdist"

  #precompute distance matrix so it only needs to be computed once
  distMat = do.call(Dist.fun, c(list(x), Dist.args))
  
  #If covariance function has support for evaluation over upper triangle of 
  #distance matrix only, convert distance matrix from compact to true matrix 
  #if necessary
  #
  if(onlyUpper && class(distMat) == 'dist')
    distMat = compactToMat(distMat)
  
  # these are all the arguments needed to call mKrigMAGMA except lambda
  mKrig.args <- c(list(x = x, y = y, weights = weights, Z = Z, distMat=distMat, 
                       nGPUs = nGPUs, MAGMArank = MAGMArank, singlePrecision=singlePrecision),
                  list(...))
  lnProfileLike.max <- -1e+20
  # find NG --  number of parameters to try
  par.grid <- data.frame(par.grid)
  if (nrow(par.grid) == 0) {
    if (is.null(lambda)) {
      NG <- 1
    }
    else {
      NG <- length(lambda)
    }
  }
  else {
    NG <- nrow(par.grid)
  }
  # output matrix to summarize results
  summary <- matrix(NA, nrow = NG, ncol = 8)
  dimnames(summary) <- list(NULL, c("EffDf", "lnProfLike", 
                                    "GCV", "sigma.MLE", "rho.MLE", "llambda.MLE", "counts eval", 
                                    "counts grad"))
  lambda.best <- NA
  # default for lambda is 1.0 for first value and exp(llambda.opt) for subsequent ones
  # this is controlled by NAs for lambda starting values.
  if (is.null(lambda)) {
    lambda <- rep(NA, NG)
  }
  # default starting value for lambda is 1 or log lambda is 0
  llambda.opt <- 0
  optim.counts <- c(NA, NA)
  lnLike.eval <- list()
  # Define the objective function as a tricksy call to mKrigMAGMA 
  # if Y is a matrix of replicated data sets use the log likelihood for the complete data sets
  temp.fn <- function(x) {
    # NOTE: FULL refers to estimates collapsed across the replicates if Y is a matrix
    # assign to hold only a few components returned by mKrigMAGMA
    hold <- do.call("mKrigMAGMA", c(mKrig.args, list(find.trA = FALSE), onlyUpper=TRUE, 
                                    list(lambda = exp(x)), cov.args.temp))[c("lambda.fixed", 
                                                                             "rho.MLE.FULL", "sigma.MLE.FULL", "lnProfileLike.FULL")]
    # add this evalution to an  object (i.e. here a matrix) in the calling frame
    temp.eval <- get("capture.evaluations")
    assign("capture.evaluations", rbind(temp.eval, unlist(hold)), 
           envir = capture.env)
    return(hold$lnProfileLike.FULL)
  }
  #
  # begin loop over covariance arguments
  for (k in 1:NG) {
    llambda.start <- ifelse(is.na(lambda[k]), llambda.opt, 
                            log(lambda[k]))
    # list of covariance arguments from par.grid with right names (some R arcania!)
    # note that this only works because 1) temp.fn will search in this frame for this object
    # par.grid has been coerced to a data frame so one has a concept of a row subscript.
    cov.args.temp <- as.list(par.grid[k, ])
    names(cov.args.temp) <- names(par.grid)
    if (lambda.profile) {
      # set up matrix to store evaluations from within optim
      capture.evaluations <- matrix(NA, ncol = 4, nrow = 1, 
                                    dimnames = list(NULL, c("lambda", "rho.MLE", 
                                                            "sigma.MLE", "lnProfileLike.FULL")))
      capture.env <- environment()
      # call to optim
      look <- optim(llambda.start, temp.fn, method = "BFGS", 
                    control = list(fnscale = -1, parscale = 0.1, 
                                   ndeps = 0.05, reltol = relative.tolerance))
      llambda.opt <- look$par
      optim.counts <- look$counts
      # call to 1-d search
      #            opt.summary     <- optimize(temp.fn, interval= llambda.start + c(-8,8), maximum=TRUE)
      #            llambda.opt <- opt.summary$maximum
      #            optim.counts<- c(nrow(capture.evaluations)-1, NA)
      # accumulate the new matrix of lnlambda and ln likelihoods (omitting first row of NAs)
      lnLike.eval <- c(lnLike.eval, list(capture.evaluations[-1, 
                                                             ]))
    }
    else {
      # no refinement for lambda so just save the the 'start' value as final one.
      llambda.opt <- llambda.start
    }
    # final fit at optimal value (or starting value if not refinement/maximization for lambda)
    obj <- do.call("mKrigMAGMA", c(mKrig.args, cov.args.temp, onlyUpper=TRUE, 
                                   list(lambda = exp(llambda.opt))))
    if (obj$lnProfileLike.FULL > lnProfileLike.max) {
      lnProfileLike.max <- obj$lnProfileLike.FULL
      cov.args.MLE <- cov.args.temp
      lambda.best <- exp(llambda.opt)
    }
    # save results of the kth covariance model evaluation
    summary[k, 1:8] <- c(obj$eff.df, obj$lnProfileLike.FULL, 
                         obj$GCV, obj$sigma.MLE.FULL, obj$rho.MLE.FULL, llambda.opt, 
                         optim.counts)
    if (verbose) {
      cat("Summary: ", k, summary[k, 1:8], fill = TRUE)
    }
  }
  return(list(summary = summary, par.grid = par.grid, cov.args.MLE = cov.args.MLE, 
              mKrig.args = list(...), lambda.best = lambda.best, lambda.MLE = lambda.best, 
              call = match.call(), lnLike.eval = lnLike.eval))
}