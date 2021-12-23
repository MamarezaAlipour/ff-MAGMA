
mKrigMAGMA <- function(x, y, weights = rep(1, nrow(x)), 
    Z = NULL, lambda = 0, cov.function = "stationary.fast.cov", m = 2, 
    chol.args = NULL, cov.args = NULL, find.trA = TRUE, NtrA = 20, 
    iseed = 123, llambda = NULL,
    nGPUs = 0, MAGMArank = 0, singlePrecision=FALSE, ...) {
    
    #nGPUs: number of GPUs to use on this problem (only 0, 1, or 2 works)
    #MAGMArank: if nGPUs is 1, this determines which GPU to use (from 0 to n
    #assuming n GPUs in total).
    doMAGMA = FALSE
    if(nGPUs < 0)
      stop('Must use nonnegative number of GPUs')
    else if(nGPUs >= 1) { #user wants to use MAGMA 
      doMAGMA = TRUE
    }
    
    # grab other arguments for covariance function
    cov.args <- c(cov.args, list(...))
    #
    #If cov.args$find.trA is true, set onlyUpper to FALSE (onlyUpper doesn't
    #play nice with predict.mKrig, called by mKrig.trace)
    #
    if(find.trA == TRUE && supportsOnlyUpper(cov.function, cov.args))
      cov.args$onlyUpper= FALSE
    
    if (!is.null(llambda)) {
        lambda <- exp(llambda)
    }
    # see comments in Krig.engine.fixed for algorithmic commentary
    #
    # check for duplicate x's.
    # stop if there are any
    if (any(duplicated(cat.matrix(x)))) 
        stop("locations are not unique see help(mKrigMAGMA) ")
    if (any(is.na(y))) 
        stop("Missing values in y should be removed")
    if (!is.null(Z)) {
        Z <- as.matrix(Z)
    }
    
    # create fixed part of model as m-1 order polynomial
    Tmatrix <- cbind(fields.mkpoly(x, m), Z)
    # set some dimensions
    np <- nrow(x)
    nt <- ncol(Tmatrix)
    nZ <- ifelse(is.null(Z), 0, ncol(Z))
    ind.drift <- c(rep(TRUE, (nt - nZ)), rep(FALSE, nZ)) 
    # as a place holder for reduced rank Kriging, distinguish between
    # observations locations and  the locations to evaluate covariance.
    # (this is will also allow predict.mKrigMAGMA to handle a Krig object)
    knots <- x
    # covariance matrix at observation locations
    # NOTE: if cov.function is a sparse constuct then Mc will be sparse.
    # see e.g. wendland.cov
    Mc <- do.call(cov.function, c(cov.args, list(x1 = x, x2 = x)))
    #
    # decide how to handle the pivoting.
    # one wants to do pivoting if the matrix is sparse.
    # if Mc is not a matrix assume that it is in sparse format.
    #
    sparse.flag <- !is.matrix(Mc)
    #
    # set arguments that are passed to cholesky
    #
    if (is.null(chol.args)) {
        chol.args <- list(pivot = sparse.flag)
    }
    else {
        chol.args <- chol.args
    }
    # quantify sparsity of Mc for the mKrig object
    nzero <- ifelse(sparse.flag, length(Mc@entries), np^2)
    # add diagonal matrix that is the observation error Variance
    # NOTE: diag must be a overloaded function to handle  sparse format.
    if (lambda != 0) {
        diag(Mc) <- (lambda/weights) + diag(Mc)
    }
    # At this point Mc is proportional to the covariance matrix of the
    # observation vector, y.
    #
    # cholesky decoposition of Mc
    if(! doMAGMA) {
        # do.call used to supply other arguments to the function
        # especially for sparse applications.
        # If chol.args is NULL then this is the same as
        #              Mc<-chol(Mc), chol.args))
        Mc <- do.call("chol", c(list(x = Mc), chol.args))
    } else {
        # Use MAGMA to perform Cholesky decomposition
        # Check to make sure no unused arguments are passed:
        if(! is.null(chol.args$pivot) && chol.args$pivot == TRUE)
            stop('MAGMA does not support pivoting')
        if(! is.null(chol.args$LINPACK) && chol.args$LINPACK == TRUE)
            stop('LINPACK is deprecated, and cannot be used concurrently with MAGMA')
        
        # Call MAGMA to perform Cholesky decomposition on Mc
        if(!singlePrecision) {
          if(nGPUs == 1)
            .Call("magmaCholeskyFinal",Mc,nrow(Mc),0,MAGMArank,0,0)
          else
            .Call("magmaCholeskyFinal_m",Mc,nrow(Mc),0,0,nGPUs,0)
        } else {
          if(nGPUs == 1)
            .Call("smagmaCholeskyFinal",Mc,nrow(Mc),0,MAGMArank,0,0)
          else
            .Call("smagmaCholeskyFinal_m",Mc,nrow(Mc),0,0,nGPUs,0)
        }
        
        #magmaChol(Mc, MAGMArank=MAGMArank, nGPUs=nGPUs, deepCopy=FALSE, 
        #          zeroTri=FALSE, singlePrecision=singlePrecision)
    }
    
    
    lnDetCov <- 2 * sum(log(diag(Mc)))
    
    # Efficent way to multply inverse of Mc times the Tmatrix
    VT <- forwardsolve(Mc, x = Tmatrix, transpose = TRUE, upper.tri = TRUE)
    qr.VT <- qr(VT)
    # start linear algebra to find solution
    # Note that all these expressions make sense if y is a matrix
    # of several data sets and one is solving for the coefficients
    # of all of these at once. In this case d.coef and c.coef are matrices
    #
    # now do generalized least squares for d
    d.coef <- as.matrix(qr.coef(qr.VT, forwardsolve(Mc, transpose = TRUE, 
        y, upper.tri = TRUE)))
    # and then find c.
    # find the coefficents for the spatial part.
    c.coef <- as.matrix(forwardsolve(Mc, transpose = TRUE, y - 
        Tmatrix %*% d.coef, upper.tri = TRUE))
    # save intermediate result this is   t(y- T d.coef)( M^{-1}) ( y- T d.coef)
    quad.form <- c(colSums(as.matrix(c.coef^2)))
    # find c coefficients
    c.coef <- as.matrix(backsolve(Mc, c.coef))
    # GLS covariance matrix for fixed part.
    Rinv <- solve(qr.R(qr.VT))
    Omega <- Rinv %*% t(Rinv)
    # MLE estimate of rho and sigma
    #    rhohat <- c(colSums(as.matrix(c.coef * y)))/(np - nt)
    # NOTE if y is a matrix then each of these are vectors of parameters.
    rho.MLE <- quad.form/np
    rhohat <- c(colSums(as.matrix(c.coef * y)))/np
    shat.MLE <- sigma.MLE <- sqrt(lambda * rho.MLE)
    # the  log profile likehood with  rhohat  and  dhat substituted
    # leaving a profile for just lambda.
    # NOTE if y is a matrix then each of this is a vector of log profile
    # likelihood values.
    lnProfileLike <- (-np/2 - log(2 * pi) * (np/2) - (np/2) * 
        log(rho.MLE) - (1/2) * lnDetCov)
    rho.MLE.FULL <- mean(rho.MLE)
    sigma.MLE.FULL <- sqrt(lambda * rho.MLE.FULL)
    # if y is a matrix then compute the combined likelihood
    # under the assumption that the columns of y are replicated
    # fields
    lnProfileLike.FULL <- sum((-np/2 - log(2 * pi) * (np/2) - 
        (np/2) * log(rho.MLE.FULL) - (1/2) * lnDetCov))
    #
    # return coefficients and   include lambda as a check because
    # results are meaningless for other values of lambda
    # returned list is an 'object' of class mKrig (micro Krig)
    # also save the matrix decompositions so coefficients can be
    # recalculated for new y values.
    out <- list(d = (d.coef), c = (c.coef), nt = nt, np = np, 
        lambda.fixed = lambda, x = x, y=y, weights=weights, knots = knots,
        cov.function.name = cov.function, 
        args = cov.args, m = m, chol.args = chol.args, call = match.call(), 
        nonzero.entries = nzero, shat.MLE = sigma.MLE, sigma.MLE = sigma.MLE, 
        rho.MLE = rho.MLE, rhohat = rho.MLE, lnProfileLike = lnProfileLike, 
        rho.MLE.FULL = rho.MLE.FULL, sigma.MLE.FULL = sigma.MLE.FULL, 
        lnProfileLike.FULL = lnProfileLike.FULL, lnDetCov = lnDetCov, 
        quad.form = quad.form, Omega = Omega, qr.VT = qr.VT, 
        Mc = Mc, Tmatrix = Tmatrix, ind.drift = ind.drift, nZ = nZ)
    #
    # find the residuals directly from solution
    # to avoid a call to predict
    out$residuals <- lambda * c.coef/weights
    out$fitted.values <- y - out$residuals
    # estimate effective degrees of freedom using Monte Carlo trace method.
    if (find.trA) {
        out2 <- mKrig.trace(out, iseed, NtrA)
        out$eff.df <- out2$eff.df
        out$trA.info <- out2$trA.info
        out$GCV <- (sum(out$residuals^2)/np)/(1 - out2$eff.df/np)^2
        if (NtrA < np) {
            out$GCV.info <- (sum(out$residuals^2)/np)/(1 - out2$trA.info/np)^2
        }
        else {
            out$GCV.info <- NA
        }
    }
    else {
        out$eff.df <- NA
        out$trA.info <- NA
        out$GCV <- NA
    }
    class(out) <- "mKrig"
    return(out)
}
