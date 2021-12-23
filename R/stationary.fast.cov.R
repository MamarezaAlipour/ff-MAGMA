"stationary.fast.cov" <- function(x1, x2=NULL, Covariance = "Exponential", Distance = "dist", 
                                  Dist.args = NULL, theta = 1, V = NULL, C = NA, marginal = FALSE, 
                                  derivative = 0, distMat = NA, onlyUpper = FALSE, ...) {
  
  # get covariance function arguments from call
  cov.args <- list(...)
  # coerce x1 and x2 to matrices
  if (is.data.frame(x1)) 
    x1 <- as.matrix(x1)
  if (!is.matrix(x1)) 
    x1 <- matrix(c(x1), ncol = 1)
  
  if(is.null(x2))
    x2 = x1
  if (is.data.frame(x2)) 
    x2 <- as.matrix(x2)
  if (!is.matrix(x2)& !is.null(x2)) 
    x2 <- matrix(c(x2), ncol = 1)
  d <- ncol(x1)
  n1 <- nrow(x1)
  n2 <- nrow(x2)
  #
  # separate out a single scalar transformation and a
  # more complicated scaling and rotation.
  # this is done partly to have the use of great circle distance make sense
  # by applying the scaling  _after_ finding the distance.
  #
  # tranform coordinates if theta not a scalar
  if (length(theta) > 1) {
    stop("theta as a vector matrix has been depreciated use the V argument")
  }
  #
  # following now treats V as a full matrix for scaling and rotation.
  #
  # try to catch incorrect conbination  of great circle distance and V
  if (Distance == "rdist.earth" & !is.null(V)) {
    stop("V not supported with great circle distance")
  }
  if (!is.null(V)) {
    if (theta != 1) {
      stop("can't specify both theta and V!")
    }
    x1 <- x1 %*% t(solve(V))
    x2 <- x2 %*% t(solve(V))
  }
  #
  # locations are now scaled and rotated correctly
  # now apply covariance function to pairwise distance matrix, or multiply
  # by C vector or just find marginal variance
  # this if block finds the cross covariance matrix
  if (is.na(C[1]) & !marginal) {
    #
    # if distMat is supplied, evaluate covariance for upper triangular part only
    #
    if(is.na(distMat[1])) {
      # distMat not supplied so must compute it along with covariance matrix
      # note overall scalling by theta (which is just theta under isotropic case)
      if(is.null(x2))
        distMat <- do.call(Distance, c(list(x1), Dist.args))
      else
        distMat <- do.call(Distance, c(list(x1=x1, x2=x2), Dist.args))
      
    }
    
    #
    # now convert distance matrix to covariance matrix
    #
    
    if(class(distMat) == "dist") {
      #
      # distMat is in compact form, so evaluate covariance over all distMat and convert to matrix form
      
      diagVal = do.call(Covariance, c(list(d=0), cov.args))
      
      if(onlyUpper)
        return(compactToMat(do.call(Covariance, c(list(d=distMat/theta), cov.args)), diagVal))
      else {
        #NOTE: this is a very inefficient way to fill in the lower triangle of the output matrix
        warning(paste0('For best performance, onlyUpper should be set to TRUE for ', "'", 
                       "dist", "'",  ' distance function and compact distance matrices.'))
        #
        # There is probably a more efficient way to perform the below operation.  Ideally, 
        # we would want to create a new class or use a symmetric matrix class in the Matrix
        # package.
        #
        
        ans = compactToMat(do.call(Covariance, c(list(d=distMat/theta), cov.args)), diagVal)
        ans[lower.tri(ans)] = t(ans[upper.tri(ans)])
        return(ans)
      }
    } else if(onlyUpper && nrow(distMat) == ncol(distMat)) {
      #
      #Eventually, it would be nice to create a CovarianceUpper function to evaluate the
      #input covariance function, but only over the upper triangle of the matrix like so:
      # return(CovarianceUpper(distMat, Covariance, cov.args))
      # 
      
      if(supportsOnlyUpper(stationary.fast.cov, list(Covariance=Covariance)))
        return(do.call("ExponentialUpper", c(list(distMat=distMat/theta), cov.args)))
      
      warning(paste0("onlyUpper option for more efficient code not yet supported for non-compact matrices with 
                     the given covariance.  Use 'supportsOnlyUpper' to see which covariance functions support 
                     the onlyUpper option."))
      
      return(do.call(Covariance, c(list(d = distMat/theta), cov.args)))
    } else {
      return(do.call(Covariance, c(list(d = distMat/theta), cov.args)))
    }
  }
  # or multiply cross covariance by C
  # as coded below this is not particularly efficient of memory
  if (!is.na(C[1])) {
    if(onlyUpper) {
      #the onlyUpper option is not compatible with the C option
      onlyUpper = FALSE
      
    }
    
    if(is.null(x2))
      bigD <- do.call(Distance, c(list(x1, x1), Dist.args))
    else
      bigD <- do.call(Distance, c(list(x1=x1, x2=x2), Dist.args))
    
    if (derivative == 0) {
      return(do.call(Covariance, c(list(d = bigD/theta), cov.args)) %*% C)
    }
    else {
      # find partial derivatives
      tempW <- do.call(Covariance, c(list(d = bigD/theta), 
                                     cov.args, derivative = derivative))
      # loop over dimensions and knock out each partial accumulate these in
      # in temp
      temp <- matrix(NA, ncol = d, nrow = n1)
      for (kd in 1:d) {
        # Be careful if the distance (tempD) is close to zero.
        # Note that the x1 and x2 are in transformed ( V inverse) scale
        sM <- ifelse(bigD == 0, 0, (tempW * (outer(x1[, kd], x2[, kd], "-"))/(theta * bigD)))
        # accumlate the new partial
        temp[, kd] <- sM %*% C
      }
      # transform back to original coordinates.
      if (!is.null(V)) {
        temp <- temp %*% t(solve(V))
      }
      return(temp)
    }
  }
  # or find marginal variance and return  a vector.
  if (marginal) {
    sigma2 <- do.call(Covariance, c(list(d = 0), cov.args))
    return(rep(sigma2, nrow(x1)))
  }
  
  # should not get here based on sequence of conditional if statements above.
}
