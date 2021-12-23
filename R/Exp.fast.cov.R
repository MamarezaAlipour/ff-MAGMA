
"Exp.fast.cov" <- function(x1, x2=NULL, theta = rep(1, ncol(x1)), 
                      distMat = NA, C = NA, marginal = FALSE, onlyUpper=FALSE) {
  
  if (!is.matrix(x1)) 
    x1 <- as.matrix(x1)
  if (is.null(x2)) 
    x2 <- x1
  if (!is.matrix(x2)) 
    x2 <- as.matrix(x2)
  if (length(theta) == 1) 
    theta <- rep(theta, ncol(x1))
  d <- ncol(x1)
  n1 <- nrow(x1)
  n2 <- nrow(x2)
  # scale the coordinates by theta
  # a more general scaling by a matrix is done in stationary.cov
  x1 <- scale(x1, center = FALSE, scale = theta)
  x2 <- scale(x2, center = FALSE, scale = theta)
  #
  # there are five possible actions listed below:
  # if no cross covariance matrix and marginal variance not desired
  if (is.na(C[1]) && !marginal) {
    #
    # if distMat is supplied, evaluate covariance for upper triangular part only
    if(!is.na(distMat[1]) && class(distMat) == "dist") {
      #distMat is in compact form, so evaluate over all distMat and convert to matrix form
      
      return(compactToMat(exp(-distMat), diagVal=1))
    } else if(!is.na(distMat[1]) && onlyUpper && nrow(distMat) == ncol(distMat)) {
      #distMat is an actual matrix, so evaluate covariance over upper triangle only
      
      return(ExponentialUpper(distMat, range=theta))
    }
    else
      #distMat not supplied so must compute it along with covariance matrix
      return(exp(-rdist(x1, x2)))
  }
  #
  # multiply cross covariance matrix by C
  # in this case implemented in FORTRAN
  #
  if (!is.na(C[1])) {
#     return(.Fortran("multeb", PACKAGE="fields",
#                     nd = as.integer(d),
#                     x1 = as.double(x1), 
#                     n1 = as.integer(n1),
#                     x2 = as.double(x2),
#                     n2 = as.integer(n2), 
#                     par = as.double(1),
#                     c = as.double(C),
#                     h = as.double(rep(0, n1)),
#                     work = as.double(rep(0, n2)))$h)
    return(.Call("multeb", 
                 nd = as.integer(d),
                 x1 = as.double(x1), 
                 n1 = as.integer(n1),
                 x2 = as.double(x2),
                 n2 = as.integer(n2), 
                 par = as.double(1),
                 c = as.double(C),
                 work = as.double(rep(0, n2))))
  }
  #
  # return marginal variance ( 1.0 in this case)
  if (marginal) {
    return(rep(1, nrow(x1)))
  }
}
