
ExponentialUpper = function(distMat, range = 1, alpha = 1/range, phi = 1) {
  # Evaluates the exponential covariance function over the upper triangle of the distance matrix
  
  if(nrow(distMat) != ncol(distMat))
    stop('distance matrix is non-symmetric.  Should not be calling ExponentialUpper.')
  
  #ans = .Fortran("ExponentialUpper", as.double(distMat), as.integer(nrow(distMat)), as.double(alpha), as.double(phi))
  return(.Call("ExponentialUpperC", as.double(distMat), as.integer(nrow(distMat)), as.double(alpha), as.double(phi)))
  
  #convert ans to standard matrix
  #ans = ans[[1]]
  #dim(ans) = dim(distMat)
  
  #return(ans)
}