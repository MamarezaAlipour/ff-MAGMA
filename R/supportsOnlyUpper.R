supportsOnlyUpper = function(cov.fun, cov.args=NULL) {
  
  goodCovFun = function(fun) {
    identical(match.fun(fun), actualCovFun)
  }
  
  
  #List of supported covariance functions:
  supported.cov.fun = c(stationary.fast.cov, Exp.fast.cov)
  supportedStationaryCovFuns = c(Exponential) #soon "Matern" will also be supported
  
  if(is.null(cov.fun)) {
    #set cov.fun to the default covariance function if not specified
    cov.fun = stationary.fast.cov
  }
  
  #cov.fun must be stationary.fast.cov or Exp.fast.cov
  actualCovFun = match.fun(cov.fun)
  if(!any(sapply(supported.cov.fun, goodCovFun)))
    return(FALSE)
  
  else {
    #if cov.fun is stationary.fast.cov, or exponential.fast.cov. Check to see if cov.args$Covariance 
    #is a supported function, if it exists.
    #
    
    #stationary.fast.cov by default uses Exponential covariance, which is supported.  Also other 
    #supported covariance functions don't have cov.args$Covariance, so return true in that case
    #as well.
    #
    if(is.null(cov.args) || is.null(cov.args$Covariance))
      return(TRUE)
    
    #check if covariance function used in stationary.fast.cov is supported
    actualCovFun = match.fun(cov.args$Covariance)
    return(any(sapply(supportedStationaryCovFuns, goodCovFun)))
    
  }
}