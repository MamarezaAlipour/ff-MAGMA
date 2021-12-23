
dist = function(x1, x2 = NULL, ...) {
  #If only one argument is input, return compact distance matrix using default dist
  if(is.null(x2))
    return(do.call(stats::dist, c(list(x1), list(...))))
  
  #otherwise, if two arguments are input, use rdist
  return(rdist(x1, x2))
}