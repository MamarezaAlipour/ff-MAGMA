
compactToMat = function(compactMat, diagVal=0, lower.tri=FALSE) {
  #compactMat: a symmetric matrix stored as a vector containing elements for the upper triangle
  #portion of the true matrix
  #diagVal: a number to put in the diagonal entries of the output matrix
  #lower.tri: if TRUE, returns a lower-triangular matrix.  Otherwise returns an upper-
  #triangular matrix.
  
  storage.mode(compactMat) <- "double"
  
  if(class(compactMat) == 'dist') {
    n = attr(compactMat, "Size")
  } else { # (n^2 - n)/2 = length(compactMat)
    stop("input matrix is not compact or is not of class \"dist\"")
    
    #or if class is not dist but input matrix is still compact, use:
    #n = (1 + sqrt(1 + 8*length(compactMat)))/2
  }
  
  return(.Call("compactToMatC", compactMat, length(compactMat), as.integer(n), as.double(diagVal), 
               as.integer(lower.tri)))
}
