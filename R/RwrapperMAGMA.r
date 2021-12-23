magmaChol = function(A, nGPUs = 1, MAGMArank = 0, deepCopy = TRUE, zeroTri = TRUE, lower.tri = FALSE, singlePrecision=FALSE) {

  if(!is.matrix(A)){
    stop("Input must be a matrix")
  }
  n = dim(A)
  if(n[1] != n[2]){
    stop("Matrix is not square")
  }
  
  #nb paramter not used at this point but needed for function call
  nb = 0

  #ensure deep copy of A:
  if(deepCopy) {
    A[1] = A[1]
  }
  
  #if nGPUs == 0, use default chol
  if(nGPUs == 0) {
    warning("since nGPUs == 0, using default `chol' function for Cholesky decomposition")
    return(chol(A))
  }
  
  if(!singlePrecision) {
    if(nGPUs == 1)
      .Call("magmaCholeskyFinal",A,n[1],nb,MAGMArank,as.integer(zeroTri),as.integer(lower.tri))
    else
      .Call("magmaCholeskyFinal_m",A,n[1],nb,as.integer(zeroTri),nGPUs,as.integer(lower.tri))
  } else {
    if(nGPUs == 1)
      .Call("smagmaCholeskyFinal",A,n[1],nb,MAGMArank,as.integer(zeroTri),as.integer(lower.tri))
    else
      .Call("smagmaCholeskyFinal_m",A,n[1],nb,as.integer(zeroTri),nGPUs,as.integer(lower.tri))
  }

  #only return result if A was deep-copied:
  if(deepCopy) {
    return(A)
  }
  else {
    return(invisible(NULL))
  }

}
