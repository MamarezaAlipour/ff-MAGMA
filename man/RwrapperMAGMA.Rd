\name{RwrapperMAGMA}
\alias{magmaChol}

\title{
%%  ~~function to do ... ~~
MAGMA-accelerated Cholesky decomposition
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
\code{magmaChol} performs a Cholesky decomposition using the 
  MAGMA multi-GPU, multi-CPU library.  The user can specify 
  several options to customize the performance of the function 
  to increase computation speed or to suite the user's needs.
}
\usage{
magmaChol(A, nGPUs = 1, MAGMArank = 0, deepCopy = TRUE, 
  zeroTri = TRUE, lower.tri = FALSE, singlePrecision=FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{A}{
The matrix to perform the Cholesky decomposition on.  Note that 
  A must be real symmetric positive-definite.
}
  \item{nGPUs}{
The number of GPUs to use with MAGMA
}
  \item{MAGMArank}{
The GPU ID (rank) of the GPU to use with MAGMA.  This is 
  only used if \code{nGPUs = 1}.  If \code{nGPUs = n}, then 
  ranks 0 to n-1 are used.
}
  \item{deepCopy}{
If \code{FALSE}, the input matrix \code{A} is overwritten with 
  the result of the Cholesky decomposition to save time and 
  space.  Otherwise \code{A} is copied and the returned matrix 
  is an overwritten copy of \code{A}.
}
  \item{zeroTri}{
If \code{TRUE}, the returned matrix sets the correct triangle 
  of the matrix to zero.  Which triangle is set to zero 
  depends on which factor of the Cholesky decomposition 
  \code{magmaChol} is set to return.
}
  \item{lower.tri}{
If \code{TRUE}, the lower triangular Cholesky factor is 
  returned.  Otherwise, the upper triangular factor is 
  returned.
}
  \item{singlePrecision}{
If \code{TRUE}, Cholesky decomposition is performed in MAGMA
  using single precision.  Otherwise uses MAGMA's double 
  precision Cholesky decompostion.  Using single precision 
  will increase computation speed but reduce accuracy.
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
The Cholesky decompsition of a real symmetric positive-
  definite matrix $A$ is $U^T U$, where $U$ is upper-
  triangular (including a diagonal).  \code{magmaChol} 
  by default is used exactly the same as the \code{chol} 
  function.  While \code{magmaChol} is more customizable 
  than \code{chol} in many ways, it does not support 
  pivoting like \code{chol}.
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
The Cholesky factor matrix resulting from performing 
  the Cholesky decomposition on \code{A}.
}
\references{
%% ~put references to the literature/web site here ~
http://icl.cs.utk.edu/magma/
}
\author{
%%  ~~who you are~~
John Paige
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
\code{\link{chol}}
}
\examples{
#generate a covariance matrix that is real symmetric positive-definite
n = 5
x = matrix(runif(2*n), ncol=2)
distMat = rdist(x, x)
covMat = exp(-distMat/.2)
covMat = covMat + diag(.1, nrow=n)
U = magmaChol(covMat, nGPUs=1)

#compare covariance matrix to reconstructed covariance matrix using
#Cholesky factor
covMat
t(U) %*% U

#Compare different uses of the magmaChol function
U
magmaChol(covMat, nGPUs=1, zeroTri=FALSE)
magmaChol(covMat, nGPUs=1, lower.tri=TRUE)
magmaChol(covMat, nGPUs=1, singlePrecision=TRUE)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ Cholesky }
\keyword{ Choleski }
\keyword{ GPU }
\keyword{ MAGMA }% __ONLY ONE__ keyword per line
