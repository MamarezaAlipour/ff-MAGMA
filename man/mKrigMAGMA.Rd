\name{mKrigMAGMA}
\alias{mKrigMAGMA}

\title{
%%  ~~function to do ... ~~
Fast Kriging with MAGMA
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
This function is an accelerated version of \code{mKrig} from \code{fields}. 
  It is able to use the MAGMA library and thereby use GPUs to accelerate 
  the Cholesky decomposition required in \code{mKrig}.
}
\usage{
mKrigMAGMA(x, y, weights = rep(1, nrow(x)), Z = NULL, lambda = 0,
  cov.function = "stationary.fast.cov", m = 2, chol.args =
  NULL, cov.args = NULL, find.trA = TRUE, NtrA = 20,
  iseed = 123, llambda = NULL, nGPUs = 0, MAGMArank = 0,
  singlePrecision = FALSE, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
%%     ~~Describe \code{x} here~~
Matrix of unique spatial locations (or in print or surface 
  the returned mKrig object.)
}
  \item{y}{
%%     ~~Describe \code{y} here~~
Vector or matrix of observations at spatial locations, 
  missing values are not allowed! Or in mKrig.coef a new 
  vector of observations. If y is a matrix the columns are 
  assumed to be independent observations vectors generated 
  from the same covariance and measurment error model.
}
  \item{weights}{
%%     ~~Describe \code{weights} here~~
Precision ( 1/variance) of each observation
}
  \item{Z}{
%%     ~~Describe \code{Z} here~~
Linear covariates to be included in fixed part of the 
  model that are distinct from the default low order 
  polynomial in \code{x}
}
  \item{lambda}{
%%     ~~Describe \code{lambda} here~~
Smoothing parameter or equivalently the ratio between 
  the nugget and process varainces.
}
  \item{cov.function}{
%%     ~~Describe \code{cov.function} here~~
The name, a text string of the covariance function.
}
  \item{m}{
%%     ~~Describe \code{m} here~~
The degree of the polynomial used in teh fixed part is (m-1)
}
  \item{chol.args}{
%%     ~~Describe \code{chol.args} here~~
A list of optional arguments (pivot, nnzR) that will be used 
  with the call to the cholesky decomposition. Pivoting is done 
  by default to make use of sparse matrices when they are 
  generated. This argument is useful in some cases for sparse 
  covariance functions to reset the memory parameter nnzR. 
  (See example below.)
}
  \item{find.trA}{
%%     ~~Describe \code{find.trA} here~~
If TRUE will estimate the effective degrees of freedom using 
  a simple Monte Carlo method. This will add to the computational 
  burden by approximately NtrA solutions of the linear system but 
  the cholesky decomposition is reused.
}
  \item{NtrA}{
%%     ~~Describe \code{NtrA} here~~
Number of Monte Carlo samples for the trace. But if NtrA is 
  greater than or equal to the number of observations the trace 
  is computed exactly.
}
  \item{iseed}{
%%     ~~Describe \code{iseed} here~~
Random seed ( using \code{set.seed(iseed)}) used to generate 
  iid normals for Monte Carlo estimate of the trace.
}
  \item{cov.args}{
%%     ~~Describe \code{cov.args} here~~
A list of optional arguments that will be used in calls to the 
  covariance function.
}
  \item{llambda}{
%%     ~~Describe \code{llambda} here~~
If not \code{NULL} then \code{lambda = exp( llambda)}
}
  \item{nGPUs}{
%%     ~~Describe \code{nGPUs} here~~
The number of GPUs to use with MAGMA
}
  \item{MAGMArank}{
%%     ~~Describe \code{MAGMArank} here~~
The GPU ID (rank) of the GPU to use with MAGMA.  This is 
  only used if \code{nGPUs = 1}.  If \code{nGPUs = n}, then 
  ranks 0 to n-1 are used.
}
  \item{singlePrecision}{
%%     ~~Describe \code{singlePrecision} here~~
If \code{TRUE}, Cholesky decomposition is performed in MAGMA
  using single precision.  Otherwise uses MAGMA's double 
  precision Cholesky decompostion.  Using single precision 
  will increase computation speed but reduce accuracy.
}
  \item{\dots}{
%%     ~~Describe \code{\dots} here~~
In mKrig and predict additional arguments that will be passed 
  to the covariance function.
}
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
An \code{mKrig} object.  See \code{?fields::mKrig} for details.
}
\references{
%% ~put references to the literature/web site here ~
http://cran.r-project.org/web/packages/fields/fields.pdf
http://www.image.ucar.edu/~nychka/Fields/
http://icl.cs.utk.edu/magma/
}
\author{
%%  ~~who you are~~
John Paige
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
\code{\link{fields}}
\code{\link[fields]{mKrig}}
\code{\link{stationary.fast.cov}}
\code{\link{mKrigMAGMA.MLE}}
}
\examples{
#use examples from ?fields:::mKrig
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ Kriging }
\keyword{ MAGMA }
\keyword{ spatial }% __ONLY ONE__ keyword per line
