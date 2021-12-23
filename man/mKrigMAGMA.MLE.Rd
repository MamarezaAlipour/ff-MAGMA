\name{mKrigMAGMA.MLE}
\alias{mKrigMAGMA.MLE}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
Accelerated mKrig.MLE
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
An accelerated version of \code{mKrig.MLE} from the \code{fields} package.  
  It is able to use the MAGMA library and thereby use GPUs to accelerate 
  the Cholesky decomposition required in \code{mKrig}.  It has also been 
  reorganized to precompute and store the distance matrix, reusing it 
  for each likelihood evaluation.
}
\usage{
mKrigMAGMA.MLE(x, y, weights = rep(1, nrow(x)), 
  Z = NULL, ..., par.grid = NULL, lambda = NULL, lambda.profile = TRUE, 
  verbose = FALSE, relative.tolerance = 1e-04,
  nGPUs = 0, MAGMArank = 0, singlePrecision=FALSE)
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
  \item{\dots}{
Additional arguments that would also be included in a call to
  \code{mKrig} to specify the covariance model and fixed model
  covariables.
}
  \item{par.grid}{
A list or data frame with components being parameters for
  different covariance models. A typical component is "theta"
  comprising a vector of scale parameters to try. If par.grid
  is "NULL" then the covariance model is fixed at values that
  are given in \dots.
}
  \item{lambda}{
%%     ~~Describe \code{lambda.guess} here~~
If \code{lambda.profile=FALSE} the values of lambda to evaluate
  the likelihood if "TRUE" the starting values for the
  optimization. If lambda is NA then the optimum value from
  previous search is used as the starting value. If lambda is
  NA and it is the first value the starting value defaults to
  1.0.
}
  \item{lambda.profile}{
  If \code{TRUE} maximize likelihood over lambda.
}
  \item{verbose}{
%%     ~~Describe \code{verbose} here~~
If \code{TRUE} print out interesting intermediate results.
}
  \item{relative.tolerance}{
Tolerance used to declare convergence when
  maximizing likelihood over lambda.
}
  \item{nGPUs}{
The number of GPUs to use with MAGMA
}
  \item{MAGMArank}{
The GPU ID (rank) of the GPU to use with MAGMA.  This is 
  only used if \code{nGPUs = 1}.  If \code{nGPUs = n}, then 
  ranks 0 to n-1 are used.
}
  \item{singlePrecision}{
If \code{TRUE}, Cholesky decomposition is performed in MAGMA
  using single precision.  Otherwise uses MAGMA double 
  precision Cholesky decompostion.  Using single precision 
  will increase computation speed but reduce accuracy.
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
See \code{?fields:::mKrig.MLE}.
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
The same as \code{mKrig.MLE} in \code{fields}.  See \code{?fields::mKrig.MLE}
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
\code{\link[fields]{mKrig.MLE}}
\code{\link{mKrigMAGMA}}
\code{\link{stationary.fast.cov}}
}
\examples{
#use examples from ?fields:::mKrig.MLE
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ Kriging }
\keyword{ MAGMA }% __ONLY ONE__ keyword per line
\keyword{ MLE }
\keyword{ spatial }
