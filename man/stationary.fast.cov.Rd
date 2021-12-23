\name{stationary.fast.cov}
\alias{stationary.fast.cov}
\alias{Exp.fast.cov}
\alias{matern.fast.cov}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
Fast Covariance Functions
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
These covariance functions are similar to covariance functions in \code{fields} 
  except they have additional arguments for precomputed distance 
  matrices and for calculating only the upper triangle and diagonal 
  of the output covariance matrix to save time.  Also, the \code{dist} 
  function for computing pairwise distances between observations is supported, 
  which only calculates the upper triangle of the distance matrix to save time.
}
\usage{
stationary.fast.cov(x1, x2=NULL, Covariance = "Exponential", Distance = "dist", 
  Dist.args = NULL, theta = 1, V = NULL, C = NA, marginal = FALSE, 
  derivative = 0, distMat = NA, onlyUpper = FALSE, ...)

Exp.fast.cov(x1, x2, theta = rep(1, ncol(x1)), distMat = NA, 
  C = NA, marginal = FALSE, onlyUpper=FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x1}{
%%     ~~Describe \code{x1} here~~
Matrix of first set of locations where each row gives the 
  coordinates of a particular point.
}
  \item{x2}{
%%     ~~Describe \code{x2} here~~
Matrix of second set of locations where each row gives the 
  coordinatesof a particular point. If this is missing x1 
  is used.
}
  \item{Covariance}{
%%     ~~Describe \code{Covariance} here~~
Character string that is the name of the covariance shape 
  function for the distance between locations. Choices in 
  fields are \code{Exponential}, \code{Matern}.
}
  \item{Distance}{
%%     ~~Describe \code{Distance} here~~
Character string that is the name of the distance function 
to use. Choices in \code{fields} are \code{rdist}, \code{rdist.earth}
}
  \item{Dist.args}{
%%     ~~Describe \code{Dist.args} here~~
A list of optional arguments to pass to the Distance function.
}
  \item{theta}{
%%     ~~Describe \code{theta} here~~
Range (or scale) parameter. This should be a scalar (use 
  the \code{V} argument for other scaling options). Any distance 
  calculated for a covariance function is divided theta 
  before the covariance function is evaluated.
}
  \item{V}{
%%     ~~Describe \code{V} here~~
A matrix that describes the inverse linear transformation 
  of the coordinates before distances are found. In R code 
  this transformation is: \code{x1 \%*\% t(solve(V))} Default 
  is \code{NULL}, that is the transformation is just dividing distance by 
  the scalar value theta. See Details below. If one has a 
  vector of "theta's" that are the scaling for each 
  coordinate then just express this as \code{V = diag(theta)} in 
  the call to this function.
}
  \item{C}{
%%     ~~Describe \code{C} here~~
A vector with the same length as the number of rows of \code{x2}. 
  If specified the covariance matrix will be multiplied by 
  this vector.
}
  \item{marginal}{
%%     ~~Describe \code{marginal} here~~
If \code{TRUE} returns just the diagonal elements of the covariance 
  matrix using the \code{x1} locations. In this case this is just 1.0. 
  The marginal argument will trivial for this function is a 
  required argument and capability for all covariance functions 
  used with \code{Krig}.
}
  \item{derivative}{
%%     ~~Describe \code{derivative} here~~
If nonzero evaluates the partials of the covariance function at locations 
  \code{x1}. This must be used with the ``\code{C}'' option and is mainly 
  called from within a predict function. The partial derivative is taken 
  with respect to \code{x1}.
}
  \item{distMat}{
%%     ~~Describe \code{distMat} here~~
If the distance matrix between \code{x1} and \code{x2} has already been 
  computed, it can be passed via this argument so it won't need to be 
  recomputed.
}
  \item{onlyUpper}{
%%     ~~Describe \code{onlyUpper} here~~
For internal use only, not meant to be set by the user.  Automatically 
  set to \code{TRUE} by \code{mKrigMAGMA.MLE} if \code{lambda.profile} 
  is set to \code{TRUE}, but set to \code{FALSE} for the final 
  parameter fit so output is compatible with rest of \code{fields}.
  
If \code{TRUE}, only the upper triangle and diagonal of the covariance 
  matrix is computed to save time (although if a non-compact distance 
  matrix is used, such as when not using \code{dist} with only one 
  argument, the onlyUpper argument is set to FALSE).  If \code{FALSE}, 
  the entire covariance matrix is computed.  In general, it should 
  only be set to \code{TRUE} for \code{mKrigMAGMA.MLE}, and the default 
  is set to \code{FALSE} so it is compatible with all \code{fields}.
}
  \item{\dots}{
%%     ~~Describe \code{\dots} here~~
Any other arguments that will be passed to the covariance function. e.g.
  \code{smoothness} for the Matern.
}
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
A covariance matrix, or, depending on the covariance function, 
  a matrix whose upper triangle and diagonal are those of a 
  covariance matrix.  The diagonal and upper triangle are the 
  same as those of the output of \code{fields::stationary.cov}.
}
\references{
%% ~put references to the literature/web site here ~
http://cran.r-project.org/web/packages/fields/fields.pdf
http://www.image.ucar.edu/~nychka/Fields/
}
\author{
%%  ~~who you are~~
John Paige
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
\code{\link[fields]{stationary.cov}}
\code{\link{dist}}
\code{\link[fields]{rdist}}
}
\examples{
################
#for data uniformly distributed on [0,1]x[0,1] square
################

#generate 5 random points on [0,1]x[0,1] square
x = matrix(runif(10), nrow=5)

#get compact form distance matrix
distMat = dist(x)

#generate covariance matrices
covMat = stationary.fast.cov(x, distMat=distMat)
covMatUpper = stationary.fast.cov(x, distMat=distMat, onlyUpper=TRUE)

#print results
print(distMat)
print(covMatUpper)
print(covMat)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ Kriging }
\keyword{ covariance }% __ONLY ONE__ keyword per line
