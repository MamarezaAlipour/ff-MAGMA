\name{CovarianceUpper}
\alias{ExponentialUpper}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
Evaluate covariance over upper triangle of distance matrix
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
Evaluates the covariance over the upper triangle of a distance matrix 
  rather than over the entire matrix to reduce computation time.
}
\usage{
ExponentialUpper(distMat, range = 1, alpha = 1/range, phi = 1)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{distMat}{
The distance matrix to evaluate the covariance over.
}
  \item{range}{
Range parameter default is one. Note that the scale can also
  be specified through the "theta" scaling argument used in
  fields covariance functions)
}
  \item{alpha}{
1/range
}
  \item{phi}{
Marginal variance.
}
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
The covariance matrix, where only the upper triangle is calculated.
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
\code{\link[fields]{Exponential}}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ covariance }
