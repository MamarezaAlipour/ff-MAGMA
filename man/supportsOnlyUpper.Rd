\name{supportsOnlyUpper}
\alias{supportsOnlyUpper}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Tests if covariance function supports evaluation of upper triangle of covariance matrix only
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
Tests if the given covariance function supports evaluation over the 
  upper triangle (and diagonal) of the covariance matrix only.  If 
  this function returns \code{TRUE}, setting the \code{onlyUpper} 
  argument to \code{TRUE} will speep up the code by avoiding 
  unnecessary covariance function evaluations.
}
\usage{
supportsOnlyUpper(cov.fun, cov.args=NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{cov.fun}{
The covariance function tested for support for evaluation of only 
  the upper triangle of the covariance matrix.
}
  \item{cov.args}{
  A list of arguments of cov.fun
}
}
\details{
Currently only \code{stationary.fast.cov} with the \code{Exponential} 
  covariance function supports evaluation of the covariance matrix 
  over the upper triangle (and diagonal) only.  Soon \code{stationary.fast.cov} 
  when used with the \code{Matern} covariance function will support 
  this as well.  \code{onlyUpper} can be set to \code{TRUE} even if 
  this function returns false, but a warning will be returned stating 
  that this argument is not yet supported for the given covariance function.
}

\value{
A logical indicating whether the given covariance function with the 
  given arguments supports evaluation of the upper triangle of the 
  covariance matrix only.
}
\author{
%%  ~~who you are~~
John Paige
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
\code{stationary.fast.cov}.  This has an \code{onlyUpper} option allowing 
  the user to evaluate the covariance matrix over the upper triangle and 
  diagonal only.
}
\examples{
################
#Test covariance function to see if it supports evaluation of 
#covariance matrix over upper triangle only
################

supportsOnlyUpper("Exponential")
supportsOnlyUpper(Exponential)
supportsOnlyUpper(Matern)
supportsOnlyUpper(stationary.fast.cov, cov.args=list(Covariance=Exponential))
supportsOnlyUpper(stationary.fast.cov, cov.args=list(Covariance="Matern"))
}

