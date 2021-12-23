\name{ff-MAGMA-package}
\alias{ff-MAGMA-package}
\alias{ff-MAGMA}
\docType{package}
\title{
\code{ff-MAGMA}: Accelerated \code{fields} Using MAGMA
}
\description{
This package has written new versions of several functions from \code{fields}, 
  most notably \code{mKrig} and \code{mKrig.MLE}.  The corresponding new 
  functions are \code{mKrigMAGMA} and \code{mKrigMAGMA.MLE}.  These output 
  results in the same format as \code{mKrig} and \code{mKrig.MLE}, but they have 
  been accelerated using the MAGMA library for heterogeneous architectures (those 
  with both CPUs and GPUs) and they have been reorganized so the distance 
  matrix does not need to be recomputed for each likelihood evaluation and 
  symmetry in distance and covariance matrices are exploited with the help
  of new covariance functions and utlity methods provided by this package.
}
\details{
\tabular{ll}{
Package: \tab ff-MAGMA\cr
Type: \tab Package\cr
Version: \tab 1.0\cr
Date: \tab 2015-05-21\cr
License: \tab Free\cr
}
See \code{fields}
}
\author{
John Paige

Maintainer: John Paige <paigejo@uw.edu>
}
\references{
http://cran.r-project.org/web/packages/fields/fields.pdf
http://www.image.ucar.edu/~nychka/Fields/
http://icl.cs.utk.edu/magma/
}

\keyword{ package }
\seealso{
\code{\link{fields}}
\code{\link{mKrigMAGMA.MLE}}
\code{\link{mKrigMAGMA}}
}
\examples{
#use examples from ?fields::mKrig and ?fields::mKrig.MLE
#but substitude mKrigMAGMA for mKrig and mKrigMAGMA.MLE 
#for mKrig.MLE
}
