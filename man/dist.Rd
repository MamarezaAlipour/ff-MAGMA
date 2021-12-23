\name{dist}
\alias{dist}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Create pariwise distance matrix for a set of spatial locations
}
\description{
This function behaves exactly like \code{stats::dist} when there is
  only one input argument, creating a compact matrix of pairwise 
  distances between the input spatial locations (only half of the 
  matrix is stored to save space due to symmetry).  In the case when 
  two inputs are supplied, this function behaves exactly like 
  \code{fields::rdist}, creating a matrix where its elements are the 
  pairwise distances between elements of the first vector and 
  elements of the second vector.
}
\usage{
dist(x1, x2 = NULL, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x1}{
Matrix of first set of locations where each row gives the coordinate 
  of a particular point.
}
  
  \item{x2}{
Matrix of second set of locations where each row gives the coordinate 
  of a particular point.  If this is not passed or given as NULL, \code{x1} 
  is used.
  }
  
  \item{...}{
Other arguments to pass to \code{stats::dist}.  These arguments will 
  not be compatible with \code{fields::rdist}.
  }
}
\details{
If the returned matrix is compact (i.e. if only one location vector is 
  supplied as input), it is stored as the columns of the lower-triangle 
  of the true distance matrix in a vector, but the \code{print} method 
  is overridden to display different sectons of the matrix depending 
  on extra arguments supplied.

Although \code{dist} and \code{rdist} each support extra arguments, 
  these arguments are not compatible.
}
\value{
If one non-NULL input is supplied then an object of class "\code{dist}" 
  is returned.  Otherwise a matrix is returned.
}
\author{
John Paige
}
\seealso{
\code{\link{rdist}}, \code{\link{dist}}, \code{\link{rdist.earth}}, \code{\link{fields.rdist.near}}
}
\examples{

#if one argument is supplied, the output matrix is stored in compact form:
distOut = dist(1:10)
distOut
as.vector(distOut)

#if two arguments are given, the entire output matrix is stored
distOut = dist(1:10, 1:10)
distOut

}
