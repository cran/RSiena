\name{coCovar}
\alias{coCovar}
\title{Function to create a constant covariate object}
\description{
  This function creates a constant covariate object from a vector.
  }
\usage{
coCovar(val, nodeSet='Actors')
}
\arguments{
  \item{val}{Vector of covariate values}
  \item{nodeSet}{Name of node set: character string.
  If the entire data set contains more than one node set,
  then the node sets must be specified in all data objects.}
}
\details{
  When part of a Siena data object,
  the covariate is associated with the node set nodeSet
  of the Siena data object.
}
\value{
  Returns the covariate as an object of class "coCovar", in which
  form it can be used as an argument to \code{SienaDataCreate}.
}
\references{See \url{http://www.stats.ox.ac.uk/~snijders/siena/}}
\author{Ruth Ripley}
\seealso{\code{\link{sienaDataCreate}}, \code{\link{varCovar}},
\code{\link{coDyadCovar}}, \code{\link{varDyadCovar}}}
\examples{

myconstCovar <- coCovar(s50a[,1])

}
\keyword{classes}
