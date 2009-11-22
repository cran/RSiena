\name{sienaNet}
\alias{sienaNet}
\title{Function to create a Siena network object }
\description{
  Creates a Siena network object from a matrix or array or list of sparse matrix
  of triples.
}
\usage{
sienaNet(netarray, type = c("oneMode", "bipartite", "behavior"),
nodeSet = "Actors", sparse = is.list(netarray))
}
\arguments{
  \item{netarray}{matrix (type="behavior" only)
    or array of values or list of sparse matrices of type "dgTMatrix"}
  \item{type}{type of network, default "oneMode"}
  \item{nodeSet}{character string naming the appropriate node set. A
	vector containing 2 character strings for a bipartite network:
	"rows" first, then "columns".}
  \item{sparse}{logical: set to TRUE if the data is in sparse matrix
    format, FALSE otherwise }
}
\details{
  Adds attributes so that the array or list of matrices can be used in a
  Siena model fit.
}
\value{
An object of class "sienaNet".
}
\references{See \url{http://www.stats.ox.ac.uk/siena/}}
\author{Ruth Ripley}
\seealso{\code{\link{sienaDataCreate}}}
\examples{

mynet1 <- sienaNet(array(c(s501,s502,s503),dim=c(50,50,3)))

}
\keyword{classes}
