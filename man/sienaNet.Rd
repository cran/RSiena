\name{sienaNet}
\alias{sienaNet}
\title{Function to create a Siena network object }
\description{
  Creates a Siena network object from a matrix or array or list of sparse matrix
  of triples.
}
\usage{
sienaNet(netarray, type=c("oneMode", "bipartite", "behavior"),
nodeSet="Actors", sparse=is.list(netarray))
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
An object of class "sienaNet". An array or (networks only) a list of sparse
matrices with attributes:
  \item{netdims}{Dimensions of the network or behavior
	variable. Senders, receivers (1 for behavior), periods}
  \item{type}{oneMode, bipartite or behavior}
  \item{sparse}{Boolean: whether a list of sparse matrices or not}
  \item{nodeSet}{Character string with name(s) of node set(s)}
}
\references{See \url{http://www.stats.ox.ac.uk/siena/}}
\author{Ruth Ripley}
\seealso{\code{\link{sienaDataCreate}}}
\examples{

mynet1 <- sienaNet(array(c(s501, s502, s503), dim=c(50, 50, 3)))
mynet2 <- sienaNet(s50a, type="behavior")
## note that the following example works although the node sets do not exist!
mynet3 <- sienaNet(array(c(s501, s502, s503), dim=c(50, 50, 3)),
type="bipartite", nodeSet=c("senders", "receivers"))
## sparse matrix input - create some RSiena edgelists first
library(Matrix)
tmp1 <- as(Matrix(s501), "dgTMatrix")
tmp2 <- as(Matrix(s502), "dgTMatrix")
tmp3 <- as(Matrix(s503), "dgTMatrix")
mymat1 <- cbind(tmp1@i + 1, tmp1@j + 1, 1, 1)
mymat2 <- cbind(tmp2@i + 1, tmp2@j + 1, 1, 2)
mymat3 <- cbind(tmp3@i + 1, tmp3@j + 1, 1, 3)
mymat <- rbind(mymat1, mymat2, mymat3)
library(Matrix)
## mymat includes all 3 waves
mymatlist <- by(mymat, mymat[, 4], function(x)
    spMatrix(50, 50, x[, 1], x[, 2], x[, 3]))
mynet4 <- sienaNet(mymatlist)
## or alternatively
mymat1 <- mymat[mymat[, 4] == 1, ]
mymat2 <- mymat[mymat[, 4] == 2, ]
mymat3 <- mymat[mymat[, 4] == 3, ]
mymat1s <- spMatrix(50, 50, mymat1[, 1], mymat1[, 2], mymat1[, 3])
mymat2s <- spMatrix(50, 50, mymat2[, 1], mymat2[, 2], mymat2[, 3])
mymat3s <- spMatrix(50, 50, mymat3[, 1], mymat3[, 2], mymat3[, 3])
mynet4 <- sienaNet(list(mymat1s, mymat2s, mymat3s))
}
\keyword{classes}
