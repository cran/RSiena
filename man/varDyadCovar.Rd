\name{varDyadCovar}
\alias{varDyadCovar}
\title{Function to create a changing dyadic covariate object.}
\description{
  This function creates a changing dyadic covariate object from an array.
}
\usage{
varDyadCovar(val, centered=TRUE, nodeSets=c("Actors", "Actors"),
    sparse=is.list(val), type=c("oneMode", "bipartite"))
}
%- maybe also "usage" for other objects documented here.
\arguments{
  \item{val}{Array of covariate values, third dimension is the
    time. Alternatively, a list of sparse matrices of type "dgTMatrix".}
  \item{centered}{Boolean: if TRUE, then the overall mean value is subtracted.}
  \item{nodeSets}{Names (character string) of the associated node sets.
   If the entire data set contains more than one node set,
   then the node sets must be specified in all data objects.}
  \item{sparse}{Boolean: whether sparse matrices or not.}
  \item{type}{oneMode or bipartite: whether the matrix refers to a one-mode
  or a bipartite (two-mode) network.}
}
\details{
  When part of a Siena data object, the covariate is assumed to be
  associated with the node sets named NodeSets
  of the Siena data object. The names of the associated node sets will only
  be checked when the Siena data object is created.
}
\value{
  Returns the covariate as an object of class "varDyadCovar", in which
  form it can be used as an argument to \code{SienaDataCreate}.
}
\references{See \url{http://www.stats.ox.ac.uk/~snijders/siena/}}
\author{Ruth Ripley}
\seealso{\code{\link{sienaDataCreate}}, \code{\link{coDyadCovar}},
\code{\link{coCovar}}, \code{\link{varCovar}}}
\examples{

mydyadvar <- varDyadCovar(array(c(s501, s502), dim=c(50, 50, 2)))

}
% Add one or more standard keywords, see file "KEYWORDS" in the
% R documentation directory.
\keyword{classes}
