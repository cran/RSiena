\name{coDyadCovar}
\alias{coDyadCovar}
\title{Function to create a constant dyadic covariate object.}
\description{
  This function creates a constant dyadic covariate object from a matrix.
}
\usage{
coDyadCovar(val, nodeSets=c("Actors","Actors"))
}
\arguments{
  \item{val}{Matrix of covariate values. }
  \item{nodeSets}{The name of the node sets with which this
    covariate is associated.}
}
\details{
  When part of a Siena data object, the covariate is assumed to be
  associated with the node sets named in nodeSets
  of the Siena data object. The name of the associated node sets will only
  be checked when the Siena data object is created.
}
\value{
   Returns the covariate as an object of class "coDyadCovar", in which
  form it can be used as an argument to SienaData.create.
}
\references{See \url{http://www.stats.ox.ac.uk/~snijders/siena/}}
\author{Ruth Ripley}
\seealso{\code{\link{sienaDataCreate}}}
\examples{

mydyadvar <- coDyadCovar(s503)

}
\keyword{classes}