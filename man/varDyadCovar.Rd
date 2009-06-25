\name{varDyadCovar}
\alias{varDyadCovar}
\title{Function to create a changing dyadic covariate object.}
\description{
  This function creates a changing dyadic covariate object from an array.
}
\usage{
varDyadCovar(val, nodeSets=c('Actors','Actors'))
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{val}{Array of covariate values, Third dimension is the time }
  \item{nodeSets}{Names (character string) of the associated node sets}
}
\details{
  When part of a Siena data object, the covariate is assumed to be
  associated with the node sets named NodeSets
  of the Siena data object.  The names of the associated node sets will only
  be checked when the Siena data object is created.
}
\value{
  Returns the covariate as an object of class 'varDyadCovar', in which
  form it can be used as an argument to SienaData.create.
}
\references{See \url{http://www.stats.ox.ac.uk/~snijders/siena/}}
\author{Ruth Ripley}
\seealso{\code{\link{sienaDataCreate}}}
\examples{

mydyadvar <- varDyadCovar(array(c(s501,s502),dim=c(50,50,2)))

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{classes}
