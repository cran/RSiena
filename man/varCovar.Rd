\name{varCovar}
\alias{varCovar}
\title{Function to create a changing covariate object.}
\description{
  This function creates a changing covariate object from a matrix.
}
\usage{
varCovar(val, centered=TRUE, nodeSet="Actors", imputationValues=NULL)
}
\arguments{
  \item{val}{ Matrix of covariate values, one row for each actor, one
    column for each period. }
  \item{centered}{Boolean: if TRUE, then the overall mean value is subtracted.}
  \item{nodeSet}{Character string containing the name of the associated
    node set. If the entire data set contains more than one node set,
    then the node sets must be specified in all data objects.}
  \item{imputationValues}{Mtrix of covariate values of same
   dimensions as \code{val}, to be used for imputation of
   \code{NA} values (if any) in \code{val}. Must not contain
   any \code{NA}.}
}
\details{
  When part of a Siena data object, the covariate is assumed to be
  associated with node set nodeSet of the Siena data object.\cr
  If there are any \code{NA} values in \code{val}, and
  \code{imputationValues} is given, then the corresponding elements
  of \code{imputationValues} are used for imputation.
  If \code{imputationValues} is \code{NULL}, imputation
  is by the overall mean value. In both cases, cases with imputed values
  are not used for calculating target statistics (see the manual).
}
\value{
  Returns the covariate as an object of class "varCovar", in which
  form it can be used as an argument to \code{\link{sienaDataCreate}}.
}
\references{See \url{http://www.stats.ox.ac.uk/~snijders/siena/}}
\author{Ruth Ripley}
\seealso{\code{\link{sienaDataCreate}}, \code{\link{coCovar}},
\code{\link{coDyadCovar}}, \code{\link{varDyadCovar}}}
\examples{

myvarCovar <- varCovar(s50a)

}
\keyword{classes}
