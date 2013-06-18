\name{sienaDataCreate}
\alias{sienaDataCreate}
\alias{siena}
\title{Function to create a Siena data object}
\description{
Creates a Siena data object from input networks, covariates, and
composition change objects.}
\usage{
sienaDataCreate(..., nodeSets=NULL, getDocumentation=FALSE)
}
\arguments{
  \item{\dots}{objects of class \code{\link{sienaDependent}},
    \code{\link{coCovar}},\code{\link{varCovar}},\code{\link{coDyadCovar}},\cr
  \code{\link{varDyadCovar}}, \code{\link{sienaCompositionChange}}}
  \item{nodeSets}{list of Siena node sets. Default is the single node
    set named 'Actors', length equal to the number of rows in the first
    object of class "sienaDependent". If the entire data set contains more than
    one node set, then the node sets must have been specified in
    the creation of all data objects.}
  \item{getDocumentation}{Flag to allow documentation of internal
    functions, not for use by users.}
}
\details{
Checks that the objects fit, that there is at least one network,
and adds various attributes to each dependent variable describing the
data. If there is more than one nodeSet they must all be specified.
}
\value{
An object of class "siena" which is designed to be used in a siena
model fit. The components of the object are.
\item{nodeSets}{List of node sets involved}
\item{observations}{Integer indicating number of waves of data}
\item{depvars}{List of networks and behavior variables}
\item{cCovars}{List of constant covariates}
\item{vCovars}{List of changing covariates}
\item{dycCovars}{List of constant dyadic covariates}
\item{dyvCovars}{List of changing dyadic covariates}
\item{compositionChange}{List of composition change objects
  corresponding to the node sets}
}
\references{See \url{http://www.stats.ox.ac.uk/~snijders/siena/}}
\author{Ruth Ripley}
\seealso{
\code{\link{sienaDependent}}, \code{\link{coCovar}}, \code{\link{varCovar}},
\code{\link{coDyadCovar}}, \code{\link{varDyadCovar}},\cr
\code{\link{sienaCompositionChange}}, \code{\link{sienaGroupCreate}}}
\examples{
mynet <- sienaDependent(array(c(s501, s502, s503), dim=c(50, 50, 3)))
mybeh <- sienaDependent(s50a, type='behavior')
mydata <- sienaDataCreate(mynet, mybeh)}

\keyword{ classes }
