\name{siena07}
\alias{siena07}
\title{Function to fit a model}
\description{
 Fits a model using method of moments, based on straightforward
 simulation, conditional or otherwise, or on an MCMC
 simulation. Estimation is done using Robbins-Monro algorithm. Note that
 the particular model to be used is passed in on the model object, and
 data for the model must be passed in using named arguments as the \code{...}.
(See examples)}
\usage{
siena07(x, batch=FALSE, verbose=FALSE, silent=FALSE, useCluster=FALSE,
nbrNodes=2, initC=FALSE, clusterString=rep("localhost", nbrNodes),
tt=NULL, parallelTesting=FALSE, ...)
}
\arguments{
  \item{x}{ A model object}
  \item{batch}{ Desired interface: 'batch' is a small amount of printout
	to the  console}
  \item{verbose}{Produces various output to the console if TRUE}
  \item{silent}{Produces no output to the console if TRUE, even if batch
  mode}
  \item{useCluster}{Boolean: whether to use a cluster of processes}
  \item{nbrNodes}{Number of processes to use if useCluster is TRUE}
  \item{initC}{Boolean: set to TRUE if the simulation will use C
	routines (currently always needed). Only for use if using
	multiple processors, to ensure all copies are initialised
	correctly.}
  \item{clusterString}{Definitions of clusters. Default set up to use
	the localmachine only.}
  \item{tt}{A \code{tcltk} toplevel window. Used if called from the
	model options screen.}
  \item{parallelTesting}{Boolean. If TRUE, sets up random numbers to
	parallel those in Siena 3.}
  \item{\dots}{Arguments for the simulation function}
}
\details{
 Runs a three-phase model estimation. Phase 1 does a few iterations to
 estimate the derivative matrix of the targets with respect to the
 parameter vector. Phase 2 does the estimation. Phase 3 runs a
 simulation to estimate standard errors and check the model. The
 simulation function is called once for each iteration in these phases
 and also once to initialise the model fitting and once to complete it.
 Unless in batch mode, displays a tcl/tk screen to allow interruption
 and to show progress.
}
\value{
  Returns an object of class "sienaFit", some parts of which are:
  \item{theta }{Fitted value of theta}
  \item{covtheta}{Estimated covariance matrix of theta}
  \item{dfra}{Matrix of estimated derivatives}
  \item{sf}{Matrix of deviations from target in phase 3}
  \item{sf2}{Array of statistics from simulations in phase 3}
  \item{targets}{Observed statistics}
  \item{targets2}{Observed statistics for each wave}
  \item{ssc}{Score function contributions for each wave for each
	simulation in phase 3}
  \item{sims}{List of simulated networks, currently just one entry for
  the first network. Each a list of edgelists, one for each period.}
}
\references{See \url{http://www.stats.ox.ac.uk/~snijders/siena/}}
\author{Ruth Ripley }
\seealso{There are print, summary and xtable methods for \code{\link{sienaFit}}
  objects
  \code{\link{sienaModelCreate}}, \code{\link{print.sienaFit}}}
\examples{
mymodel <- sienaModelCreate(fn=simstats0c, nsub=2, n3=100)
mynet1 <- sienaNet(array(c(tmp3, tmp4), dim=c(32, 32, 2)))
mydata <- sienaDataCreate(mynet1)
myeff <- getEffects(mydata)
ans <- siena07(mymodel, data=mydata, effects=myeff, batch=TRUE)

#or for conditional estimation
\dontrun{
mymodel$condname <- 'mynet1'
mymodel$cconditional <- TRUE
ans <- siena07(mymodel, data=mydata, effects=myeff, batch=TRUE)}
}
\keyword{models}
