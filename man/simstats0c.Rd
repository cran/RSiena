\name{simstats0c}
\alias{simstats0c}
\title{A version of FRAN}
\description{
 The function to be called as 'FRAN' by \code{\link{siena07}}
 in case of estimation by Method of Moments. Calls compiled C++.
}
\usage{
simstats0c(z, x, INIT = FALSE, TERM = FALSE, initC=FALSE,
           data=NULL, effects=NULL, fromFiniteDiff=FALSE,
           profileData=FALSE, prevAns=NULL, returnDeps=FALSE)
}
\arguments{
  \item{z}{Control object, passed in automatically in  \code{\link{siena07}}.}
  \item{x}{Model object, passed in automatically in  \code{\link{siena07}}.}
  \item{INIT}{If TRUE, do initial processing. May be required to set up z.}
  \item{TERM}{If TRUE, do end processing. }
  \item{initC}{If TRUE, call is to setup the data and model in C++. For
	use with multiple processes only.}
  \item{data}{A sienaData object as returned by \code{\link{sienaDataCreate}}.}
  \item{effects}{List of sienaEffects objects as
	returned by \code{\link{getEffects}}.}
  \item{fromFiniteDiff}{Boolean used during calculation of derivatives
	by finite differences. Not for user use.}
  \item{profileData}{Boolean to force dumping of the data for profiling
	with \command{sienaProfile.exe}.}
  \item{prevAns}{An object of class "sienaFit" as returned by
   \code{\link{siena07}}, from which scaling
	information (derivative matrix and standard deviation of the
	deviations) will be extracted along with the latest version of the
	parameters which will be used as the initial values, unless the
	model requests the use of standard initial values. If the previous
	model is the same as the current one, Phase 1 will be omitted. If
	not, any parameter estimates for effects which are included in
	the new model will be used as initial
	values, but phase1 will still be carried out. 
	If the results used as \code{prevAns} are a reasonable starting point,
	this will increase efficiency of the algorithm. }
  \item{returnDeps}{Boolean. Whether to return the simulated networks in
	Phase 3.}
}
\details{
  This should be used for the element FRAN of the model object in case of
  estimation by Method of Moments.
  The arguments with no defaults must be passed in on the call to
  \code{\link{siena07}}. 
}
\value{Depends on the call. If \code{INIT} or \code{initC} or
  \code{TERM}
  are true, returns \code{z}, the control object. Otherwise,
  returns a list containing:
  \item{fra}{Simulated statistics.}
  \item{sc}{Scores with which to calculate the derivative (not phase 2
    or if using finite differences or maximum likelihood).}
  \item{dff}{For maximum likelihood, 2nd deriv??}
  \item{ntim}{For conditional processing, time taken.}
  \item{feasible}{Currently set to TRUE.}
  \item{OK}{Could be set to FALSE if serious error has occurred.}
  \item{sims}{A list of simulation results, one for each iteration. Each
	list consists of a list for each data object, each of which consists of
	a list for each network, each of which consists of a list for each period,
	each component of which is an edgelist in matrix form (the columns
	are from, to, value) (or vector for behavior variables).
	Only if \code{returnDeps} is TRUE.}
}
\references{See \url{http://www.stats.ox.ac.uk/~snijders/siena/}}
\author{Ruth Ripley }
\seealso{\code{\link{siena07}} }
\examples{

mynet1 <- sienaNet(array(c(tmp3, tmp4), dim=c(32, 32, 2)))
mydata <- sienaDataCreate(mynet1)
myeff<- getEffects(mydata)
myeff[myeff$effectName=='transitive triplets'&
            myeff$type=='eval','include']<- TRUE
mymodel<- sienaModelCreate(fn = simstats0c, nsub=2, n3=100)
ans<- siena07(mymodel, data=mydata, effects=myeff,  batch=TRUE)
}
\keyword{models}
