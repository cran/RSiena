\name{simstats0c}
\alias{simstats0c}
\title{A version of FRAN}
\description{
 The function to be called as 'FRAN'. Calls compiled C++.
}
\usage{
simstats0c(z, x, INIT = FALSE, TERM = FALSE, initC=FALSE,
    data=NULL, effects=NULL, fromFiniteDiff=FALSE,
    profileData=FALSE, prevAns=NULL)
}
\arguments{
  \item{z}{control object, passed in automatically in Siena07}
  \item{x}{model object, passed in automatically in Siena07}
  \item{INIT}{if TRUE, do initial processing. May be required to set up z}
  \item{TERM}{if TRUE, do end processing. }
  \item{initC}{ if TRUE, call is to setup the data and model in C++. For
	use with multiple processes only}
  \item{data}{A siena object }
  \item{effects}{list of data frames as returned by getEffects}
  \item{fromFiniteDiff}{Boolean used during calculation of derivatives
	by finite differences. Not for user use.}
  \item{profileData}{boolean to force dumping of the data for profiling
	with \command{sienaProfile.exe}}
  \item{prevAns}{An object of class "sienaFit", from which scaling
	information (derivative matrix and standard deviation of the
	deviations) will be extracted. Phase 1 will then be omitted.}
}
\details{
  This should be used for the element FRAN of the model object.
  The arguments with no defaults must be passed in on the call to siena07.
}
\value{Returns a list containing:
  \item{fra}{Simulated statistics}
  \item{sc}{Scores with which to calculate the derivative (not phase 2
    or if using finite differences or maximum likelihood)}
  \item{dff}{For maximum likelihood, 2nd deriv??}
  \item{ntim}{For conditional processing, time taken}
  \item{feasible}{currently set to TRUE}
  \item{OK}{could be set to FALSE if serious error has occurred}
}
\references{See \url{http://www.stats.ox.ac.uk/~snijders/siena/}}
\author{Ruth Ripley }
\seealso{\code{\link{siena07}} }
\examples{

mynet1 <- sienaNet(array(c(tmp3,tmp4),dim=c(32,32,2)))
mydata <- sienaDataCreate(mynet1)
myeff<- getEffects(mydata)
myeff[myeff$effectName=='transitive triplets'&
            myeff$type=='eval','include']<- TRUE
mymodel<- model.create(findiff=TRUE, fn = simstats0c)
ans<- siena07(mymodel, data=mydata, effects=myeff,  batch=TRUE)
}
\keyword{models}
