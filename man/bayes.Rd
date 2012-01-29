\name{bayes}
\alias{bayes}
\title{A function for fitting Bayesian models}
\description{A function to fit a Bayesian model to Siena Data
  objects. Uses the function \code{\link{maxlikec}} for the MCMC part,
  the Bayesian part is performed in R.
}
\usage{
bayes(data, effects, model, nwarm=100, nmain=100, nrunMHBatches=20,
      plotit=FALSE, nbrNodes=1, dfra=NULL, n=10,
      priorSigma=NULL, prevAns=NULL, clusterType=c("PSOCK", "FORK"),
      getDocumentation=FALSE)
}
\arguments{
  \item{data}{A sienaData object as returned by \code{\link{sienaDataCreate}}.}
  \item{effects}{List of sienaEffects objects as
	returned by \code{\link{getEffects}}.}
  \item{model}{Model object, as created by
	\code{\link{sienaModelCreate}}. Should contain all options required
	for the MCMC scheme, and a random seed if required.}
  \item{nwarm}{Number of iterations in the warm up phase.}
  \item{nmain}{Number of iterations in the main phase.}
  \item{nrunMHBatches}{Number of iterations of MCMC between Bayesian
	update steps.}
  \item{plotit}{Boolean: whether to plot parameters during the run}
  \item{nbrNodes}{Number of processes to be used. Cannot be more than
	the number of waves.}
  \item{dfra}{Hessian matrix. Optional}
  \item{n}{Number of iterations used to estimate Hessian}
  \item{priorSigma}{Covariance matrix to use for the prior distribution
	of the thetas.}
  \item{prevAns}{An object of class "sienaFit" as returned by
	\code{\link{siena07}}, from which scaling
	information (derivative matrix and standard deviation of the
	deviations) will be extracted along with the latest version of the
	parameters which will be used as the initial values, unless the
	model requests the use of standard initial values. If the previous
	model is exactly the same as the current one, no estimate of the
	Hessian will be made. If not, any parameter estimates for effects
	which are included in the new model will be used as initial
	values, but an estimate of the Hessian will still be made.
	If the results used as \code{prevAns} are a reasonable starting point,
	this will increase the efficiency of the algorithm. }
   \item{clusterType}{If using multiple processes, whether to use
	 forking processes or not. (Only "PSOCK" can be used on Windows.)}
   \item{getDocumentation}{Flag to allow documentation of internal
    functions, not for use by users.}
}
\details{
  This function wraps Bayesian sampling of parameters around calls to
  \code{\link{maxlikec}}. Unless a Hessian is supplied, one will be
  estimated from \code{n} batches of MH steps using the initial
  parameters. It is then transformed so the basic rate parameters are on
  a log scale. And inverted.
  It then attempts to scale the sampling
  covariance matrix to achieve about 40 out of 100 acceptances of Bayes
  proposals after single MH steps. Then a warming phase is done of
  \code{nwarm} Bayesian proposals each with 4 MH steps. Finally
  \code{nmain} repeats of (\code{nrunMHBatches} of \code{nrunMH} steps
  plus 1 Bayesian propsal) are performed. If \code{plotit} is TRUE,
  plots are produced at intervals during the run showing progress. They
  can be memory hogs, but the code may be useful for use on the returned
  values.
}
\value{
  Returns a list containing, among other things:
  \item{posteriorTot}{?}
  \item{posteriorMII}{?}
  \item{candidates}{array of parameters}
  \item{acceptances}{matrix of booleans: whether the corresponding
	change to the parameters was accepted. By group.}
  \item{MHacceptances}{array of acceptances of the MH steps, by step
	type and group but summed over dependent variables.}
  \item{MHrejections}{array of rejections of the ML steps}
  \item{MHproportions}{array of proportions of the MH steps accepted}
}
\references{See \url{http://www.stats.ox.ac.uk/~snijders/siena/}}
\author{Ruth Ripley }
\seealso{\code{\link{siena07}} }
\examples{
\donttest{
mynet1 <- sienaNet(array(c(tmp3, tmp4), dim=c(32, 32, 2)))
mydata <- sienaDataCreate(mynet1)
myeff <- getEffects(mydata)
myeff <- includeEffects(myeff, transTrip)
mymodel <- sienaModelCreate(maxlike=TRUE)
ans <- bayes(mydata, myeff, mymodel)
}
}
\keyword{models}
