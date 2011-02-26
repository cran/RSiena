\name{siena08}
\alias{siena08}
\title{
Function to perform a meta analysis of a collection of Siena fits.
}
\description{
Estimates a meta analysis based on  a collection of Siena fits.
}
\usage{
siena08(..., projname = "sienaMeta", bound = 5)
}
\arguments{
  \item{\dots}{names of Siena fit objects, returned from
	\code{\link{siena07}}. They will be renamed if entered in format
	\code{newname=oldname}.
}
\item{projname}{
  Base name of report file if required
}
\item{bound}{
  Upper limit of standard error for inclusion in the meta analysis.
}
}
\details{
  A meta analysis as described in the reference below is performed. This
  consists of two parts: an iterated weighted least squares and some
  calculations of Fisher combinations of p-values.
}
\value{An object of class \code{sienaMeta}.  There are \code{print},
  \code{summary} and \code{plot} methods for this class.

  An object of class \code{sienaMeta} is a list containing at least the
  following. (Items \code{cor.est} to
  \code{ns} appear once for each effect.)
  \item{cor.est}{Pearson correlation coefficient between estimates ans their standard errors.}
  \item{cor.pval}{p-value for above}
  \item{regfit}{Part of the result of the fit of \code{\link{iwlsm}}.}
  \item{regsummary}{The summary of the fit, which includes the
  coefficient table.}
  \item{Tsq}{test statistic for effect zero in every model}
  \item{pTsq}{p-value for above}
  \item{tratio}{test statistics that mean effect is 0}
  \item{ptratio}{p-value for above}
  \item{Qstat}{Test statistic for variance of effects is zero}
  \item{pttilde}{p-value for above}
  \item{cjplus}{Test statistic for at least one theta strictly greater than 0}
  \item{cjminus}{Test statistic for at least one theta strictly less than 0}
  \item{cjplusp}{p-value for \code{cjplus}}
  \item{cjminusp}{p-value for \code{cjminus}}
  \item{n1}{Number of fits on which the meta analysis is based}
  \item{scoreplus}{Test statistic for combination of right one-sided
  p-values from score tests}
  \item{scoreminus}{Test statistic for combination of left one-sided
  p-values from score tests}
  \item{scoreplusp}{pvalue for \code{scoreplus}}
  \item{scoreminusp}{p-value for \code{scoreminus}}
  \item{ns}{Number of fits on which the score test analysis is based}
\item{thetadf}{Data frame containing the coefficients, standard errors
  and score test results}
\item{projname}{Name for any output file to be produced by the print method}
\item{bound}{Estimates with standard error above this value were
  excluded from the calculations.}
\item{scores}{Object of class \code{by} indicating, for each effect in
  the models, whether score test  information was present. }
}
\references{
  T. A. B. Snijders and Chris Baerveldt.
  Multilevel network study of the effects
  of delinquent behavior on friendship
  evolution. Journal of Mathematical Sociology, 27: 123--151, 2003.
  See also \url{http://www.stats.ox.ac.uk/~snijders/siena/}
}
\author{Ruth Ripley}

\seealso{\code{\link{sienaMeta}}, \code{\link{siena07}}
}
\examples{
\dontrun{
mymodel <- sienaModelCreate(fn=simstats0c, nsub=2, n3=100)
mynet1 <- sienaNet(array(c(s501, s502), dim=c(50, 50, 2)))
mynet2 <- sienaNet(array(c(s502, s503), dim=c(50, 50, 2)))
mydata1 <- sienaDataCreate(mynet1)
mydata2 <- sienaDataCreate(mynet2)
myeff1 <- getEffects(mydata1)
myeff2 <- getEffects(mydata2)
myeff1 <- setEffect(myeff1, transTrip, fix=TRUE, test=TRUE)
myeff2 <- setEffect(myeff2, transTrip, fix=TRUE, test=TRUE)
myeff1 <- setEffect(myeff1, cycle3, fix=TRUE, test=TRUE)
myeff2 <- setEffect(myeff2, cycle3, fix=TRUE, test=TRUE)
ans1 <- siena07(mymodel, data=mydata1, effects=myeff1, batch=TRUE)
ans2 <- siena07(mymodel, data=mydata2, effects=myeff2, batch=TRUE)
meta <- siena08(ans1, ans2)
}
}
\keyword{models}
