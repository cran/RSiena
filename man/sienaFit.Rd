\name{sienaFit.methods}
\alias{sienaFit.methods}
\alias{print.summary.sienaFit}
\alias{sienaFit}
\alias{xtable.sienaFit}
\alias{print.sienaFit}
\alias{summary.sienaFit}
\title{Methods for processing sienaFit objects
}
\usage{
\method{print}{sienaFit}(x, tstat=TRUE, ...)

\method{summary}{sienaFit}(object,  ...)

\method{print}{summary.sienaFit}(x,  ...)

\method{xtable}{sienaFit}(x, caption = NULL, label = NULL, align = NULL,
digits = NULL, display = NULL, ...)
}
\arguments{
  \item{object}{An object of class \code{sienaFit}}
  \item{x}{An object of class \code{sienaFit}, or
	\code{summary.sienaFit} as appropriate}
  \item{tstat}{Boolean: add the t-statistics for convergence to the report}
  \item{caption}{
	See documentation for \code{\link{xtable}}}
  \item{label}{
  See documentation for \code{\link{xtable}}}
  \item{align}{
  See documentation for \code{\link{xtable}}}
  \item{digits}{
  See documentation for \code{\link{xtable}}}
  \item{display}{
  See documentation for \code{\link{xtable}}}
\item{\dots}{
  Add extra parameters for  \code{\link{print.xtable}} here. e.g.
  \code{type}, \code{file}}
}
\description{
  \code{print}, \code{summary}, and \code{xtable} methods for
   \code{sienaFit} objects.
}
\value{
 The function \code{print.sienaFit} prints a table containing estimated
 parameter values, standard errors and (optionally)
 t-statistics for convergence.

 The function \code{summary.sienaFit} prints a table containing
 estimated parameter values, standard errors and t-statistics for
 convergence together with the covariance matrix of the estimates, the
 derivative matrix of expected statistics \code{X} by parameters, and the
 covariance matrix of the expected statistics \code{X}.

 The function \code{xtable.sienaFit} creates an object of class
 \code{xtable.sienaFit} which inherits from class \code{xtable} and
 passes an extra arguments to the \code{print.xtable}.
}
\references{See \url{http://www.stats.ox.ac.uk/~snijders/siena/}
}
\author{Ruth Ripley}

\seealso{\code{\link{xtable}}, \code{\link{print.xtable}},
  \code{\link{siena07}}}
\examples{
mymodel <- sienaModelCreate(fn=simstats0c, nsub=2, n3=100)
mynet1 <- sienaNet(array(c(tmp3, tmp4), dim=c(32, 32, 2)))
mydata <- sienaDataCreate(mynet1)
myeff <- getEffects(mydata)
ans <- siena07(mymodel, data=mydata, effects=myeff, batch=TRUE)
ans
summary(ans)
xtable(ans, type='html', file='ans.html')
}
\keyword{method}
