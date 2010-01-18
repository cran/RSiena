\name{RSiena-package}
\alias{RSiena-package}
\alias{RSiena}
\docType{package}
\title{
Simulation Investigation for Empirical Network Analysis
}
\description{
Fits models to longitudinal sets of networks.
}
\details{

  Use \code{\link{siena07}} to fit a model.

  Data objects can be created via the
  Gui displayed by \code{\link{siena01Gui}}, via \cr
  \code{\link{sienaDataCreateFromSession}} or directly from matrices and
  vectors using \code{\link{sienaNet}}, \code{\link{coCovar}} etc.\ and finally
  \code{\link{sienaDataCreate}}.

  Effects are selected using an \emph{effects} object, which can be created
  using \code{\link{getEffects}}.

  Control of the fitting requires a \emph{model} object, which can be
  created by \code{\link{sienaModelCreate}}.

  More detailed help is available in the manual which you can display
  using \cr
  \code{RShowDoc("s_man400", package="RSiena")}
\tabular{ll}{
Package: \tab RSiena\cr
Type: \tab Package\cr
Version: \tab 1.0.9\cr
Date: \tab 2010-01-18\cr
License: \tab GPL-2 \cr
LazyLoad: \tab yes\cr
}
}
\author{
Who wrote it

Maintainer: Who to complain to <ruth@stats.ox.ac.uk>
}
\seealso{\code{\link{siena07}}}
\references{See \url{http://www.stats.ox.ac.uk/~snijders/siena/}}
\keyword{ package }
\examples{
mynet1 <- sienaNet(array(c(tmp3, tmp4), dim=c(32, 32, 2)))
mydata <- sienaDataCreate(mynet1)
myeff <- getEffects(mydata)
mymodel <- sienaModelCreate(findiff=FALSE, fn=simstats0c, nsub=2, n3=100)
ans <- siena07(mymodel, data=mydata, effects=myeff, batch=TRUE)
}
