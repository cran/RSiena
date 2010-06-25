\name{includeTimeDummy}
\alias{includeTimeDummy}
\title{Functions to include time dummy effects in a Siena model}
\description{
  This function provides an interface to set the TimeDummy column on
  selected rows of a Siena effects object.
 }

\usage{
includeTimeDummy(myeff, ..., timeDummy="all", name=myeff$name[1],
type="eval", interaction1="", interaction2="", include=TRUE,
character=FALSE)
}
\arguments{
  \item{myeff}{a Siena effects object as created by \code{\link{getEffects}}
}
  \item{\dots}{
short names to identify the effects which should be included or excluded.
}
\item{timeDummy}{character string. Either "all" or the periods for which
  to create dummies, space delimited.}
  \item{include}{Boolean. default TRUE, but can be switched to FALSE to
	turn off an effect.
  }
  \item{name}{Name of network for which effects are being
	included. Defaults to the first in the effects object.
}
  \item{type}{
  Type of effects to be included.
}
  \item{interaction1}{
  Name of siena object where needed to completely identify the
  effects e.g. covariate name or behavior variable name.
}
  \item{interaction2}{
  Name of siena object where needed to completely identify the
  effects e.g. covariate name or behavior variable name.
}
\item{character}{Boolean: are the effect names character strings or not}
}
\details{
The arguments should identify the effects completely. The
  include column will be set to the value requested (TRUE or FALSE)}
\value{
	An updated version of the input effects object, with the include
	columns for one or more rows
	updated.
	Details of the rows altered will be printed.
}

\references{See \url{http://www.stats.ox.ac.uk/~snijders/siena/}
  for general information on RSiena.

}
\author{Josh Lospinoso}
\seealso{\code{\link{sienaTimeTest}}, \code{\link{getEffects}},
  \code{\link{siena07}}}
\examples{
## Estimate a restricted model
mymodel <- sienaModelCreate(fn=simstats0c, nsub=2, n3=100)
mynet1 <- sienaNet(array(c(s501, s502, s503), dim=c(50, 50, 3)))
mydata <- sienaDataCreate(mynet1)
myeff <- getEffects(mydata)
myeff <- includeEffects(myeff, transTrip, balance)
ans <- siena07(mymodel, data=mydata, effects=myeff, batch=TRUE)

## Conduct the score type test to assess whether heterogeneity is present.
tt <- sienaTimeTest(ans)

## Suppose that we wish to include two time dummies.
## Add them in the following way:
myeff <- includeTimeDummy(myeff, recip, balance, timeDummy="2")
ans2 <- siena07(mymodel, data=mydata, effects=myeff, batch=TRUE)

## Re-assess the time heterogeneity
tt2 <- sienaTimeTest(ans2)

## And so on..

## A demonstration of RateX heterogeneity. Note that rate
## interactions are not implemented in general, just for
## Rate x cCovar.
\dontrun{
mymodel <- sienaModelCreate(fn=simstats0c, nsub=2, n3=100)
mynet1 <- sienaNet(array(c(s501, s502, s503), dim=c(50, 50, 3)))
myccov <- coCovar(s50a[,1])
mydata <- sienaDataCreate(mynet1, myccov)
myeff <- getEffects(mydata)
myeff <- includeEffects(myeff, transTrip, balance)
myeff <- includeTimeDummy(myeff, RateX, type="rate",
			interaction1="myccov")
ans <- siena07(mymodel, data=mydata, effects=myeff, batch=TRUE)
}
}
\keyword{models}
