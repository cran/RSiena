\name{sienaTimeTest}
\alias{sienaTimeTest}
\title{Functions to assess and account for time heterogeneity of parameters}
\description{
 Takes a \code{sienaFit} object from a \code{\link{siena07}} estimation
 and tests the addition of dummy parameters at waves \code{m=2...(M-1)}
 through the score type test of Schweinberger (2007b).
 Tests for joint significance,
 parameter-wise significance, individual significance, and one step
 estimation of the unrestricted model parameters are returned in a list.

 After assessing time heterogeneity, effects objects can be modified via
 the \code{timeDummy} column. Simply type the periods which you would
 like to have time dummied separated by commas into this column. Using
 the \code{timeDummy} column within the effects object, you may specify
 which time dummy interacted parameters are to be estimated.

 If you wish to use this function with sienaFit objects which use the finite
 differences method of derivative estimation, or which use maximum
 likelihood estimation, you must request the derivatives to be returned
 by wave using the \code{byWave=TRUE} option to \code{\link{siena07}}.
 }

 \usage{
sienaTimeTest(sienaFit, effects=NULL, condition=FALSE)
}
\arguments{
  \item{sienaFit}{ A \code{sienaFit} object returned by \code{siena07}. }
  \item{effects}{Optional vector of effect numbers to test. Use the
	number on the print of the sienaFit object.}
  \item{condition}{ Whether to orthogonalize individual score tests
  against base effects and un-estimated dummy terms or just base effects.}
}
\details{
  This test follows the score type test of Schweinburger (2007b) as
 implemented by Lospinoso et. al. (2010) by using statistics already
 calculated at each wave to populate vectors of partitioned moment
 functions corresponding to a restricted model (the model that has been
 fit by \code{sienaFit} Object) and an unrestricted model (which contains
 dummies for \code{m=2...(M-1)}). A covariance matrix of these statistics
 is calculated, and a delta matrix is constructed through the score
 functions' outer products with these statistics. Through an
 orthogonalization and a Delta method of approximation, the
 variance-covariance structure of the dummy statistics is calculated,
 and appropriate statistical tests can be used.

 If it is determined that a time heterogeneity occurs for any number of
 time periods or effects (or any combination therein), the
 \code{timeDummy} column provides facilities within a \code{sienaFit} object to
 quickly adjust the model and re-estimate in the usual way.
 The \code{includeTimeDummy} function can be used to add the desired
 dummies entries to the effects object.
}
\value{
  \code{sienaTimeTest} Returns a list containing many items,
  including the following:
  \item{JointTest }{A \code{chi^2} test for joint significance of the dummies.}
  \item{ParameterTest }{ A \code{chi^2} test for joint significance
    across dummies for each separate effect.}
  \item{IndividualTest }{ A matrix displaying initial estimates, one
	step estimates, and a p-value for \code{H0:} the unrestricted
	parameters are equal to zero.}
}

\references{See \url{http://www.stats.ox.ac.uk/~snijders/siena/}
  for general information on RSiena.

Lospinoso, J.A., Schweinberger, M., Snijders, T.A.B, and Ripley,
R.M. "Assessing and Accounting for Time Heterogeneity in Stochastic
Actor Oriented Models". Advances in Data Analysis and Computation.
Special Issue on Social Networks. Submitted. Available from
\url{http://www.stats.ox.ac.uk/~lospinos/}.

M. Schweinberger and T.A.B. Snijders. Markov models for digraph panel
data: Monte carlo-based derivative estimation. Comput. Stat. Data Anal.,
51(9):4465-4483, 2007.  ISSN 0167-9473.

M. Schweinberger. Chapter 4: Statistical Modeling of Network Panel Data:
Goodness of Fit. PhD thesis, University of Groningen, 2007.
}
\author{Josh Lospinoso}
\seealso{\code{\link{siena07}}, \code{\link{plot.sienaTimeTest}},
  \code{\link{includeTimeDummy}}}
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

## A demonstration of the plotting facilities, on a larger dataset:
\dontrun{
mymodel <- sienaModelCreate(fn=simstats0c, nsub=4, n3=500)
mynet1 <- sienaNet(array(c(s501, s502, s503, s501, s503, s502), dim=c(50, 50, 6)))
mydata <- sienaDataCreate(mynet1)
myeff <- getEffects(mydata)
myeff <- includeEffects(myeff, transTrip, balance)
myeff <- includeTimeDummy(myeff, recip, timeDummy="2,3,5")
myeff <- includeTimeDummy(myeff, balance, timeDummy="4")
myeff <- includeTimeDummy(myeff, density, timeDummy="all")
ansp <- siena07(mymodel, data=mydata, effects=myeff, batch=FALSE)
ttp <- sienaTimeTest(ansp)

## Pairwise plots show
plot(ttp, pairwise=TRUE)

## Time test plots show
plot(ttp, effects=1:4, dims=c(2,2))
}
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
