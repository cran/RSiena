\name{updateTheta}
\alias{updateTheta}
\title{
Function to update the initial values of theta.
}
\description{
  This function copies the final values of any matching selected effects
  from a \code{\link{sienaFit}} object to a Siena effects object.
}
\usage{
updateTheta(effects, prevAns)
}
\arguments{
  \item{effects}{Object of class \code{\link{sienaEffects}}
}
\item{prevAns}{Object of class \code{\link{sienaFit}}
  as returned by \code{\link{siena07}}.
}
}
\details{
  The initial values of any selected effects
  in the input effects object which match an effect estimated in
  \code{prevAns} will be updated.
  If the previous run was conditional, the estimated rate parameters for
  the dependent variable on which the run was conditioned are added to
  the final value of theta.
}
\value{
  The effects object with initial value column updated.
}
\note{Using this function explicitly before calling \code{\link{siena07}}
  rather than using it via the argument
  \code{prevAns} of \code{\link{siena07}} will not permit the use of
  the previous derivative matrix. This will be inefficient if
  the new model has the same selected effects as the previous one,
  and the initial values in the \code{prevAns} object are close
  to the final estimates.
}
\references{See \url{http://www.stats.ox.ac.uk/~snijders/siena/}}
\author{Ruth Ripley }


\seealso{\code{\link{siena07}}, \code{\link{getEffects}}
}
\examples{
mynet1 <- sienaDependent(array(c(tmp3, tmp4), dim=c(32, 32, 2)))
mydata <- sienaDataCreate(mynet1)
myeff <- getEffects(mydata)
myeff <- includeEffects(myeff, transTrip)
myalgorithm <- sienaAlgorithmCreate(nsub=2, n3=100)
ans <- siena07(myalgorithm, data=mydata, effects=myeff,  batch=TRUE)
myeff <- updateTheta(myeff, ans)
}
\keyword{ models }
