\name{includeEffects}
\alias{includeEffects}
\title{
  Function to include effects in a Siena model
}
\description{
  This function provides an interface to set the include column on
  selected rows of a Siena effects object
}
\usage{
includeEffects(myeff, ..., include = TRUE, name = myeff$name[1],
 type = "eval", interaction1 = "", interaction2 = "", character=FALSE)
}
\arguments{
  \item{myeff}{a Siena effects object as created by \code{\link{getEffects}}
}
  \item{\dots}{
short names to identify the effects which should be included or excluded.
}
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
\details{The arguments should identify the effects completely. The
  include column will be set to the value requested (TRUE or FALSE).
}
\value{
	An updated version of the input effects object, with the include
	columns for one or more rows
	updated.
	Details of the rows altered will be printed.
}
\references{See \url{http://www.stats.ox.ac.uk/~snijders/siena/}}
\author{Ruth Ripley }
\seealso{
  \code{\link{getEffects}}
}

\examples{
mynet1 <- sienaNet(array(c(s501, s502, s503), dim=c(50, 50, 3)))
mynet2 <- sienaNet(s50a, type="behavior")
mydata <- sienaDataCreate(mynet1, mynet2)
myeff <- getEffects(mydata)
myeff <- includeEffects(myeff, transTrip, balance)
}
\keyword{ classes }
