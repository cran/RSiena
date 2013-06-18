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
  Type of effects to be included: "eval", "endow", "creation", or "rate".
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
\details{The effects indicated by the arguments \ldots, \code{type},
  and (if necessary) \code{interaction1} and \code{interaction2}
  are included or excluded
  from the model specified by the effects object.
  The names \code{interaction1} and \code{interaction2} do not refer
  to interactions between effects, but to dependence of effects on
  other variables in the data set.
  The arguments should identify the effects completely.

  A list of all effects available in a given effects object (e.g.,
  \code{myeff}), including their names of dependent variables,
  effect names, short names,
  and values of interaction1 and interaction2 (if any),
  is obtained by executing \code{\link{effectsDocumentation}(myeff)}.

  The function \code{}includeEffects operates by providing an interface
  to set the "include" column on
  selected rows of the effects object, to the value requested (TRUE or FALSE).
}
\value{
	An updated version of the input effects object, with the include
	columns for one or more rows updated.
	Details of the rows altered will be printed.
}
\references{See \url{http://www.stats.ox.ac.uk/~snijders/siena/}}
\author{Ruth Ripley }
\seealso{
  \code{\link{getEffects}}, \code{\link{setEffect}},
  \code{\link{includeInteraction}}, \code{\link{effectsDocumentation}}
}

\examples{
mynet1 <- sienaDependent(array(c(s501, s502, s503), dim=c(50, 50, 3)))
mybeh  <- sienaDependent(s50a, type="behavior")
mydata <- sienaDataCreate(mynet1, mybeh)
myeff <- getEffects(mydata)
myeff <- includeEffects(myeff, transTrip, balance)
myeff <- includeEffects(myeff, avAlt, name="mybeh", interaction1="mynet1")
myeff
}
\keyword{ classes }
