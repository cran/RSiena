\name{includeInteraction}
\alias{includeInteraction}
\title{
Function to create user-specified interactions for a Siena model.
}
\description{This function provides an interface to allow easy update of
  an unspecified interaction row in a Siena effects object.
}
\usage{
includeInteraction(myeff, ..., include = TRUE,
name = myeff$name[1], type = "eval",
interaction1 = rep("", 3), interaction2 = rep("", 3),
character=FALSE)
}
\arguments{
  \item{myeff}{a Siena effects object as created by \code{\link{getEffects}}}
  \item{\dots}{
2 or 3 short names to identify the effects which should be interacted.
}
  \item{include}{Boolean. default TRUE, but can be switched to FALSE to
	turn off an interaction.
}
  \item{name}{Name of network for which interactions are being
	defined. Defaults to the first in the effects object.
}
\item{type}{
  Type of effects to be interacted.
}
\item{interaction1}{
  Vector of siena objects where needed to completely identify the
  effect e.g. covariate name or behavior variable name.
  Trailing blanks may be omitted.
}
  \item{interaction2}{
  Vector of siena objects where needed to completely identify the
  effect e.g. covariate name or behavior variable name.
  Trailing blanks may be omitted.
}
\item{character}{Boolean: are the effect names character strings or not}
}
\details{The details provided should uniquely identify up to three
  effects. If so, an interaction effect will be created and included or
  not in the model.
  }
  \value{
	An updated version of the input effects object, with the include
	column and the effect1 and effect2 and possibly effect3 columns of one row
	updated.
	Details of the fields altered will be printed.
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
myeff <- includeInteraction(myeff, transTrip, egoX,
interaction1=c("", "mynet2"))
}
\keyword{ classes }
