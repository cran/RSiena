\name{setEffect}
\alias{setEffect}
\title{
  Function to set various columns in an effects object in a Siena model
}
\description{
  This function provides an interface to change various columns of a
  selected row of a Siena effects object
}
\usage{
setEffect(myeff, shortName, parameter = 0, fix = FALSE,
test = FALSE, initialValue = 0, include = TRUE,
name = myeff$name[1], type = "eval", interaction1 = "",
interaction2 = "")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{myeff}{a Siena effects object as created by \code{\link{getEffects}}
}
  \item{shortName}{A short name (no quotes) to identify the effect which should be changed.
}
  \item{parameter}{Integer value required. Default 0.
}
  \item{fix}{Boolean required. Default FALSE.
}
  \item{test}{Boolean required. Default FALSE.
}
  \item{initialValue}{Integer value required. Default 0.
}
  \item{include}{Boolean. default TRUE, but can be switched to FALSE to
	turn off an effect.}
  \item{name}{Name of network for which effects are being
	included. Defaults to the first in the effects object.
}
  \item{type}{Character string indicating the
  type of the effect to be changed : currently "rate", "eval" or
  "endow". Default "eval".
}
  \item{interaction1}{
  Name of siena object where needed to completely identify the
  effect e.g. covariate name or behavior variable name.
}
  \item{interaction2}{
  Name of siena object where needed to completely identify the
  effect e.g. covariate name or behavior variable name.
}
}
\details{
The arguments should identify the effects completely. The parm column
will be set to the parameter value requested. The
  include column will be set to the include value requested (TRUE or FALSE)}
\value{
	An updated version of the input effects object, with one rows
	updated.
	Details of the row altered will be printed.
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
myeff <- setEffect(myeff, outInv, 3) ##parameter
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ classes }
