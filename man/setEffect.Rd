\name{setEffect}
\alias{setEffect}
\title{
  Function to set various columns in an effects object in a Siena model.
}
\description{
  This function provides an interface to change various columns of a
  selected row of a Siena effects object.
}
\usage{
setEffect(myeff, shortName, parameter = NULL, fix = FALSE,
test = FALSE, random=FALSE, initialValue = 0, timeDummy = ",", include = TRUE,
name = myeff$name[1], type = "eval", interaction1 = "",
interaction2 = "", effect1=0, effect2=0, effect3=0,
period=1, group=1, character=FALSE)
}
%- maybe also "usage" for other objects documented here.
\arguments{
  \item{myeff}{a Siena effects object as created by \code{\link{getEffects}}
}
\item{shortName}{A short name (all with or all without quotes)
  to identify the effect which should be changed.
}
  \item{parameter}{Value of internal effect parameter. If NULL, no change is made.
}
  \item{fix}{For fixing effects. Boolean required. Default FALSE.
}
  \item{test}{For testing effects by score-type tests. Boolean required.
  Default FALSE.
}
  \item{random}{For specifying that effects will vary randomly;
  not relevant for \code{RSiena} at this moment.
  Boolean required. Default FALSE.
}
  \item{initialValue}{Initial value required. Default 0.
}
\item{timeDummy}{string: Comma delimited string of which periods to
dummy. Alternatively, use \code{\link{includeTimeDummy}}.
}
  \item{include}{Boolean. default TRUE, but can be switched to FALSE to
     turn off an effect.}
  \item{name}{Name of dependent variable (network or behavior)
    for which effects are being modified.
    Defaults to the first in the effects object.
}
  \item{type}{Character string indicating the
  type of the effect to be changed : currently "rate", "eval",
  "endow", or "creation". Default "eval".
}
  \item{interaction1}{
  Name of siena object where needed to completely identify the
  effect e.g. covariate name or behavior variable name.
}
  \item{interaction2}{
  Name of siena object where needed to completely identify the
  effect e.g. covariate name or behavior variable name.
}
\item{effect1}{
  Only for \code{shortName=unspInt}, which means this is a
  user-defined interaction effect: \code{effect1} is
  a natural number indicating the first component of the interaction effect;
  the number is the one listed when applying \code{print()} to \code{myeff}.
}
\item{effect2}{
  Only for \code{shortName=unspInt}: second component of interaction effect
  (see \code{effect1}).
}
\item{effect3}{
  Only for \code{shortName=unspInt}: third component of interaction effect,
  if any (see \code{effect1}).
}
\item{period}{
  Number of period if basic rate. Use numbering within groups.
}
\item{group}{
  Number of group if basic rate.
  }
\item{character}{Boolean: is the short name a character string or not.}
}
\details{
The arguments \code{shortName}, \code{name}, \code{type},
\code{interaction1}, \code{interaction2}, \code{effect1}, \code{effect2},
and \code{effect3} should identify the effects completely.
The column elements of the output effects object for \code{parm},
\code{fix}, \code{test}, \code{randomEffects},
\code{initialValue}, and \code{timeDummy}
will be set to the values requested.}
\value{
  An object of class \code{\link{sienaEffects}} or
  \code{\link{sienaGroupEffects}}. This will be an updated version of the
  input effects object, with one row updated.
  Details of the row altered will be printed.
}
\references{See \url{http://www.stats.ox.ac.uk/~snijders/siena/}}
\author{Ruth Ripley }
\seealso{
  \code{\link{getEffects}}, \code{\link{includeEffects}},
  \code{\link{includeInteraction}}.
}

\examples{
mynet <- sienaDependent(array(c(s501, s502, s503), dim=c(50, 50, 3)))
mybeh <- sienaDependent(s50a, type="behavior")
mydata <- sienaDataCreate(mynet, mybeh)
myeff <- getEffects(mydata)
# Set the initial rate parameter for one period:
myeff <- setEffect(myeff, Rate, initialValue=1.5, name="mybeh",
                   type="rate", period=2)
# Specify an effect parameter:
myeff <- setEffect(myeff, outTrunc, parameter=1)
myeff
# Now request the interaction between reciprocity and alcohol similarity
# to be tested without estimation.
myeff <- includeInteraction(myeff, recip, simX, interaction1=c("", "mybeh"))
myeff
# The following code uses the fact that the interaction effect is number 6.
# In practice one would just look up the numbers of effect1 and effect2
# from the response to \code{myeff} and fill them in.
(eff1 <- myeff[myeff$include,]$effect1[6])
(eff2 <- myeff[myeff$include,]$effect2[6])
myeff <- setEffect(myeff, unspInt, fix=TRUE, test=TRUE,
                                    effect1=eff1, effect2=eff2)
myeff
}
% Add one or more standard keywords, see file "KEYWORDS" in the
% R documentation directory.
\keyword{ classes }
