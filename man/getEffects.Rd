\name{getEffects}
\alias{getEffects}
\alias{sienaEffects}
\alias{sienaGroupEffects}
\title{Function to create a Siena effects object}
\description{
  Creates a basic list of effects for each of the dependent
  variables in the input siena object.
}
\usage{
getEffects(x, nintn = 10, behNintn=4, getDocumentation=FALSE)
}
\arguments{
  \item{x}{an object of class `siena' or `sienaGroup'}
  \item{nintn}{Number of lines for user defined network interactions. }
  \item{behNintn}{Number of lines for user defined behavior interactions. }
  \item{getDocumentation}{Flag to allow documentation of internal
  functions, not for use by users.}
}
\details{
Considers all the elements of the input siena data object and creates
lists of effects for use in siena model fits. Note that the class of the
return object may be lost if the data.frame is edited using \command{fix}.
}
\value{An object of class "sienaEffects" or "sienaGroupEffects": this is a data
  frame, each part of which relates to one dependent variable in the input
  object, with columns
  \item{name}{name of the dependent variable}
  \item{effectName}{name of the effect}
   \item{functionName}{name of the function}
   \item{shortName}{short name for the effect}
   \item{interaction1}{second variable in interaction, if any}
   \item{interaction2}{third variable in interaction, if any}
 \item{type}{"eval","endow", or "rate"}
  \item{basicRate}{boolean: whether a basic rate parameter}
  \item{include}{boolean: include or not}
  \item{randomEffects}{boolean: random or fixed effect}
  \item{fix}{boolean: fix value or not}
  \item{test}{boolean: test required or not}
  \item{timeDummy}{comma separated list of periods, or "all", or ',' for none --
	which time dummy interacted parameters should be included?}
  \item{initialValue}{starting value for estimation}
  \item{parm}{parameter values}
    \item{functionType}{"objective", "rate"}
  \item{period}{period for basic rate parameters}
  \item{rateType}{"Structural", "covariate"}
  \item{effectFn}{here NULL, but could be replaced by a function later}
  \item{statisticFn}{here NULL, but could be replaced by a function
	later}
  \item{untrimmedValue}{Used to store initial values which could be trimmed}
  \item{netType}{"oneMode", "behavior", "bipartite"}
  \item{groupName}{name of relevant data object}
  \item{group}{sequential number of relevant data object in total}
}
\references{See \url{http://www.stats.ox.ac.uk/~snijders/siena/}}
\author{Ruth Ripley}
\seealso{\code{\link{sienaDataCreate}}}
\examples{
mynet1 <- sienaNet(array(c(s501, s502, s503), dim=c(50, 50, 3)))
mynet2 <- sienaNet(s50a, type='behavior')
mycovar <- coCovar(rnorm(50))
mydyadcovar <- coDyadCovar(matrix(as.numeric(rnorm(2500) > 2), nrow=50))
mydata <- sienaDataCreate(mynet1, mynet2, mycovar, mydyadcovar)
myeff <- getEffects(mydata)
}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{classes}
