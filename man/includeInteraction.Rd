\name{includeInteraction}
\alias{includeInteraction}
\title{
Function to create user-specified interactions for a Siena model.
}
\description{This function allows the user to include or exclude
  an interaction effect in a Siena effects object.
}
\usage{
includeInteraction(myeff, ..., include = TRUE, name = myeff$name[1],
    type = "eval", interaction1 = rep("", 3), interaction2 = rep("", 3),
    character = FALSE, verbose = TRUE)
}
\arguments{
  \item{myeff}{a Siena effects object as created by \code{\link{getEffects}}}
  \item{\dots}{2 or 3 short names to identify the effects
               which should be interacted.}
  \item{include}{Boolean. default TRUE, but can be switched to FALSE to
	turn off an interaction.
}
  \item{name}{Name of dependent variable (network or behavior)
    for which interactions are being defined.
	Defaults to the first in the effects object.
}
\item{type}{
  Type of effects to be interacted.
}
\item{interaction1}{
  Vector of Siena objects where needed to completely identify the
  effect e.g. covariate name or behavior variable name.
  Trailing blanks may be omitted.
}
  \item{interaction2}{
  Vector of siena objects where needed to completely identify the
  effect e.g. covariate name or behavior variable name.
  Trailing blanks may be omitted.
}
\item{character}{Boolean: are the effect names character strings or not}
\item{verbose}{Boolean: should the print of altered effects be produced.}
}
\details{The details provided should uniquely identify up to three
  effects. If so, an interaction effect will be created and included or
  not in the model. \cr
  Whether interactions between two or three given effects
  can be created depends on their \code{interactionType} (which can be, for
  dependent network variables, empty, ego, or dyadic; and for dependent
  behavioral variables, empty or OK). Consult the section on Interaction
  Effects in the manual for this. The \code{interactionType} is shown
  in the list of effects obtained from the function
  \code{\link{effectsDocumentation}}.\cr
  From the point of view of model building it is usually advisable, when
  including an interaction effect in a model, also to include the
  corresponding main effects. This is however not enforced
  by \code{includeInteraction()}.

  The input names \code{interaction1} and \code{interaction2} do not themselves
  refer to a created interactions, but to dependence of the base effects on
  other variables in the data set.
  They are used to completely identify the effects.

  A list of all effects in a given effects object (e.g., \code{myeff}),
  including their names of dependent variables, effect names, short names,
  and values of interaction1 and interaction2 (if any),
  is obtained by executing \code{\link{effectsDocumentation}(myeff)}.

  }
\value{
  An updated version of the input effects object, with the "include"
  column and the "effect1" and "effect2" and possibly "effect3" columns
  of one row updated. \cr
  If \code{verbose=TRUE}, details of the fields altered will be printed.
}
\references{See \url{http://www.stats.ox.ac.uk/~snijders/siena/}}
\author{Ruth Ripley }
\seealso{
  \code{\link{getEffects}}, \code{\link{includeEffects}},
  \code{\link{effectsDocumentation}
}}
\examples{
mynet1 <- sienaDependent(array(c(s501, s502, s503), dim=c(50, 50, 3)))
mybeh  <- sienaDependent(s50a, type="behavior")
mydata <- sienaDataCreate(mynet1, mybeh)
myeff <- getEffects(mydata)
myeff <- includeEffects(myeff, transTrip)
myeff <- includeEffects(myeff, egoX, interaction1="mybeh")
myeff <- includeInteraction(myeff, transTrip, egoX,
         interaction1=c("", "mybeh"))
myeff
}
\keyword{ classes }
