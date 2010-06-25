\name{edit.sienaEffects}
\alias{edit.sienaEffects}
\title{
  Allow editing of a sienaEffects object if a gui is available.
}
\description{
  Interactive editor for an effects object. A wrapper to edit.data.frame.
}
\usage{
\method{edit}{sienaEffects}(name, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{name}{An object of class \code{sienaEffects}}
  \item{\dots}{For extra arguments (none used at present)
}}
\details{Will be invoked by \code{fix(name)} for an object of class
  \code{sienaeffects}.
}
\value{
  The updated object. There is no backup copy, and the edits cannot be undone.
}
\references{
  See \url{http://www.stats.ox.ac.uk/~snijders/siena/}
}
\author{Ruth Ripley
}
\examples{mynet1 <- sienaNet(array(c(s501, s502, s503), dim=c(50, 50, 3)))
mynet2 <- sienaNet(s50a, type='behavior')
mycovar <- coCovar(rnorm(50))
mydyadcovar <- coDyadCovar(matrix(as.numeric(rnorm(2500) > 2), nrow=50))
mydata <- sienaDataCreate(mynet1, mynet2, mycovar, mydyadcovar)
myeff <- getEffects(mydata)
\dontrun{
fix(myeff)
}
}
\keyword{ methods }
