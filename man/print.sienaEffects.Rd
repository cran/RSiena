\name{print.sienaEffects}
\alias{print.sienaEffects}
\alias{summary.sienaEffects}
\alias{print.summary.sienaEffects}
\title{
  Print methods for Siena effects objects
}
\description{
  Print the major columns of the effects object. Or all, with any non
  atomic columns listed separately.
}
\usage{
\method{print}{sienaEffects}(x, fileName = NULL, includeOnly=TRUE,
expandDummies = FALSE, ...)
\method{summary}{sienaEffects}(object, fileName = NULL,
includeOnly=TRUE, expandDummies = FALSE, ...)
\method{print}{summary.sienaEffects}(x, fileName = NULL, ...)
}
\arguments{
  \item{object}{An object of class \code{sienaEffects}}
  \item{x}{An object of class \code{sienaEffects} or
	\code{summary.sienaEffects} as appropriate.
  }
  \item{fileName}{Character string denoting file name if file output desired.
}
  \item{includeOnly}{Boolean. If TRUE, only effects with the include
	flag TRUE will be printed.}
  \item{expandDummies}{Interpret the \code{timeDummy} column and show
	any effects which would be added by \code{sienaTimeFix}}
  \item{\dots}{For extra arguments (none used at present)}
}
\value{
  The function \code{print.sienaEffects} prints details of the main
  columns of the selected rows of the effects object.

 The function \code{summary.sienaEffects} checks the rows for valid
 printing via \code{print.data.frame} and excludes any that will
 fail. The OK columns are printed first, followed by any others.

 Output from either can be directed to a file by using the argument
 \code{filename}.
}
\references{
  See \url{http://www.stats.ox.ac.uk/~snijders/siena/}
}
\author{Ruth Ripley
}
\seealso{\code{\link{sienaTimeTest}}}
\examples{
mynet1 <- sienaNet(array(c(s501, s502, s503), dim=c(50, 50, 3)))
mynet2 <- sienaNet(s50a, type='behavior')
mycovar <- coCovar(rnorm(50))
mydyadcovar <- coDyadCovar(matrix(as.numeric(rnorm(2500) > 2), nrow=50))
mydata <- sienaDataCreate(mynet1, mynet2, mycovar, mydyadcovar)
myeff <- getEffects(mydata)
myeff
summary(myeff)
}
\keyword{methods}
