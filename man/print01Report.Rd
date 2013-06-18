\name{print01Report}
\Rdversion{1.1}
\alias{print01Report}
\title{
Function to produce the Siena01 report from R objects}
\description{
Prints a report of a Siena data object and its default effects.
}
\usage{
print01Report(data, myeff, modelname = "Siena", session = NULL,
getDocumentation=FALSE)
}
\arguments{
  \item{data}{a Siena data object
}
\item{myeff}{
  a Siena Effects object
}
\item{modelname}{
  Character string used to name the output file "modelname.out"
}
\item{session}{
  Used to pass in a Siena01Gui() style session object so that data file names
  can be printed.
}
  \item{getDocumentation}{Flag to allow documentation of internal
  functions, not for use by users.}
}
\details{
  First deletes any file of the name "modelname.out", then prints a new one.
  }
\value{
  No value returned.
}
\references{See \url{http://www.stats.ox.ac.uk/~snijders/siena/}}
\author{Ruth Ripley }
\seealso{
 \code{\link{siena01Gui}}
}

\examples{
mynet1 <- sienaDependent(array(c(s501, s502, s503), dim=c(50, 50, 3)))
mydata <- sienaDataCreate(mynet1)
myeff <- getEffects(mydata)
\dontrun{print01Report(mydata, myeff)}

}
\keyword{print}
