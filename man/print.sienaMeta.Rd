\name{print.sienaMeta}
\alias{print.sienaMeta}
\alias{sienaMeta}
\alias{summary.sienaMeta}
\alias{plot.sienaMeta}
\alias{print.summary.sienaMeta}
\title{Methods for processing sienaMeta objects}
\usage{
\method{print}{sienaMeta}(x, file=FALSE, ...)

\method{summary}{sienaMeta}(object, file=FALSE, extra=TRUE, ...)

\method{print}{summary.sienaMeta}(x, file=FALSE, extra=TRUE, ...)

\method{plot}{sienaMeta}(x,  ..., layout=c(2,2))
}
\description{
  \code{print}, \code{summary}, and \code{plot} methods for
   \code{sienaMeta} objects.
}
\arguments{
  \item{object}{An object of class \code{sienaMeta}}
  \item{x}{An object of class \code{sienaMeta}, or
	\code{summary.sienaMeta} as appropriate}
  \item{file}{Boolean: if TRUE, sends output to file named
	\code{x$projname.out}. If FALSE, output is to the terminal.}
	\item{extra}{Boolean: if TRUE, prints more information}
  \item{layout}{the vector giving number of rows and columns
       in the arrangement of the several panels in a 
       rectangular array, possibly spanning multiple pages}
  \item{\dots}{For extra arguments (none used at present)
}
}
\value{
 The function \code{print.sienaMeta} prints details of the merged
 estimates of the meta-analysis, with test statistics.

 The function \code{summary.sienaMeta} prints details as for the
 \code{print} method, but also details of the \code{sienaFit} objects
 included.

 Output from either can be directed to a file by using the argument
 \code{file}. It will be appended to any existing file of the same
 name: \code{projname.out}  where \code{projname} is the value of the
 argument to \code{siena08}.

 The function \code{plot.sienaMeta} plots estimates against standard
 errors for each effect, with reference lines added at the two-sided 
 significance threshold 0.05. It returns an object of class \code{trellis},
 of the \code{\link{lattice}}.package. Effects for which a score test
 was requested are not plotted.
}
\references{
  T. A. B. Snijders and Chris Baerveldt.
  Multilevel network study of the effects
  of delinquent behavior on friendship
  evolution. Journal of Mathematical Sociology, 27: 123--151, 2003.

  See also \url{http://www.stats.ox.ac.uk/~snijders/siena/}
}
\author{Ruth Ripley, Tom Snijders
}

\examples{
\dontrun{
mymodel <- sienaModelCreate(fn=simstats0c, nsub=2, n3=100)
mynet1 <- sienaNet(array(c(s501, s502), dim=c(50, 50, 2)))
mynet2 <- sienaNet(array(c(s502, s503), dim=c(50, 50, 2)))
mydata1 <- sienaDataCreate(mynet1)
mydata2 <- sienaDataCreate(mynet2)
myeff1 <- getEffects(mydata1)
myeff2 <- getEffects(mydata2)
myeff1 <- setEffect(myeff1, transTrip)
myeff2 <- setEffect(myeff2, transTrip)
myeff1 <- setEffect(myeff1, cycle3, fix=TRUE, test=TRUE)
myeff2 <- setEffect(myeff2, cycle3, fix=TRUE, test=TRUE)
ans1 <- siena07(mymodel, data=mydata1, effects=myeff1, batch=TRUE)
ans2 <- siena07(mymodel, data=mydata2, effects=myeff2, batch=TRUE)
meta <- siena08(ans1, ans2)
meta
summary(meta)
plo <- plot(meta)
plo
plo <- plot(meta, layout = c(1,1))
plo[3]

}}
\keyword{ methods }
