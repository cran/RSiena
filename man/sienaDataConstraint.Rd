\name{sienaDataConstraint}
\alias{sienaDataConstraint}
\title{
Function to change the values of the constraints between networks.
}
\description{This function allows the user to change the constraints of
  "higher", "disjoint" and "atLeastOne" for a specified pair of networks
  in a Siena data object.
}
\usage{
sienaDataConstraint(x, net1, net2,
type = c("higher", "disjoint", "atLeastOne"), value = FALSE)
}
\arguments{
  \item{x}{Siena Data Object Maybe a group object?}
  \item{net1}{name of first network}
  \item{net2}{name of second network}
  \item{type}{one of "higher", "disjoint", "atleastOne". Default is "higher".
}
  \item{value}{Boolean giving the value
}
}
\details{
  The value of the appropriate attribute is set to the value requested.
}
\value{
  Updated Siena data object.
}
\references{See \url{http://www.stats.ox.ac.uk/~snijders/siena/}}
\author{Ruth Ripley}
\seealso{\code{\link{sienaDataCreate}}}

\examples{
mynet1 <- sienaNet(array(c(s501,s502,s503),dim=c(50,50,3)))
mynet2 <- sienaNet(s50a,type='behavior')
mydata <- sienaDataCreate(mynet1, mynet2)
mydata <- sienaDataConstraint(mydata, mynet2, mynet1, "higher", FALSE)
}
\keyword{ classes }
