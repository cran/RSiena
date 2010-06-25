\name{sienaNodeSet}
\alias{sienaNodeSet}
\title{Function to create a node set }
\description{Creates a Siena node set which can be used as the nodes in a
  siena network.
}
\usage{
sienaNodeSet(n, nodeSetName="Actors", names=NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{n}{integer, size of set.}
  \item{nodeSetName}{character string naming the node set.}
  \item{names}{optional character string vector of length n of the names of the
 nodes.}
}
\value{
  Returns a Siena node set, an integer vector, possibly with names, plus the
  attributes, class equal to 'sienaNodeSet', and nodeSetName equal to
  the argument nodeSetName.
}
\references{See \url{http://www.stats.ox.ac.uk/~snijders/siena/ }}
\author{Ruth Ripley }
\examples{
students <- sienaNodeSet(50, "student")
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{classes }
