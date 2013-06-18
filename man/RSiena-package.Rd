\name{RSiena-package}
\alias{RSiena-package}
\alias{RSiena}
\docType{package}
\title{
Simulation Investigation for Empirical Network Analysis
}
\description{
Fits statistical models to longitudinal sets of networks, and to
longitudinal sets of networks and behavioral variables.
Not only one-mode networks but also two-mode networks
and multivariate networks are allowed.
The models are stochastic actor-oriented models.

Package \code{"RSienaTest"} is the development version, and
is distributed through R-Forge, see
\url{http://r-forge.r-project.org/R/?group_id=461}.
Package \code{"RSiena"} is the official release.
}
\details{

  Use \code{\link{siena07}} to fit a model.

  Data objects can be created from matrices and
  vectors using \code{\link{sienaDependent}}, \code{\link{coCovar}} etc.,
  and finally \code{\link{sienaDataCreate}}.
  Another possibility (but with less flexibility) is via the
  Gui displayed by \code{\link{siena01Gui}}, or via
  \code{\link{sienaDataCreateFromSession}}.

  Effects are selected using an \emph{effects} object, which can be created
  using \code{\link{getEffects}}.

  Control of the estimation algorithm requires a
  \code{sienaAlgorithm} object that
  defines the settings (parameters) of the algorithm,\cr
  and which can be created by \code{\link{sienaAlgorithmCreate}}
  (alias \code{\link{sienaModelCreate}}).

  More detailed help is available in the manual which you can display
  using \cr
  \code{RShowDoc("RSiena_Manual", package="RSiena")}
\tabular{ll}{
Package: \tab RSiena\cr
Type: \tab Package\cr
Version: \tab 1.1-232\cr
Date: \tab 2013-06-18\cr
Depends: \tab R (>= 2.15.0)\cr
Imports: \tab Matrix\cr
Suggests: \tab tcltk, network, codetools, lattice, MASS, parallel,
		  xtable, tools\cr
SystemRequirements: \tab GNU make, tcl/tk 8.5, Tktable\cr
License: \tab GPL-2 \cr
LazyData: \tab yes\cr
NeedsCompilation: yes\cr
BuildResaveData: \tab no\cr
}
}
\author{
  Ruth Ripley, Krists Boitmanis, Tom Snijders.
  Contributions by Josh Lospinoso, Charlotte Greenan, Christian Steglich,
  Johan Koskinen, Mark Ortmann, and Nynke Niezink.

  Maintainer: Tom A.B. Snijders <tom.snijders@nuffield.ox.ac.uk>
}
\references{
  \itemize{
   \item Schweinberger, Michael, and Snijders,
    Tom A.B. (2007). Markov models for digraph panel data:
    Monte Carlo-based derivative estimation.
    \emph{Computational Statistics and Data Analysis} 51, 4465-4483.

    \item Snijders, Tom A.B. (2001).
    The statistical evaluation of social network dynamics.
    \emph{Sociological Methodology}, 31, 361-395.
   \item Snijders, Tom A.B., Steglich, Christian E.G., and Schweinberger,
    Michael (2007). Modeling the co-evolution of networks and behavior.
    Pp. 41-71 in \emph{Longitudinal models in the behavioral
    and related sciences},
    edited by Kees van Montfort, Han Oud and Albert Satorra; Lawrence Erlbaum.
   \item Steglich, Christian E. G., Snijders, Tom A. B., and
    Pearson, Michael A. (2010). Dynamic networks and behavior:
    Separating selection from influence. Sociological Methodology, 40, 329-393.
   \item Further see the extensive manual accessible by the command\cr
    \code{RShowDoc("RSiena_Manual", package="RSiena")}\cr
	and the website \url{http://www.stats.ox.ac.uk/~snijders/siena/} .
           }
}
\seealso{\code{\link{siena07}}}

\examples{
mynet1 <- sienaDependent(array(c(tmp3, tmp4), dim=c(32, 32, 2)))
mydata <- sienaDataCreate(mynet1)
myeff <- getEffects(mydata)
myeff <- includeEffects(myeff, transTrip)
myeff
myalgorithm <- sienaAlgorithmCreate(nsub=3, n3=200)
ans <- siena07(myalgorithm, data=mydata, effects=myeff, batch=TRUE)
summary(ans)
}
\keyword{ package }
