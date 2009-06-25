\name{model.create}
\alias{model.create}
\alias{sienaModel}
\title{Function to create a model object for parameters to Siena07}
\description{
 Creates a basic object which can be used as an argument for Siena07
}
\usage{
model.create(fn = simstats0c, projname = "Siena",
    MaxDegree = 0, useStdInits = FALSE, n3 = 1000,
    nsub = 4, maxlike = FALSE, diag = TRUE, condvarno = 0,
    condname = "", firstg = 0.2, cond = FALSE, findiff = FALSE,
    seed = NULL)
}
\arguments{
  \item{fn}{function to do one simulation in the Robbins-Monro
    algorithm}
  \item{projname}{character string name of project. No embedded spaces}
  \item{MaxDegree}{Named vector of maximum degree values for
	corresponding networks}
  \item{useStdInits}{Boolean. If TRUE, the initial values in the effects
	object will be ignored and default values used instead.}
  \item{n3}{Number of iterations in phase 3.}
  \item{nsub}{Number of subphases in phase 2.}
  \item{maxlike}{whether to use maximum likelihood method or
	straightforward simulation}
  \item{diag}{Boolean: if FALSE,
	use the complete estimated derivative matrix in}
  the Robbins-Munro procedure; if TRUE just use the diagonal entries.
  \item{condvarno}{If conditional, the sequential number of the network
  or behavior variable on which to condition}
  \item{condname}{If conditional, the name of the dependent variable on
  which to condition. Use one or other of \code{condname} or
  \code{condvarno} to specify the variable.}
  \item{firstg}{initial value of scaling parameter in the Robbins-Monro
  procedure.}
  \item{cond}{Boolean. If TRUE, use conditional simulation}
  \item{findiff}{ Boolean: If TRUE, estimate derivatives using finite
  differences. If FALSE use scores.}
  \item{seed}{Integer. Starting value of random seed. Not used if
  parallel testing.}
}
\details{Model specification is done via this object for Siena7. This
  function creates an object with the correct elements. Those not
  available as arguments can be changed manually where desired.
}
\value{Returns a model object of class "sienaModel" containing:
  \item{projname}{String value of name of project}
  \item{useStdInits}{Boolean, see above}
  \item{checktime}{Boolean, set to TRUE, report time in the phases or not}
  \item{n3}{number of iterations in Phase 3}
  \item{firstg}{Initial value of the scale in the Robbins-Munro algorithm}
  \item{maxrat}{value used to control the maximum size of the jumps}
  \item{maxmaxrat}{value used to control the maximum size of the jumps}
  \item{FRAN}{Simulation fn}
  \item{maxlike}{is FRAN using maximum likelihood?}
  \item{cconditional}{is FRAN using conditional estimation?}
  \item{condvarno}{number of dependent variable on which to condition}
  \item{condname}{name of dependent variable on which to condition}
  \item{FinDiff.method}{Boolean, should we calculate derivatives using
    finite differences?}
  \item{nsub}{Number of subphases in phase 2}
  \item{diag}{Boolean - use just the diagonal of the derivative matrix?}
  \item{ModelType}{always 1}
  \item{MaxDegree}{named vector of maximum degree values, or NULL}
}
\references{See \url{http://www.stats.ox.ac.uk/~snijders/siena/}}
\author{Ruth Ripley }
\examples{

mymodel <- model.create(findiff=TRUE, fn = simstats0c)


}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ classes }
