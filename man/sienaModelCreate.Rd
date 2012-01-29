\name{sienaModelCreate}
\alias{sienaModelCreate}
\alias{sienaModel}
\alias{model.create}
\title{Function to create an object containing the algorithm
  specifications
  for parameter estimation in RSiena}
\description{
  Creates an object with specifications for the algorithm for parameter
  estimation in RSiena.
}
\usage{
sienaModelCreate(fn, projname = "Siena", MaxDegree = 0,
     useStdInits = FALSE, n3 = 1000, nsub = 4, maxlike = FALSE,
     diag = !maxlike, condvarno = 0, condname = "", firstg = 0.2,
     cond = NA, findiff = FALSE, seed = NULL, pridg=0.05,
     prcdg=0.05, prper=0.2, pripr=0.3, prdpr=0.3, prirms=0.05,
     prdrms=0.05, maximumPermutationLength=40,
     minimumPermutationLength=2, initialPermutationLength=20,
     modelType=1, mult=5)
}
\arguments{
  \item{fn}{Function to do one simulation in the Robbins-Monro
    algorithm. Not to be touched.}
  \item{projname}{Character string name of project; the output file will be
    called projname.out. No embedded spaces!!!}
  \item{MaxDegree}{Named vector of maximum degree values for
	corresponding networks. Allows to restrict the model to networks
	with degrees not higher than this maximum.}
  \item{useStdInits}{Boolean. If TRUE, the initial values in the effects
	object will be ignored and default values used instead.
	If FALSE, the initial values in the effects object will be used.}
  \item{n3}{Number of iterations in phase 3.}
  \item{nsub}{Number of subphases in phase 2.}
  \item{maxlike}{Whether to use maximum likelihood method or
	Method of Moments estimation.}
  \item{diag}{Boolean: if FALSE,
	use the complete estimated derivative matrix in
  the Robbins-Monro procedure; if TRUE just use the diagonal entries.}
\item{condvarno}{If \code{cond} (conditional simulation), the
  sequential number of the network
  or behavior variable on which to condition.}
  \item{condname}{If conditional, the name of the dependent variable on
  which to condition. Use one or other of \code{condname} or
  \code{condvarno} to specify the variable.}
\item{firstg}{Initial value of scaling ('gain') parameter for updates
  in the Robbins-Monro procedure.}
  \item{cond}{Boolean. If TRUE, use conditional simulation. If missing,
	decision is deferred until \code{\link{siena07}}, when it is set to
	TRUE if there is only one dependent variable, FALSE otherwise.}
  \item{findiff}{ Boolean: If TRUE, estimate derivatives using finite
  differences. If FALSE use scores.}
  \item{seed}{Integer. Starting value of random seed. Not used if
	parallel testing.}
  \item{pridg}{Real number. Probability used in Metropolis-Hastings
	routine in ML estimation.}
  \item{prcdg}{Real number. Probability used in Metropolis-Hastings
	routine in ML estimation.}
  \item{prper}{Real number. Probability used in Metropolis-Hastings
	routine in ML estimation.}
  \item{pripr}{Real number. Probability used in Metropolis-Hastings
	routine in ML estimation.}
  \item{prdpr}{Real number. Probability used in Metropolis-Hastings
	routine in ML estimation.}
  \item{prirms}{Real number. Probability used in Metropolis-Hastings
	routine in ML estimation.}
  \item{prdrms}{Real number. Probability used in Metropolis-Hastings
	routine in ML estimation.}
  \item{maximumPermutationLength}{Maximum length of permutation in steps in ML
	estimation}
  \item{minimumPermutationLength}{Minimum length of permutation in steps in ML
	estimation}
  \item{initialPermutationLength}{Initial length of permutation in steps in ML
	estimation}
  \item{modelType}{Type of model to be fitted: 1=directed, 2:6 for
	symmetric networks: 2=forcing, 3=Initiative model,
	4=Pairwise forcing model, 5=Pairwise mutual model, 6=Pairwise joint model}
  \item{mult}{Multiplication factor for maximum likelihood. Number of
	steps per iteration is set to this multiple of the total distance
	between the observations at start and finish of the wave.}
}
\details{Model specification is done via this object for
  \code{\link{siena07}}.
  This function creates an object with the elements required to control the
  Robbins-Monro algorithm. Those not
  available as arguments can be changed manually where desired.
}
\value{Returns a model object of class "sienaModel" containing:
  \item{projname}{String value of name of project.}
  \item{useStdInits}{Boolean, see above.}
  \item{checktime}{Boolean, set to TRUE: report time in the phases or not.}
  \item{n3}{number of iterations in Phase 3}
  \item{firstg}{Initial value of the scaling ('gain') parameter in the
	Robbins-Monro algorithm.}
  \item{maxrat}{Value used to control the maximum size of the jumps.}
  \item{maxmaxrat}{Value used to control the maximum size of the jumps.}
  \item{maxlike}{Boolean: is FRAN using maximum likelihood?}
  \item{FRANname}{Name of simulation function FRAN. Is derived by
	\code{sienaModelCreate} from
                  \code{fn} and \code{maxlike}.}
  \item{cconditional}{Boolean: is FRAN using conditional estimation?}
  \item{condvarno}{Number of dependent variable on which to condition.}
  \item{condname}{Name of dependent variable on which to condition.}
  \item{FinDiff.method}{Boolean: are derivatives calculated using
    finite differences?}
  \item{nsub}{Number of subphases in phase 2.}
  \item{diag}{Boolean: use only the diagonal of the derivative matrix?}
  \item{modelType}{Type of model to be fitted: 1=directed, 2:6 for
	symmetric networks: 2=forcing, 3=Initiative model,
	4=Pairwise forcing model, 5=Pairwise mutual model, 6=Pairwise joint model}
  \item{MaxDegree}{Named vector of maximum degree values, or NULL.}
  \item{randomSeed}{Integer. Starting value of random seed. Not present
	unless given in call.}
  \item{pridg}{Real number. Probability used in Metropolis-Hastings
	routine in ML estimation.}
  \item{prcdg}{Real number. Probability used in Metropolis-Hastings
	routine in ML estimation.}
  \item{prper}{Real number. Probability used in Metropolis-Hastings
	routine in ML estimation.}
  \item{pripr}{Real number. Probability used in Metropolis-Hastings
	routine in ML estimation.}
  \item{prdpr}{Real number. Probability used in Metropolis-Hastings
	routine in ML estimation.}
  \item{prirms}{Real number. Probability used in Metropolis-Hastings
	routine in ML estimation.}
  \item{prdrms}{Real number. Probability used in Metropolis-Hastings
	routine in ML estimation.}
  \item{maximumPermutationLength}{Maximum length of permutation in steps in ML
	estimation}
  \item{minimumPermutationLength}{Minimum length of permutation in steps in ML
	estimation}
  \item{initialPermutationLength}{Initial length of permutation in steps in ML
	estimation}
  \item{mult}{Multiplication factor for maximum likelihood. Number of
	steps per iteration is set to this multiple of the total distance
	between the observations at start and finish of the wave.}
}
\references{See \url{http://www.stats.ox.ac.uk/~snijders/siena/}}
\author{Ruth Ripley }
\seealso{This is for use in \code{\link{siena07}}.}
\examples{

mymodel <- sienaModelCreate(projname="NetworkDyn")
StdModel <- sienaModelCreate(projname="NetworkDyn",useStdInits=TRUE)
CondModel <- sienaModelCreate(projname="NetworkDyn",condvarno=1,cond=TRUE)
Max10Model <- sienaModelCreate(projname="NetworkDyn",MaxDegree=c(mynet=10))
# where mynet is the name of the network object created by sienaNet().

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ classes }
