\name{sienaAlgorithmCreate}
\alias{sienaModelCreate}
\alias{sienaAlgorithmCreate}
\alias{sienaModel}
\alias{sienaAlgorithm}
\alias{model.create}
\title{Function to create an object containing the algorithm
  specifications for parameter estimation in RSiena}
\description{
  Creates an object with specifications for the algorithm for parameter
  estimation in RSiena.\cr
  \code{sienaAlgorithmCreate()} and \code{sienaModelCreate()}
  are identical functions; the second name was
  used from the start of the \code{RSiena} package, but the first name
  indicates more precisely the purpose of this function.
}
\usage{
sienaAlgorithmCreate(fn, projname = "Siena", MaxDegree = 0,
     useStdInits = FALSE, n3 = 1000, nsub = 4,
     dolby=TRUE, maxlike = FALSE, diagonalize=1.0*!maxlike,
     condvarno = 0, condname = "", firstg = 0.2,
     cond = NA, findiff = FALSE, seed = NULL, pridg=0.05,
     prcdg=0.05, prper=0.2, pripr=0.3, prdpr=0.3, prirms=0.05,
     prdrms=0.05, maximumPermutationLength=40,
     minimumPermutationLength=2, initialPermutationLength=20,
     modelType=1, mult=5, simOnly=FALSE)

sienaModelCreate(fn, projname = "Siena", MaxDegree = 0,
     useStdInits = FALSE, n3 = 1000, nsub = 4,
     dolby=TRUE, maxlike = FALSE, diagonalize=1.0*!maxlike,
     condvarno = 0, condname = "", firstg = 0.2,
     cond = NA, findiff = FALSE, seed = NULL, pridg=0.05,
     prcdg=0.05, prper=0.2, pripr=0.3, prdpr=0.3, prirms=0.05,
     prdrms=0.05, maximumPermutationLength=40,
     minimumPermutationLength=2, initialPermutationLength=20,
     modelType=1, mult=5, simOnly=FALSE)}
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
  \item{dolby}{Boolean. Should there be noise reduction by regression
    on augmented data score. In most cases dolby=TRUE yields better
	convergence; if convergence is problematic, however, dolby=FALSE
	may be tried. Just use whatever works best.}
  \item{maxlike}{Whether to use maximum likelihood method or
    Method of Moments estimation.}
  \item{diagonalize}{Number between 0 and 1 (bounds included),
    values outside this interval will be truncated;
    for diagonalize=0 the complete estimated derivative matrix will be used
    for updates in the Robbins-Monro procedure;
    for diagonalize=1 only the diagonal entries will be used;
    for values between 0 and 1, the weighted average will be used
    with weight diagonalize for the diagonalized matrix.
    Has no effect for ML estimation.\cr
    Higher values are more stable, lower values potentially more efficient.
    Default: for ML estimation, diagonalize=0; for MoM estimation,
    diagonalize = 1.0.}
  \item{condvarno}{If \code{cond} (conditional simulation), the
    sequential number of the network
    or behavior variable on which to condition.}
  \item{condname}{If conditional, the name of the dependent variable on
    which to condition. Use one or other of \code{condname} or
    \code{condvarno} to specify the variable.}
  \item{firstg}{Initial value of scaling ('gain') parameter for updates
    in the Robbins-Monro procedure.}
  \item{cond}{Boolean. Only relevant for Method of Moments
    simulation/estimation.
    If TRUE, use conditional simulation; if FALSE, unconditional simulation.
	If missing, decision is deferred until \code{\link{siena07}},
	when it is set to TRUE if there is only one dependent variable,
	FALSE otherwise.}
  \item{findiff}{ Boolean: If TRUE, estimate derivatives using finite
    differences. If FALSE, use scores.}
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
    between the observations at start and finish of the wave.
    Decreasing \code{mult} below a certain value has no further effect.}
  \item{simOnly}{Logical: If TRUE, then the calculation of the covariance
    matrix and standard errors of the estimates at the end of
    Phase 3 of the estimation algorithm in function siena07 is skipped.
    This is suitable if nsub=0 and siena07 is used only for the purpose
    of simulation.}
}
\details{Model specification is done via this object for
  \code{\link{siena07}}.
  This function creates an object with the elements required to control the
  Robbins-Monro algorithm. Those not
  available as arguments can be changed manually where desired.
}
\value{Returns an object of class \code{sienaAlgorithm} containing:
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
    \code{sienaModelCreate} from \code{fn} and \code{maxlike}.}
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
  \item{simOnly}{Logical, indicating whether output of covariance matrix
   by \code{\link{siena07}} is to be skipped.}
}
\references{See \url{http://www.stats.ox.ac.uk/~snijders/siena/}}
\author{Ruth Ripley and Tom A.B. Snijders}
\seealso{\code{\link{siena07}}, \code{\link{simstats0c}}.}
\examples{

myAlgorithm <- sienaAlgorithmCreate(projname="NetworkDyn")
StdAlgorithm <- sienaAlgorithmCreate(projname="NetworkDyn", useStdInits=TRUE)
CondAlgorithm <- sienaAlgorithmCreate(projname="NetworkDyn", condvarno=1, cond=TRUE)
Max10Algorithm <- sienaAlgorithmCreate(projname="NetworkDyn", MaxDegree=c(mynet=10))
# where mynet is the name of the network object created by sienaDependent().

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ classes }
