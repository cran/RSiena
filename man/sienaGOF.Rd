\name{sienaGOF}
\alias{sienaGOF}
\alias{plot.sienaGOF}
\title{Functions to assess goodness of fit for SAOMs}
\description{
 The function \code{sienaGOF} assesses goodness of fit
 for a model specification as represented by an estimated
 \code{sienaFit} object. This is done by simulations of auxiliary
 statistics, that differ from the statistics used for estimating
 the parameters. The auxiliary statistics must be given explicitly.
 The fit is good if the average values
 of the auxiliary statistics over many simulation runs are close
 to the values observed in the data.
 A Monte Carlo test based on the Mahalanobis
 distance is used to calculate frequentist \eqn{p}-values.
 Plotting functions can be used to diagnose bad fit. There are basic
 functions for calculating auxiliary statistics available out of the
 box, and the user is also permitted to create custom functions.
 }

\usage{
sienaGOF(sienaFitObject, auxiliaryFunction,
         period=NULL, verbose=FALSE, join=TRUE, twoTailed=FALSE,
         cluster=NULL, robust=FALSE, groupName="Data1",
         varName, \dots)
\method{plot}{sienaGOF}(x, center=FALSE, scale=FALSE, violin=TRUE, key=NULL,
          perc=.05, period=1, main=main, ylab=ylab, \dots)
}
\arguments{
  \item{sienaFitObject}{Results from a call to \code{\link{siena07}}
       with \code{returnDeps = TRUE}.}
  \item{auxiliaryFunction}{Function to be used to calculate
    the auxiliary statistics;
    this can be a user-defined function, e.g. depending on the \code{sna}
    or \code{igraph} packages.

    See \code{Examples} and \code{\link{sienaGOF-auxiliary}}
    for more information on the signature of this function.
    The basic signature is\cr
    \code{function(index, data, sims, period, groupName, varName, \dots)},
    where \code{index} is the index of the simulated network,
    or \code{NULL} if the observed variable is needed;
    \code{data} is the observed data object from which
    the relevant variables are extracted; \code{sims} is the list of simulations
    returned from \code{siena07}; \code{period} is the index of the period;
    and \dots are further arguments (like \code{levls} in the examples below
	and in \code{\link{sienaGOF-auxiliary}}).}
  \item{period}{Vector of period(s) to be used (may run from 1 to
    number of waves - 1). Has an effect only if \code{join=FALSE}.}
  \item{verbose}{Whether to print intermediate results. This may give some
       peace of mind to the user because calculations can take some time.}
  \item{join}{Boolean: should sienaGOF do tests on all of the periods
       individually (\code{FALSE}), or sum across periods (\code{TRUE})?}
  \item{twoTailed}{Whether to use two tails for calculating \eqn{p}-values on the
        Monte Carlo test. Recommended for advanced users only,
        as it is probably only applicable in rare cases.}
  \item{cluster}{Optionally, a \code{snow} cluster to execute the auxiliary
  function calculations on.}
  \item{robust}{ Whether to use robust estimation of the covariance matrix.}
  \item{groupName}{ Name of group; relevant for multi-group data sets.}
  \item{varName}{ Name of dependent variable.}
  \item{x}{ Result from a call to sienaGOF. }
  \item{center}{ Whether to center the statistics by median during plotting.}
  \item{scale}{ Whether to scale the statistics by range during plotting.}
  \item{violin}{ Use violin plots (vs. box plots only)? }
  \item{key}{ Keys in the plot for the levels of the auxiliary statistic
        (as given by parameter \code{levls} in the examples). }
  \item{perc}{ 1 minus confidence level for the confidence bands (two sided). }
  \item{main}{ Main title of the plot. }
  \item{ylab}{ The y-axis label for the plot. }
  \item{\dots}{Other arguments.}
}
\details{
 This function is used to assess the goodness of fit of a stochastic actor
 oriented model for an arbitrarily defined multidimensional auxiliary statistic.
 The auxiliary statistics are calculated for the simulated dependent
 variables in Phase 3 of the estimation algorithm,
 returned in \code{sienaFitObject} because of having used
 \code{returnDeps = TRUE} in the call to \code{\link{siena07}}.
 These statistics should be chosen to represent features of the network
 that are not explicitly fit by the estimation procedure but can be
 considered important properties that the model at hand should represent well.
 Some examples are:
 \itemize{
 \item Outdegree distribution
 \item Indegree distribution
 \item Distribution of the dependent behavior variable (if any).
 \item Distribution of geodesic distances
 \item Triad census
 \item Edgewise homophily counts
 \item Edgewise shared partner counts
 \item Statistics depending on the combination of network
       and behavioral variables.
}
The function is written so that the user can easily define other functions to
capture some other relevant aspects of the network, behaviors, etc.
This is further illustrated in the help page\cr
\code{\link{sienaGOF-auxiliary}}.

We recommend the following heuristic approach to model checking:
\enumerate{
 \item Check convergence of the estimation.
 \item Assess time heterogeneity by \code{\link{sienaTimeTest}}
   and if there is evidence for time heterogeneity
   either modify the base effects or include time dummy terms.
 \item Assess goodness of fit (primarily using \code{join=TRUE}) on
  auxiliary statistics, and if necessary refine the model.
}

 The \code{print} function will display some useful information
 to help with model selection if some effects are set to FIX and
 TEST on the effects object. A rough estimator for the Mahalanobis
 distance that would be obtained at each proposed specification
 is given in the output. This can help guide model selection.
 This estimator is called the modified Mahalanobis distance (MMD).
 See Lospinoso (2012), the manual, or the references for more information.

 The following functions are pre-fabricated for ease of use, and can
 be passed in as the \cr
 \code{auxiliaryFunction} with no extra effort;
 see \code{\link{sienaGOF-auxiliary}} and the examples below.
\itemize{
 \item \code{\link{IndegreeDistribution}}
 \item \code{\link{OutdegreeDistribution}}
 \item \code{\link{BehaviorDistribution}}
}

}
\value{
  \code{sienaGOF} returns a result of class \code{sienaGOF};
    this is a list of elements of class \code{sienaGofTest};
    if \code{join=TRUE}, the list has length 1; if \code{join=FALSE},
    each list element corresponds to a period analyzed;
    the list elements are themselves lists again, including the
    following elements:
    \describe{
    \item{* Observations}{The observed values for the auxiliary statistics.}
    \item{* Simulations}{The simulated auxiliary statistics.}
    \item{* ObservedTestStat}{The observed Mahalobis distance in the data.}
    \item{* SimulatedTestStat}{The Mahalobis distance for the simulations.}
    \item{* TwoTailed}{Whether the \eqn{p}-value corresponds to a one- or
                     two-tailed Monte Carlo test.}
    \item{* p}{The \eqn{p}-value for the observed Mahalanobis distance in the
             permutation distribution of the\cr
			 simulated Mahalanobis distances.}
    \item{* Rank}{ Rank of the covariance matrix of the simulated auxiliary
              statistics.}
    }
}

\references{
\itemize{
\item See \url{http://www.stats.ox.ac.uk/~snijders/siena/}
  for general information on RSiena.
\item Lospinoso, J.A. and Snijders, T.A.B., \dQuote{Goodness of fit for
  Stochastic Actor Oriented Models.} Presentation given at Sunbelt XXXI,
  St. Pete's Beach, Fl. 2011.
\item Lospinoso, J.A. (2012).
  \dQuote{Statistical Models for Social Network Dynamics.}
  Ph.D. Thesis. University of Oxford: U.K.
}
}
\author{Josh Lospinoso, modifications by Ruth Ripley and Tom Snijders}
\seealso{\code{\link{siena07}}, \code{\link{sienaGOF-auxiliary}},
         \code{\link{sienaTimeTest}} }
\examples{
\dontrun{
   mynet1 <- sienaDependent(array(c(s501, s502, s503), dim=c(50, 50, 3)))
   mynet2 <- sienaDependent(array(c(s503, s502, s501), dim=c(50, 50, 3)))
   mybeh <- sienaDependent(s50a, type='behavior')
   mydata <- sienaDataCreate(mynet1, mynet2, mybeh)
   myeff <- getEffects(mydata)
   myeff <- includeEffects(myeff, transTrip)
   myeff <- includeEffects(myeff, recip, name="mynet2")
   myeff <- setEffect(myeff, cycle3, fix=TRUE, test=TRUE, include=TRUE)
   myeff  <- setEffect(myeff, nbrDist2, fix=TRUE, test=TRUE, include=TRUE)
   myeff <- setEffect(myeff, transTies, fix=TRUE, test=TRUE, include=TRUE)
   myalgorithm <- sienaAlgorithmCreate(n3=200) # Shorter phase 3, just for example.
   ans <- siena07(myalgorithm, data=mydata, effects=myeff, returnDeps=TRUE)
   gofi <- sienaGOF(ans, IndegreeDistribution, verbose=TRUE, join=TRUE,
                    varName="mynet1")
   gofi
   plot(gofi)

   gofi2 <- sienaGOF(ans, IndegreeDistribution, verbose=TRUE, join=TRUE,
                     varName="mynet2")
   gofi2
   plot(gofi2)

   gofb <- sienaGOF(ans, BehaviorDistribution, varName = "mybeh",
                    verbose=TRUE, join=TRUE)
   plot(gofb)

   gofo <- sienaGOF(ans, OutdegreeDistribution, verbose=TRUE, join=TRUE,
	                varName="mynet1", cumulative=FALSE)
   # cumulative is an example of "...".
   gofo
   plot(gofo)
}
}
\keyword{models}
