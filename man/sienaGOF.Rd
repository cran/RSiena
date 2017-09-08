\name{sienaGOF}
\alias{sienaGOF}
\alias{plot.sienaGOF}
\alias{descriptives.sienaGOF}
\title{Functions to assess goodness of fit for SAOMs}
\description{
 The function \code{sienaGOF} assesses goodness of fit
 for a model specification as represented by an estimated
 \code{sienaFit} object created by \code{\link{siena07}}.
 This is done by simulations of auxiliary
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
        perc=.05, period=1, \dots)
descriptives.sienaGOF(x, center=FALSE, scale=FALSE, perc=.05, key=NULL,
        period=1, showAll=FALSE)
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
  \item{scale}{ Whether to scale the statistics by range during plotting.
                scale=TRUE makes little sense without also center=TRUE.}
  \item{violin}{ Use violin plots (vs. box plots only)? }
  \item{key}{ Keys in the plot for the levels of the auxiliary statistic
        (as given by parameter \code{levls} in the examples). }
  \item{perc}{ 1 minus confidence level for the confidence bands (two sided). }
  \item{\dots}{Other arguments; e.g. (1), \code{levls} for \code{sienaGOF()} where
      \code{levls} is a parameter for the auxiliary statistic as in
      \code{\link{sienaGOF-auxiliary}}; e.g. (2), \code{main}, \code{xlab}
      and/or \code{ylab} for \code{plot.sienaGOF()}.}
  \item{showAll}{If FALSE, drops statistics with variance 0, like in the plot.}
}
\details{
 This function is used to assess the goodness of fit of a stochastic actor
 oriented model for an arbitrarily defined multidimensional auxiliary statistic.
 It operates basically by comparing the observed values, at the ends of
 the periods, with the simulated values for the ends of the periods.
 The differences are assessed by combining the components of the auxiliary
 statistic using the Mahalanobis distance.

 To achieve comparability between simulated and observed dependent variables,
 variables that are missing in the data at the start or end of a period
 are replaced by 0 (for tie variables) or NA (for behavior variables).
 Tie variables that are structurally determined at the beginning of a period
 are used to replace observed values at the end of the period;
 tie variables that are structurally determined at the end, but not
 the beginning, of a period are used to replace
 simulated values at the end of the period.

 The auxiliary statistics documented in \code{\link{sienaGOF-auxiliary}}
 are calculated for the simulated dependent
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
    \item{- sienaFitName}{The name of sienaFitObject.}
    \item{- auxiliaryStatisticName}{The name of auxiliaryFunction.}
    \item{- Observations}{The observed values for the auxiliary statistics.}
    \item{- Simulations}{The simulated auxiliary statistics.}
    \item{- ObservedTestStat}{The observed Mahalanobis distance in the data.}
    \item{- SimulatedTestStat}{The Mahalanobis distance for the simulations.}
    \item{- TwoTailed}{Whether the \eqn{p}-value corresponds to a one- or
                     two-tailed Monte Carlo test.}
    \item{- p}{The \eqn{p}-value for the observed Mahalanobis distance in the
             permutation distribution of the
             simulated Mahalanobis distances.}
    \item{- Rank}{ Rank of the covariance matrix of the simulated auxiliary
              statistics.}
    In addition there are several attributes which give, for model
    specifications with fixed-and-tested effects, approximations
    to the expected Mahalanobis distance for model specifications where
    each of these effects would be added. This is reported
    in the \code{summary} method. \cr
    The \code{plot} method makes violin plots or box plots, with superimposed
    confidence bands, for the simulated distributions of all elements
    of the \code{auxiliaryFunction},
    with the observed values indicated by red dots;
    but statistics with variance 0 are dropped.

    \code{descriptives.sienaGOF} returns a matrix giving numerical
    information about what is plotted in the \code{plot} method:
    maximum, upper percentile, mean, median, lower percentile, minimum,
    of the simulated distributions of the auxiliary statistics, and the
    observed values. If \code{center=TRUE} the median is subtracted
    from mean, median, and percentiles; if \code{scale=TRUE}
    these numbers are divided by (maximum - minimum). \cr
    If \code{showAll=FALSE}, statistics with variance 0 will be dropped.
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
   mynet <- sienaDependent(array(c(s501, s502), dim=c(50, 50, 2)))
   mybeh <- sienaDependent(s50a[,1:2], type="behavior")
   mydata <- sienaDataCreate(mynet, mybeh)
   myeff <- getEffects(mydata)
   myeff <- includeEffects(myeff, transTrip)
   myeff <- setEffect(myeff, cycle3, fix=TRUE, test=TRUE)
   myeff <- setEffect(myeff, transTies, fix=TRUE, test=TRUE)
   myalgorithm <- sienaAlgorithmCreate(nsub=1, n3=50)
                   # Shorter phases 2 and 3, just for example.
   ans <- siena07(myalgorithm, data=mydata, effects=myeff, batch=TRUE, returnDeps=TRUE)
   gofi <- sienaGOF(ans, IndegreeDistribution, verbose=TRUE, join=TRUE,
                    varName="mynet")
   summary(gofi)
   plot(gofi)

\dontrun{
# Illustration just for showing a case with two dependent networks;
# running time backwards is not meaningful!
   mynet1 <- sienaDependent(array(c(s501, s502), dim=c(50, 50, 2)))
   mynet2 <- sienaDependent(array(c(s503, s501), dim=c(50, 50, 2)))
   mybeh <- sienaDependent(s50a[,1:2], type="behavior")
   mydata <- sienaDataCreate(mynet1, mynet2, mybeh)
   myeff <- getEffects(mydata)
   myeff <- includeEffects(myeff, transTrip)
   myeff <- includeEffects(myeff, recip, name="mynet2")
   myeff <- setEffect(myeff, cycle3, fix=TRUE, test=TRUE)
   myeff <- setEffect(myeff, transTies, fix=TRUE, test=TRUE)
   ans <- siena07(myalgorithm, data=mydata, effects=myeff, batch=TRUE, returnDeps=TRUE)
   gofi <- sienaGOF(ans, IndegreeDistribution, verbose=TRUE, join=TRUE,
                    varName="mynet1")
   summary(gofi)
   plot(gofi)

   (gofi.nc <- sienaGOF(ans, IndegreeDistribution, cumulative=FALSE,
                        varName="mynet1"))
   # cumulative is an example of "...".
   plot(gofi.nc)
   descriptives.sienaGOF(gofi.nc)

   (gofi2 <- sienaGOF(ans, IndegreeDistribution, varName="mynet2"))
   plot(gofi2)

   (gofb <- sienaGOF(ans, BehaviorDistribution, varName = "mybeh"))
   plot(gofb)

   (gofo <- sienaGOF(ans, OutdegreeDistribution, varName="mynet1",
                        levls=0:6, cumulative=FALSE))
   # levls is another example of "...".
   plot(gofo)
}

## A demonstration of using multiple processes
\dontrun{
   library(parallel)
   (n.clus <- detectCores() - 1) # subtract 1 to keep time for other processes
   myalgorithm.c <- sienaAlgorithmCreate(nsub=4, n3=1000, seed=12345)
   (ans.c <- siena07(myalgorithm.c, data=mydata, effects=myeff, batch=TRUE,
                        returnDeps=TRUE, useCluster=TRUE, nbrNodes=n.clus))
   gofi.1 <- sienaGOF(ans.c, IndegreeDistribution, verbose=TRUE,
                    varName="mynet1")
   cl <- makeCluster(n.clus)
   gofi.cl <- sienaGOF(ans.c, IndegreeDistribution, varName="mynet1",
                    cluster=cl)
   # compare simulation times
   attr(gofi.1,"simTime")
   attr(gofi.cl,"simTime")
}

}
\keyword{models}
