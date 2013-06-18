\name{sienaGOF-auxiliary}
\alias{sienaGOF-auxiliary}
\alias{sparseMatrixExtraction}
\alias{networkExtraction}
\alias{behaviorExtraction}
\alias{OutdegreeDistribution}
\alias{IndegreeDistribution}
\alias{BehaviorDistribution}

\title{Auxiliary functions for goodness of fit assessment by \code{\link{sienaGOF}}}
\description{
 The functions given here are auxiliary to function \code{\link{sienaGOF}}
 which assesses goodness of fit for actor-oriented models.

 The auxiliary functions are, first, some functions of networks or behavior
 (i.e., statistics) for which the simulated values for the fitted model
 are compared to the observed value; second, some extraction functions
 to extract the observed and simulated networks and/or behavior
 from the \code{\link{sienaFit}} object produced by \code{\link{siena07}}
 with \code{returnDeps=TRUE}.

 These functions are exported here mainly to enable users to write
 their own versions. At the end of this help page some non-exported functions
 are listed. These are not exported because they depend on packages that
 are not in the R base distribution; and to show templates for readers
 wishing to contruct their own functions.
 }

\usage{
OutdegreeDistribution(i, obsData, sims, period, groupName, varName,
         levls=0:8, cumulative=TRUE)

IndegreeDistribution(i, obsData, sims, period, groupName, varName,
         levls=0:8, cumulative=TRUE)

BehaviorDistribution(i, obsData, sims, period, groupName, varName,
         levls=NULL, cumulative=TRUE)

sparseMatrixExtraction(i, obsData, sims, period, groupName, varName)

networkExtraction(i, obsData, sims, period, groupName, varName)

behaviorExtraction(i, obsData, sims, period, groupName, varName)
}
\arguments{
  \item{i}{ Index number of simulation to be extracted, ranging from 1 to
    \code{length(sims)}; if NULL, the data observation will be extracted.}
  \item{obsData}{ The observed data set to which the model was fitted; normally
    this is \code{x$f} where \code{x} is the \code{\link{sienaFit}}
    object for which the fit is being assessed.}
  \item{sims}{ The simulated data sets to be compared with the observed data;
    normally this is
    \code{x$sims} where \code{x} is the \code{\link{sienaFit}}
    object for which the fit is being assessed.}
  \item{period}{ Period for which data and simulations are used (may run from
    1 to number of waves - 1).}
  \item{groupName}{ Name of group; relevant for multi-group data sets;
    defaults in \code{\link{sienaGOF}} to "\code{Data1}".}
  \item{varName}{ Name of dependent variable. }
  \item{levls}{ Levels used as values of the auxiliary statistic.
    For \code{BehaviorDistribution}, this defaults to the observed
    range of values.}
  \item{cumulative}{ Are the distributions to be considered as raw or
    cumulative (<=) distributions?}
}
\details{
 The statistics should be chosen to represent features of the network
 that are not explicitly fit by the estimation procedure but can be
 considered important properties that the model at hand should represent well.
 The three given here are far from a complete set;
 they will be supplemented in due time by statistics depending on networks and
 behavior jointly.

 The method signature for the auxiliary statistics generally is\cr
 \code{function(i, obsData, sims, period, groupName, varName, \dots)}.
 For constructing new auxiliary statistics, it is helpful to
 study the code of \code{OutdegreeDistribution}, \code{IndegreeDistribution},
 and \code{BehaviorDistribution} and of the example functions below.
 }

\value{
 \code{OutdegreeDistribution} returns a named vector, the distribution
  of the observed or simulated outdegrees for the values in \code{levls}.

 \code{IndegreeDistribution} returns a named vector, the distribution
  of the observed or simulated indegrees for the values in \code{levls}.

 \code{BehaviorDistribution} returns a named vector, the distribution of the
  observed or simulated behavioral variable for the values in \code{levls}.

 \code{sparseMatrixExtraction} returns the simulated network as a
  \code{dgCMatrix};
  this is the "standard" class for sparse numeric matrices in the
  \code{Matrix} package. See the help file for \code{\link{dgCMatrix-class}}.
  Tie variables for ordered pairs with a missing value for wave=period or
  period+1 are zeroed;
  note that this also is done in RSiena for calculation of target statistics.\cr
  To treat the objects returned by this function as regular matrices,
  it is necessary to attach the \code{Matrix} package in your session.

\code{networkExtraction} returns the network as an edge list of class "network"
  according to the network package (used for package sna).
  Tie variables for ordered pairs with a missing value for wave=period or
  period+1 are zeroed;
  note that this also is done in RSiena for calculation of target statistics.

\code{behaviorExtraction} returns the dependent behavior
  variable as an integer vector.
  Values for actors with a missing value for wave=period or period+1 are
  transformed to NA.
}

\references{
\itemize{
\item See \url{http://www.stats.ox.ac.uk/~snijders/siena/}
  for general information on RSiena.
\item Lospinoso, J.A. and Snijders, T.A.B., \dQuote{Goodness of fit for
Stochastic Actor Oriented Models.} Presentation given at Sunbelt XXXI,
 St. Pete's Beach, Fl. 2011.
\item Lospinoso, J.A. (2012). \dQuote{Statistical Models for Social Network Dynamics.}
Ph.D. Thesis. University of Oxford: U.K.
}
}
\author{Josh Lospinoso, Tom Snijders}
\seealso{\code{\link{siena07}}, \code{\link{sienaGOF}} }
\examples{
\dontrun{
   ### For use out of the box:

   mynet1 <- sienaDependent(array(c(s501, s502, s503), dim=c(50, 50, 3)))
   mybeh <- sienaDependent(s50a, type='behavior')
   mydata <- sienaDataCreate(mynet1, mybeh)
   myeff <- getEffects(mydata)
   myeff <- includeEffects(myeff, transTies, cycle3)
   # Shorter phases 2 and 3, just for example:
   myalgorithm <- sienaAlgorithmCreate(nsub=2, n3=300)
   ans <- siena07(myalgorithm, data=mydata, effects=myeff, returnDeps=TRUE)

   OutdegreeDistribution(NULL, ans$f, ans$sims, period=1, groupName="Data1",
                    levls=0:7, varName="mynet1")
   IndegreeDistribution(5, ans$f, ans$sims, period=1, groupName="Data1",
                    varName="mynet1")
   BehaviorDistribution(20, ans$f, ans$sims, period=1, groupName="Data1",
                    varName="mybeh")
   sparseMatrixExtraction(50, ans$f, ans$sims, period=1, groupName="Data1",
                    varName="mynet1")
   networkExtraction(100, ans$f, ans$sims, period=1, groupName="Data1",
                    varName="mynet1")
   behaviorExtraction(200, ans$f, ans$sims, period=1, groupName="Data1",
                    varName="mybeh")

   gofi <- sienaGOF(ans, IndegreeDistribution, verbose=TRUE, join=TRUE,
                    varName="mynet1")
   gofi
   plot(gofi)

   (gofo <- sienaGOF(ans, OutdegreeDistribution, verbose=TRUE, join=TRUE,
                    varName="mynet1", cumulative=FALSE))
   # cumulative is an example of "\dots".
   plot(gofo)

   (gofb <- sienaGOF(ans, BehaviorDistribution, varName = "mybeh",
                    verbose=TRUE, join=TRUE, cumulative=FALSE))
   plot(gofb)

   ### Here come some useful functions for building your own auxiliary statistics:
   ### First an extraction function.

   # igraphNetworkExtraction extracts simulated and observed networks
   # from the results of a siena07 run.
   # It returns the network as an edge list of class "graph"
   # according to the igraph package.
   # Ties for ordered pairs with a missing value for wave=period or period+1
   # are zeroed;
   # note that this also is done in RSiena for calculation of target statistics.
   igraphNetworkExtraction <- function(i, data, sims, period, groupName, varName){
     require(igraph)
     dimsOfDepVar<- attr(data[[groupName]]$depvars[[varName]], "netdims")
     missings <- is.na(data[[groupName]]$depvars[[varName]][,,period]) |
                 is.na(data[[groupName]]$depvars[[varName]][,,period+1])
     if (is.null(i)) {
   # sienaGOF wants the observation:
       original <- data[[groupName]]$depvars[[varName]][,,period+1]
       original[missings] <- 0
       returnValue <- graph.adjacency(original)
     }
     else
     {
       missings <- graph.adjacency(missings)
   #sienaGOF wants the i-th simulation:
       returnValue <- graph.difference(
                graph.edgelist(sims[[i]][[groupName]][[varName]][[period]][,1:2]),
                missings)
     }
     returnValue
   }

   ### Then some auxiliary statistics.

   # GeodesicDistribution calculates the distribution of directed
   # geodesic distances; see ?sna::geodist
   # The default for \code{levls} reflects that geodesic distances larger than 5
   # do not differ appreciably with respect to interpretation.
   # Note that the levels of the result are named;
   # these names are used in the \code{plot} method.
   GeodesicDistribution <- function (i, data, sims, period, groupName,
                           varName, levls=c(1:5,Inf), cumulative=TRUE, ...) {
     x <- networkExtraction(i, data, sims, period, groupName, varName)
     require(sna)
     a <- sna::geodist(x)$gdist
     if (cumulative)
     {
       gdi <- sapply(levls, function(i){ sum(a<=i) })
     }
	 else
     {
       gdi <- sapply(levls, function(i){ sum(a==i) })
     }
     names(gdi) <- as.character(levls)
     gdi
   }

   # Holland and Leinhardt Triad Census; see ?sna::triad.census.
   TriadCensus <- function(i, data, sims, wave, groupName, varName, levls=1:16){
       unloadNamespace("igraph") # to avoid package clashes
       require(sna)
       require(network)
       x <- networkExtraction(i, data, sims, wave, groupName, varName)
       tc <- sna::triad.census(x)[1,levls]
       # names are transferred automatically
       tc
   }

  # Distribution of Burt's constraint values; see ?igraph::constraint
  # the maximum finite value is 9/8 (see Buskens and van de Rijt, AJS 2008).
  ConstraintDistribution <- function (i, data, sims, period, groupName, varName,
                           levls=c(seq(0,1.125,by=0.125)), cumulative=TRUE){
     require(igraph)
     x <- igraphNetworkExtraction(i, data, sims, period, groupName, varName)
     a <- igraph::constraint(x)
     a[is.na(a)] <- Inf
     lel <- length(levls)
     if (cumulative)
     {
       cdi <- sapply(2:lel, function(i){sum(a<=levls[i])})
     }
     else
     {
       cdi <- sapply(2:lel, function(i){
                     sum(a<=levls[i]) - sum(a <= levls[i-1])})
     }
     names(cdi) <- as.character(levls[2:lel])
     cdi
    }

  ## Finally some examples of the three auxiliary statistics constructed above.

   mynet1 <- sienaDependent(array(c(s501, s502, s503), dim=c(50, 50, 3)))
   mybeh <- sienaDependent(s50a, type='behavior')
   mydata <- sienaDataCreate(mynet1, mybeh)
   myeff <- getEffects(mydata)
   myeff <- includeEffects(myeff, transTrip, cycle3, nbrDist2)
   myeff <- includeEffects(myeff, outdeg, name="mybeh",
            interaction1="mynet1")
   myeff <- includeEffects(myeff,  outdeg, name="mybeh",
            interaction1="mynet1")
   # Shorter phases 2 and 3, just for example:
   myalgorithm <- sienaAlgorithmCreate(nsub=2, n3=300)
   (ans2 <- siena07(myalgorithm, data=mydata, effects=myeff, returnDeps=TRUE))

   gofc <- sienaGOF(ans2, ConstraintDistribution, varName="mynet1",
           verbose=TRUE, join=TRUE)
   plot(gofc)

   goftc <- sienaGOF(ans2, TriadCensus, varName="mynet1", verbose=TRUE, join=TRUE)
   plot(goftc, center=TRUE, scale=TRUE)
   # For this type of auxiliary statistics
   # it is advisable in the plot to center and scale.
   # note the keys at the x-axis; these names are given by sna::triad.census

   gofgd <- sienaGOF(ans2, GeodesicDistribution, varName="mynet1",
            verbose=TRUE, join=TRUE, cumulative=FALSE)
   plot(gofgd)
   # and without infinite distances:
   gofgdd <- sienaGOF(ans2, GeodesicDistribution, varName="mynet1",
             verbose=TRUE, join=TRUE, levls=1:7, cumulative=FALSE)
   plot(gofgdd)
  }
}
\keyword{models}
