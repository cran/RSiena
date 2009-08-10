##/*****************************************************************************
## * SIENA: Simulation Investigation for Empirical Network Analysis
## *
## * Web: http://www.stats.ox.ac.uk/~snijders/siena
## *
## * File: sienaprint.r
## *
## * Description: This file contains the print and summary modules for the
## * classes siena, sienaFit and sienaModel
## *
## ****************************************************************************/

print.siena <- function(x, ...)
{
  if (!inherits(x, "siena"))
        stop("not a legitimate Siena data object")
  cat('Dependent variables: ', names(x$depvars),'\n')
  cat('Number of waves:', x$observations)
  if (!is.null(x$nodesets))
  {
      tmp <- cbind(c('Nodesets',
                     names(x$nodesets)),
                   c('Number of nodes',
                     sapply(x$nodesets,length)))
      print(tmp)
  }
}
print.sienaGroup <- function(x, ...)
{
  if (!inherits(x, "sienaGroup"))
        stop("not a legitimate Siena group data object")
  att <- attributes(x)
  cat('Dependent variables: \n')
  cat(paste(att$netnames, ":", att$types),'\n')
  cat('Total number of periods:', att$observations)
  cat("\nmore to be added!\n")
}

print.sienaFit <- function(x, ...)
{
   if (!inherits(x, "sienaFit"))
        stop("not a legitimate Siena model fit")
   if (!x$OK)
   {
       cat("Error end of estimation algorithm")
   }
   else
       {
           cat("Estimates and standard errors\n\n")
           ses <- format(sqrt(diag(x$covtheta)),nsmall=4,digits=4,width=9)
           ses[x$fixed] <- ' (fixed) '
           cat(paste(1:length(x$effects$effectName),
                     format(x$effects$type),
                     format(x$effects$effectName),
                     format(x$theta,nsmall=4,digits=4,width=10),
                     ' (',ses,')\n',collapse=''))
           cat("\nTotal of",x$n,"iteration steps.\n\n")
           if (x$termination =="UserInterrupt")
               cat(" \n*** Warning ***",
                   "Estimation terminated early at user request.\n")
       }
}

summary.sienaFit <- function(object, ...)
{
    class(object) <- "summary.sienaFit"
    object
}

print.summary.sienaFit <- function(x, ...)
{
   if (!inherits(x, "summary.sienaFit"))
        stop("not a legitimate summary of a Siena model fit")
   if (!x$OK)
   {
       cat("Error end of estimation algorithm")
   }
   else
       {
           cat("Estimates and standard errors\n\n")
           ses <- format(sqrt(diag(x$covtheta)),nsmall=4,digits=4,width=9)
           ses[x$fixed] <- ' (fixed) '
           cat(paste(1:length(x$effects$effectName),
                     format(x$effects$type),
                     format(x$effects$effectName),
                     format(x$theta,nsmall=4,digits=4,width=10),
                     ' (',ses,')\n',collapse=''))
           cat("Total of",x$n,"iteration steps.\n\n")
           cat("Covariance matrix of estimates (correlations below diagonal)\n\n")
           covcor <- x$covtheta
           correl <- x$covtheta/sqrt(diag(x$covtheta))[row(x$covtheta)]/
               sqrt(diag(x$covtheta))[col(x$covtheta)]
           covcor[lower.tri(covcor)] <- correl[lower.tri(correl)]
           printMatrix(format(round(t(covcor),digits=3),width=12))
           cat("\nDerivative matrix of expected statistics X by parameters:\n\n")
           printMatrix(format(round(x$dfra,digits=3),width=12))
              cat("\nCovariance matrix of X (correlations below diagonal):\n\n")
           covcor <- x$msf
           correl <- x$msf/sqrt(diag(x$msf))[row(x$msf)]/sqrt(diag(x$msf))[col(x$msf)]
           covcor[lower.tri(covcor)] <- correl[lower.tri(correl)]
           printMatrix(format(round(t(covcor),digits=3),width=12))
      }

}


printMatrix<- function(mat)
{
    cat(mat,sep=c(rep.int(' ',ncol(mat)-1),'\n'))
}

print.sienaModel <- function(x, ...)
{
    cat('\n Project name:', x$projname, '\n')
    cat('\n Use standard initial values:', x$useStdInits, '\n')
    cat(' Random seed:',x$randomSeed,'\n')
    cat(' Starting value of gain parameter', x$startg, '\n')
    if (any(x$MaxDegree > 0))
    {
        cat(' Restrictions on degree in simulations:')
        cat(x$MaxDegree,'\n')
    }
    cat(' Method for calculation of Derivatives:',
        c('Finite Differences', 'Scores')[as.numeric(x$FinDiff.method)], '\n')
    cat(' Number of subphases in phase 2:', x$nsub, '\n')
    cat(' Number of iterations in phase 3:', x$n3, '\n')
    if (!x$cconditional)
    {
        cat('Unconditional simulation\n')
    }
    else
    {
        cat('Conditional simulation:')
        if (x$condname != '')
            cat('conditioned on', x$condname,'\n')
        else
            if (x$condvarno > 0)
                cat('conditioned on First variable')
    }

}
