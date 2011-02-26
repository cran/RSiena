#/******************************************************************************
# * SIENA: Simulation Investigation for Empirical Network Analysis
# *
# * Web: http://www.stats.ox.ac.uk/~snidjers/siena
# *
# * File: maxlikec.r
# *
# * Description: This module contains the code for simulating the process,
# * communicating with C++. For use with maximum likelihood method, so
# * never conditional or from finite differences, or parallel testing!
# *****************************************************************************/
##@maxlikec siena07 ML Simulation Module
maxlikec <- function(z, x, INIT=FALSE, TERM=FALSE, initC=FALSE, data=NULL,
                    effects=NULL, profileData=FALSE, prevAns=NULL,
                    returnDeps=FALSE)
{
    if (INIT || initC)  ## initC is to initialise multiple C processes in phase3
    {
        z <- initializeFRAN(z, x, data, effects, prevAns, initC,
                            profileData=profileData, returnDeps=returnDeps)
        if (initC)
        {
            return(NULL)
        }
        else
        {
            return(z)
        }
    }
    if (TERM)
    {
        z <- terminateFRAN(z, x)
        return(z)
    }
    ######################################################################
    ## iteration entry point
    ######################################################################
    ## retrieve stored information
    f <- FRANstore()
    if (z$Phase == 2)
    {
        returnDeps <- FALSE
    }
    else
    {
        returnDeps <- z$returnDeps
    }
    ## create a grid of periods with group names in case want to parallelize
    ## using this
    groupPeriods <- attr(f, "groupPeriods")
    callGrid <- cbind(rep(1:f$nGroup, groupPeriods - 1),
                      as.vector(unlist(sapply(groupPeriods - 1,
                                              function(x) 1:x))))
    ## z$int2 is the number of processors if iterating by period, so 1 means
    ## we are not. Can only parallelize by period at the moment.
    ## browser()
    if (nrow(callGrid) == 1)
    {
        ans <- .Call('mlPeriod', PACKAGE=pkgname, z$Deriv, f$pData,
                     f$pModel, f$myeffects, z$theta,
                     returnDeps, 1, 1, z$nrunMH, z$addChainToStore,
                     z$needChangeContributions)
        ans[[7]] <- list(ans[[7]])
    }
    else
    {
        if (z$int2 == 1)
        {
            anss <- apply(callGrid, 1, doMLModel, z$Deriv, z$theta,
                          returnDeps,  z$nrunMH, z$addChainToStore,
                          z$needChangeContributions)
        }
        else
        {
            use <- 1:(min(nrow(callGrid), z$int2))
            anss <- parRapply(z$cl[use], callGrid, doMLModel, z$Deriv, z$theta,
                              returnDeps, z$nrunMH, z$addChainToStore,
                              z$needChangeContributions)
        }
        ## reorganize the anss so it looks like the normal one
        ans <- NULL
        ans[[1]] <- sapply(anss, "[[", 1) ## statistics
        ans[[2]] <- NULL ## scores
        ans[[3]] <- NULL ## seeds
        ans[[4]] <- NULL ## ntim
        ans[[5]] <- NULL # randomseed
        if (returnDeps) ## chains
        {
            fff <- lapply(anss, "[[", 6)
            fff <- split(fff, callGrid[, 1])
            ans[[6]] <-
                lapply(fff, function(x)
                   {
                       lapply(1:length(f$depNames), function(x, z)
                              lapply(z, "[[", x), z=x)
                   }
                       )
        }
        else
        {
            ans[[6]] <-  NULL
        }
        ans[[7]] <- lapply(anss, "[[", 7) ## derivative
        ans[[8]] <- lapply(anss, "[[", 8)
        ans[[8]] <- do.call("+",  ans[[8]]) ## accepts
        ans[[9]] <- lapply(anss, "[[", 9)
        ans[[9]] <- do.call("+",  ans[[9]]) ## rejects
    }

    dff <- ans[[7]]
    fra <- -t(ans[[1]]) ##note sign change

    FRANstore(f)

    if (returnDeps)
        sims <- ans[[6]]
    else
        sims <- NULL

    #print(length(sims[[1]]))
    if (returnDeps) ## this is not finished yet!
    {
        ## attach the names
        names(sims) <- f$groupNames
        periodNo <- 1
        for (i in 1:length(sims))
        {
            names(sims[[i]]) <- f$depNames
            for (j in 1:length(sims[[i]]))
            {
                periodNos <- periodNo:(periodNo  + length(sims[[i]][[j]]) - 1)
                names(sims[[i]][[j]]) <- periodNos
            }
            periodNo <- periodNos[length(periodNos)] + 2
       }
    }
    if (z$Deriv) ## need to reformat the derivatives
    {
        ## current format is a list of vectors of the lower? triangle
        ## of the matrix. Need to be put into a symmetric matrix.
        ## Tricky part is getting the rates in the right place!
        ## Note that we do not yet deal with rate effects other than the basic
        dff <- matrix(0, nrow=z$pp, ncol=z$pp)
        nPeriods <- length(ans[[7]])
        dff2 <- array(0, dim=c(nPeriods, z$pp, z$pp))
        for (period in 1:nPeriods)
        {
            dffraw <- ans[[7]][[period]]
            ## start indexes row/col for first effect for the variable
            start <- 1
            ## rawsub is subscript in the vector
            rawsub <- 1
            ## f$myeffects is a list of data frames per dependent variable
            ## of selected effects.
            for (i in 1:length(f$myeffects))
            {
                dffPeriod <- matrix(0, nrow=z$pp, ncol=z$pp)
                ## rows is the number of effects for this variable
                rows <- nrow(f$myeffects[[i]])
                ## nRates is the number of rate effects for this variable
                ## at present nRates will be the number of periods.
                nRates <- sum(f$myeffects[[i]]$type == 'rate')
                ## nonRates is the number of rows/cols in the objective function
                ## part of the derivative matrix dff.
                nonRates <- rows - nRates
                ## first put the basic rate for this variable in the right place
                dffPeriod[cbind(start:(start + nRates -1),
                                start:(start + nRates - 1))] <-
                                    dffraw[rawsub:(rawsub + nRates - 1)]
                ##
                ## now the matrix of objective function effects
                ##
                start <- start + nRates
                ##
                rawsub <- rawsub + nRates
                ##
                dffPeriod[start : (start + nonRates - 1 ),
                          start : (start + nonRates - 1)] <-
                              dffraw[rawsub:(rawsub + nonRates * nonRates - 1)]
                start <- start + nonRates
                rawsub <- rawsub + nonRates * nonRates
                dffPeriod <- dffPeriod + t(dffPeriod)
                diag(dffPeriod) <- diag(dffPeriod) / 2
                dff2[period , , ] <- dff2[period, , ] - dffPeriod
                dff <- dff - dffPeriod
            }
        }
    }
    else
    {
        dff <- NULL
        dff2 <- NULL
    }
    ## browser()

    list(fra = fra, ntim0 = NULL, feasible = TRUE, OK = TRUE,
         sims=sims, dff = dff, dff2=dff2,
         chain = list(ans[[6]]), accepts=ans[[8]],
         rejects= ans[[9]])
}

dist2full<-function(dis) {

      n<-attr(dis,"Size")

      full<-matrix(0,n,n)

      full[lower.tri(full)]<-dis

      full+t(full)

}
##@doMLModel Maximum likelihood
doMLModel <- function(x, Deriv, theta, returnDeps, nrunMH, addChainToStore,
                      needChangeContributions)
{
    f <- FRANstore()
    .Call("mlPeriod", PACKAGE=pkgname, Deriv, f$pData,
          f$pModel, f$myeffects, theta, returnDeps,
          as.integer(x[1]), as.integer(x[2]), nrunMH, addChainToStore,
          needChangeContributions)
}
