#/******************************************************************************
# * SIENA: Simulation Investigation for Empirical Network Analysis
# *
# * Web: http://www.stats.ox.ac.uk/~snidjers/siena
# *
# * File: simstatsc.r
# *
# * Description: This module contains the code for simulating the process,
# * communicating with C++. Only subsidiary routines used for maximum likelihood
# *****************************************************************************/
##@simstats0c siena07 Simulation Module
simstats0c <- function(z, x, INIT=FALSE, TERM=FALSE, initC=FALSE, data=NULL,
                       effects=NULL, fromFiniteDiff=FALSE,
                       profileData=FALSE, prevAns=NULL, returnDeps=FALSE,
                       returnChains=FALSE)
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
    ## ####################################################################
    ## iteration entry point: not for ML
    ## ####################################################################
    ## retrieve stored information
    f <- FRANstore()
    ## browser()
    ## fix up the interface so can call from outside robmon framework
    if (is.null(z$Deriv))
    {
        z$Deriv <- FALSE
    }
    if (is.null(z$Phase))
    {
        z$Phase <- 1 ### nb be aware
    }
    if (fromFiniteDiff)
    {
        returnDeps <- FALSE
    }
    else
    {
        returnDeps <- z$returnDeps
    }
    if (is.null(f$seeds))
    {
        seeds <- NULL
    }
    else
    {
        seeds <- f$seeds
    }
    if (is.null(f$randomseed2))
    {
        randomseed2 <- NULL
    }
    else
    {
        if (fromFiniteDiff)
        {
            randomseed2 <- as.integer(f$storedseed)
        }
        else
        {
            randomseed2 <- as.integer(f$randomseed2)
            f$storedseed <- randomseed2
        }
        ## cat(randomseed2, '\n')
    }
    ## create a grid of periods with group names in case want to parallelize
    ## using this
    groupPeriods <- attr(f, "groupPeriods")
    callGrid <- cbind(rep(1:f$nGroup, groupPeriods - 1),
                      as.vector(unlist(sapply(groupPeriods - 1,
                                              function(x) 1:x))))
    ## z$int2 is the number of processors if iterating by period, so 1 means
    ## we are not
    if (z$int2==1 || nrow(callGrid) == 1)
    {
        ##   cat("theta", z$theta, "\n")
        ans <- .Call('model', PACKAGE=pkgname, z$Deriv, f$pData, seeds,
                     fromFiniteDiff, f$pModel, f$myeffects, z$theta,
                     randomseed2, returnDeps, z$FinDiff.method,
                     !is.null(z$cl), z$addChainToStore,
                     z$needChangeContributions, returnChains)
    }
    else
    {
        use <- 1:(min(nrow(callGrid), z$int2))
        anss <- parRapply(z$cl[use], callGrid, doModel,
                          z$Deriv, seeds, fromFiniteDiff, z$theta,
                          randomseed2, returnDeps, z$FinDiff.method, TRUE,
                          returnChains)
        ## reorganize the anss so it looks like the normal one
        ## browser()
        ans <- NULL
        ans[[1]] <- sapply(anss, "[[", 1) ## statistics
        ans[[2]] <- sapply(anss, "[[", 2) ## scores
        ans[[3]] <- split(lapply(anss, "[[", 3), callGrid[, 1]) ## seeds
        ans[[4]] <- sapply(anss, "[[", 4) # ntim
        ans[[5]] <- NULL # randomseed not sensible here
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
    ## browser()
    if (!fromFiniteDiff)
    {
        if (z$FinDiff.method)
            f$seeds <- ans[[3]]
    }
    if (z$Deriv )
    {
        sc <- t(ans[[2]])
    }
    else
    {
        sc <-  NULL
    }
    ntim <- ans[[4]]
    fra <- t(ans[[1]])
    f$randomseed2 <- ans[[5]]#[c(1,4,3,2)]
    FRANstore(f)
    if (returnDeps)
    {
        sims <- ans[[6]]
    }
    else
    {
        sims <- NULL
    }
    if (returnChains)
    {
        chain <- ans[[7]]
    }
    else
    {
        chain <- NULL
    }
    if (returnDeps)
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
                                       ## browser()
    list(sc = sc, fra = fra, ntim0 = ntim, feasible = TRUE, OK = TRUE,
         sims=sims, f$seeds, chain=chain)
}
doModel <- function(x, Deriv, seeds, fromFiniteDiff, theta, randomseed2,
                    returnDeps, FinDiff.method, useStreams, returnChains)
{
    f <- FRANstore()
    seeds <- seeds[[x[1]]][[x[2]]]
    .Call("modelPeriod", PACKAGE=pkgname, Deriv, f$pData, seeds,
          fromFiniteDiff, f$pModel, f$myeffects, theta,
          randomseed2, returnDeps, FinDiff.method, useStreams,
          as.integer(x[1]), as.integer(x[2]), returnChains)
}

##@clearData siena07 Finalizer to clear Data object in C++
clearData <- function(pData)
{
    .Call('deleteData', PACKAGE=pkgname, pData)
}
##@clearModel siena07 Finalizer to clear Model object in C++
clearModel <- function(pModel)
{
    .Call('deleteModel', PACKAGE=pkgname, pModel)
}
