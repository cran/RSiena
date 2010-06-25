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
        f <- FRANstore()
        f$pModel <- NULL
        f$pData <- NULL
        FRANstore(NULL) ## clear the stored object
        if (is.null(z$print) || z$print)
        {
            PrintReport(z, x)
        }
        if (sum(z$test))
        {
            z$fra <- colMeans(z$sf, na.rm=TRUE)
            ans <- ScoreTest(z$pp, z$dfra, z$msf, z$fra, z$test, x$maxlike)
            z <- c(z, ans)
            TestOutput(z, x)
        }
        if (!is.null(z$dfra))
        {
            dimnames(z$dfra)[[1]] <- as.list(z$requestedEffects$shortName)
        }
        return(z)
    }
    ######################################################################
    ## iteration entry point
    ######################################################################
    ## retrieve stored information
    f <- FRANstore()
    ## browser()
    ##if (z$Phase == 2)
    ##{
    ##    returnDeps <- FALSE
    ##}
    ##else
    ##{
    returnDeps <- z$returnDeps
    ##}
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
        ## for now!
        ans <- .Call('mlPeriod', PACKAGE=pkgname, z$Deriv, f$pData,
                     f$pModel, f$myeffects, f$pMLSimulation, z$theta,
                     returnDeps, 1, 1, z$nrunMH)
    }
    else
    {
        use <- 1:(min(nrow(callGrid), z$int2))
      #  anss <- parRapply(z$cl[use], callGrid, doModel,
      #                    z$Deriv, seeds, , z$theta,
      #                    randomseed2, returnDeps, z$FinDiff.method, TRUE)
      #  ##anss <- apply(callGrid, 1, doModel,
        ##            z$Deriv, fromFiniteDiff, z$theta,
        ##           returnDeps, z$FinDiff.method)
        ## reorganize the anss so it looks like the normal one
        ## browser()
        ans <- NULL
      #  ans[[1]] <- sapply(anss, "[[", 1) ## statistics
       # ans[[2]] <- sapply(anss, "[[", 2) ## scores
      #  ans[[3]] <- split(lapply(anss, "[[", 3), callGrid[, 1]) ## seeds
      #  ans[[4]] <- sapply(anss, "[[", 4) # ntim
      #  ans[[5]] <- NULL # randomseed not sensible here
      #  fff <- lapply(anss, "[[", 6)
      #  fff <- split(fff, callGrid[, 1])
     #   ans[[6]] <-
     #       lapply(fff, function(x)
      #         {
       #            lapply(1:length(f$depNames), function(x, z)
        #                  lapply(z, "[[", x), z=x)
         #      }
          #         )
    }
    ## browser()
    dff <- ans[[7]]
    fra <- -t(ans[[1]]) ##note sign change

    FRANstore(f)

    if (returnDeps)
        sims <- ans[[6]]
    else
        sims <- NULL
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
        dffraw <- ans[[7]]
        dff <- matrix(0, nrow=z$pp, ncol=z$pp)
        start <- 1
        rawsub <- 1
        for (i in 1:length(f$myeffects))
        {
            rows <- nrow(f$myeffects[[i]])
            nRates <- sum(f$myeffects[[i]]$type == 'rate')
            nonRates <- rows - nRates
            dff[cbind(start:(start + nRates -1), start:(start + nRates - 1))] <-
                dffraw[rawsub:(rawsub + nRates - 1)]
            start <- start + nRates
            rawsub <- rawsub + nRates
            dff[start : (start + nonRates - 1 ),
                           start : (start + nonRates - 1)] <-
                dffraw[rawsub:(rawsub + nonRates * nonRates - 1)]
            start <- start + nonRates
            rawsub <- rawsub + nonRates * nonRates
        }
        dff <- dff + t(dff)
        diag(dff) <- diag(dff) / 2
        dff <- -dff
    }
    else
    {
        dff <- NULL
    }
   # browser()

   list(fra = fra, ntim0 = NULL, feasible = TRUE, OK = TRUE,
         sims=sims, dff = dff, chain = list(ans[[6]]), accepts=ans[[8]],
        rejects= ans[[9]])
}

dist2full<-function(dis) {

      n<-attr(dis,"Size")

      full<-matrix(0,n,n)

      full[lower.tri(full)]<-dis

      full+t(full)

}
