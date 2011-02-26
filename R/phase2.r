##/*****************************************************************************
## * SIENA: Simulation Investigation for Empirical Network Analysis
## *
## * Web: http://www.stats.ox.ac.uk/~snidjers/siena
## *
## * File: phase2.r
## *
## * Description: This module contains the functions phase2.1, proc2subphase
## * and doIterations which together perform a robbins-monro stochastic
## * approximation algorithm.
## ****************************************************************************/
## args: z: internal control object
##       x: model object (readonly as not returned)

##@usesim siena07 Used to avoid Namespace problems with multiple processes
usesim <- function(...)
{
   simstats0c(...)
}
##@storeinFRANstore siena07 Used to avoid Namespace problems with multiple processes
storeinFRANstore <- function(...)
{
    FRANstore(...)
}
##@phase2.1 siena07 Start phase 2
phase2.1<- function(z, x, ...)
{
    #initialise phase2
    if (x$maxlike)
    {
        z$phase2fras <- array(0, dim=c(4, z$pp, 1000))
        z$rejectprops <- array(0, dim=c(4, 4, 1000))
    }
    int <- 1
    f <- FRANstore()
    z$Phase <- 2
    z$writefreq <- 1
    if (!is.batch())
    {
        tkconfigure(z$tkvars$earlyEndPhase2,state='normal')
        tkconfigure(z$tkvars$subphase,state='normal')
        tkconfigure(z$tkvars$subphaselabel,state='normal')
    }
    z$Deriv <- FALSE
    z$sd <- sqrt(apply(z$sf, 2, function(x) sum(x^2) / z$n1 - mean(x)^2))
    z$sd[z$fixed] <- 0
    ##   browser()
    Report(paste('\nPhase 2 has', x$nsub, 'subphases.\n'), cf)
    z$gain <- 2 * x$firstg
    if (x$nsub <= 0)
    {
        Report('With 0 subphases, you decided to skip phase 2!!\n')
    }
    else
    {
        if (z$maxrepeatsubphase >= 2)
        {
            Report(c('Each subphase can be repeated up to', z$maxrepeatsubphase,
                   'times\n'), cf)
        }
        z <- proc2subphase(z, x, 1, ...)
    }
    z
}
##@proc2subphase siena07 Do one subphase of phase 2
proc2subphase <- function(z, x, subphase,  useAverage=TRUE, ...)
{
    ## init subphase of phase 2
    z <- AnnouncePhase(z, x, subphase)
    z$n2min <- z$n2minimum[subphase]
    z$n2max <- z$n2maximum[subphase]
    if (x$maxlike)
    {
        z$gain <- z$gain * 0.25
    }
    else
    {
        z$gain <- z$gain * 0.5
    }
    z$repeatsubphase <- 0
    repeat
    {
        z$repeatsubphase <- z$repeatsubphase + 1
        z$truncated <- rep(FALSE, z$n2max)
        z$positivized <- matrix(FALSE, nrow = z$n2max, ncol = z$pp)
        z$ctime <- proc.time()[3]
        z$time1 <- proc.time()[3]
        z$thav <- z$theta
        z$thavn <- 1
        ## cat(z$thav, z$theta, '\n')
        z$prod0 <- rep(0, z$pp)
        z$prod1 <- rep(0, z$pp)
        ## ###############################################
        ## do the iterations for this repeat of this subphase
        ## ##############################################
        ##z <- doIterationsCopy(z, x, subphase, ...)
        z <- doIterations(z, x, subphase, ...)
        ##   if (z$nit == 50) browser()
        if (!z$OK || UserInterruptFlag() || UserRestartFlag() ||
            EarlyEndPhase2Flag())
        {
            break
        }
        ##
        ## end processing for this repeat of this subphase
        ##
        ## report truncations and positivizations
        if (any(z$truncated))
        {
            truncations<- (1 : length(z$truncated))[z$truncated]
            msg<- paste('Intervention 2.', subphase,
                        '.1: changes truncated, iterations: ',sep='')
            Report(msg, cf)
            Report(which(z$truncated), cf, fill=80)
        }
        if (any(z$positivized))
        {
            msg <- paste('Intervention 2.',subphase,
                         '.2: positivity restriction:\n ',sep='')
            Report(msg,cf)
            subs <- which(rowSums(z$positivized) > 0)
            msg<- sapply(subs, function(i, y)
                         paste('Observation:', i, 'Coordinate(s):',
                               paste((1:z$pp)[y[i,]], collapse = ' ')),
                         y = z$positivized)
            Report(msg, cf, fill=80)
        }
        if (z$maxacor >= sqrt(2.0 / (z$nit + 1)))
        {
            Report('Warning: an autocorrelation is positive at the end of',cf)
            Report(' this subphase.\n',cf)
            Report('Autocorrelations:\n',cf)
            prtmat <- z$prod1 / z$prod0
            PrtOutMat(as.matrix(prtmat), cf)
        }
        if (z$nit >= z$n2max || z$maxacor < 1e-10 ||
            z$repeatsubphase >= z$maxrepeatsubphase)
        {
            break
        }
    }
   ##finalize the subphase
    if (!z$OK || UserInterruptFlag() || UserRestartFlag())
    {
        return(z)
    }
    if (EarlyEndPhase2Flag())
    {
        Report('The user asked for early end of phase 2.\n', outf)
    }
    if (x$checktime)
    {
        time1 <- proc.time()[[3]] - z$ctime
        subphaseTime <- time1 / z$nit
        Report(paste('Time per iteration in phase ', z$Phase, '.', subphase,
                     ' = ', format(subphaseTime, nsmall=4, digits=4),
                     '\n', sep=''), lf)
    }
    if (useAverage)
    {
        z$theta <- z$thav / z$thavn #(z$nit + 1)
    }
    else
    {
        ##  use regession
        cat(z$thav / z$thavn, '\n') #(z$nit + 1)
        mylm <- lm(z$sf[1:z$nit, ] ~ z$thetaStore[1:z$nit, ])
        coefs <- coef(mylm)[-1, ]
        newvals <- solve(t(coefs), - mylm$coef[1, ])
        z$theta <- newvals
        cat(z$theta, '\n')
    }
    DisplayThetaAutocor(z)
    ##    cat('it',z$nit,'\n')
    ##recalculate autocor using -1 instead of -2 as error
    ac <- ifelse (z$prod0 > 1e-12, z$prod1 / z$prod0, -1)
    maxacor <- max(-99, ac[!z$fixed]) ##note -1 > -99
    ##   browser()
    Report(paste('Phase ', z$Phase,'.', subphase, ' ended after ', z$nit,
                 ' iterations.\n', sep = ''), cf)
    if (maxacor >= sqrt(2 / (z$nit + 1)))
    {
        Report('Warning. Autocorrelation criterion not satisfied.\n', cf)
    }
    WriteOutTheta(z)
    if (EarlyEndPhase2Flag())
    {
        return(z)
    }
    if (subphase == 2 && z$restart) ## this means we restarted in phase 1
        ##because of epsilon and need to restart the whole thing again now
    {
        Report('Restart after subphase 2.2 from current parameter values\n', cf)
        Report('because the initial values used ', cf)
        Report('led to questionable epsilon values\n', cf)
        z$fixed[z$newfixed] <- FALSE
        z$restarted <- TRUE
    }
    z
} ##end of this subphase

##@doIterations siena07 Do all iterations for 1 repeat of 1 subphase of phase 2
doIterations<- function(z, x, subphase,...)
{
    z$nit <- 0
    ac <- 0
    xsmall <- NULL
    zsmall <- makeZsmall(z)
    z$returnDeps <- FALSE
    repeat
    {
        z$n <- z$n+1
        z$nit <- z$nit + 1
        if (subphase == 1 && z$nit == 2)
            z$time1 <- proc.time()[[3]]
        if (subphase == 1 && z$nit == 11)
        {
            time1 <- proc.time()[[3]] - z$time1
            if (time1 > 1e-5)
            {
                z$writefreq <- max(1, round(20.0 / time1))
            }
            else
            {
                z$writefreq <- 20
            }
          #  if (is.batch())
          #  {
          #      z$writefreq <-  z$writefreq * 2 ##compensation for it
          #      ## running faster with no tcl/tk
          #  }
            z$writefreq <- roundfreq(z$writefreq)
        }
        if ((z$nit <= 10) || (z$nit %% z$writefreq ==0))
        {
            DisplayIteration(z)
           if (is.batch())
            {
                val <- getProgressBar(z$pb)
                increment <- ifelse(z$nit <= 10, 1, z$writefreq)
                Report(paste('Phase ', z$Phase, ' Subphase ', subphase,
                             ' Iteration ', z$nit,' Progress: ',
                             round((increment + val) /
                                   z$pb$pbmax * 100),
                             '%\n', sep = ''))
                z$pb <- setProgressBar(z$pb, val + increment)
            }
            else
            {
               if (z$nit>1)
                {
                    DisplayDeviations(z, fra)
                }
                if  (z$nit %% z$writefreq == 0)
                {
                    val <- getProgressBar(z$pb)
                    z$pb <-setProgressBar(z$pb, val + z$writefreq)
                }
            }
        }
        zsmall$nit <- z$nit
        if (z$int == 1)
        {
            zz <- x$FRAN(zsmall, xsmall)
          ##  browser()
            fra <- colSums(zz$fra) - z$targets
            if (!zz$OK)
            {
                z$OK <- zz$OK
                break
            }
        }
        else
        {
            zz <- clusterCall(z$cl, usesim, zsmall, xsmall)
            fra <- sapply(zz, function(x) colSums(x$fra)- z$targets)
            dim(fra) <- c(z$pp, z$int)
            fra <- rowMeans(fra)
            zz$OK <- sapply(zz, function(x) x$OK)
            if (!all(zz$OK))
            {
                z$OK <- FALSE
                break
            }
        }
        if (x$maxlike)
        {
         #   z$phase2fras[subphase, ,z$nit] <- fra
         #   z$rejectprops[subphase, , z$nit] <- zz$rejectprop
        }
        if (z$nit %% 2 == 1)
        {
            prev.fra <- fra
        }
        else
        {
            z$prod0 <- z$prod0 + fra * fra
            z$prod1 <- z$prod1 + fra * prev.fra
            ac <- ifelse (z$prod0 > 1e-12, z$prod1 / z$prod0, -2)
            z$maxacor <- max(-99, ac[!z$fixed]) ##note -2 > -99
            z$minacor <- min(1, ac[(!z$fixed) & ac > -1.0])
            z$ac <- ac
            if  (z$nit %% z$writefreq == 0)
            {
                DisplayThetaAutocor(z)
            }
        }
        ## limit change.  Reporting is delayed to
        ## end of phase.
     ##   browser()
        if (x$diag)## !maxlike at present
        {
            maxrat<- max(ifelse(z$sd, abs(fra)/ z$sd, 1.0))#### check this
            if (maxrat > x$maxmaxrat)
            {
                maxrat <- x$maxmaxrat / maxrat
                z$truncated[z$nit] <- TRUE
            }
            else
                maxrat <- 1.0
            fchange<- z$gain * fra * maxrat / diag(z$dfra)
        }
        else
        {
            fchange <- as.vector(z$gain * fra %*% z$dinv)
        }
        ##   browser()
        fchange[z$fixed] <- 0.0
        ## check positivity restriction
        z$positivized[z$nit, ] <- z$posj & (fchange > z$theta)
        fchange <- ifelse(z$posj & (fchange > z$theta), z$theta * 0.5, fchange)
        zsmall$theta <- zsmall$theta - fchange
        z$theta <- zsmall$theta
        z$thav <- z$thav + zsmall$theta
        z$thavn <- z$thavn + 1
        if (x$maxlike && !is.null(x$moreUpdates) && x$moreUpdates > 0)
        {
            z <- doMoreUpdates(z, x, x$moreUpdates * subphase)
            zsmall$theta <- z$theta
        }
        ##check for user interrupt
       ##   browser()
        CheckBreaks()
        if (UserInterruptFlag() || UserRestartFlag() || EarlyEndPhase2Flag())
        {
            break
        }
        ## do we stop?
        if ( (z$nit >= z$n2min && z$maxacor < 1e-10)
            || (z$nit >= z$n2max) || (z$nit >= 50 && z$minacor < -0.8 &&
                z$repeatsubphase < z$maxrepeatsubphase))
        {
            break
        }
    }
    z
}
##@doIterationsCopy siena07 Do all iterations for 1 repeat of 1 subphase of phase 2
doIterationsCopy <- function(z, x, subphase, numberIterations=10,
                             ...)
{
    ## designed to try things out!
    if (z$cconditional)
    {
        stop("cannot use this routine with conditional estimation")
    }
    z$nit <- 0
    ac <- 0
    xsmall <- NULL
    zsmall <- makeZsmall(z)
    z <- setUpPhase2Storage(z, numberIterations)
   repeat
    {
        z$n <- z$n+1
        z$nit <- z$nit + 1
        if (subphase == 1 && z$nit > 1)
            z$time1 <- proc.time()[[3]]
        if (subphase == 1 && z$nit > 10)
        {
            time1 <- proc.time()[[3]] - z$time1
            if (time1 > 1e-5)
            {
                z$writefreq <- max(1, round(20.0 / time1))
            }
            else
            {
                z$writefreq <- 20
            }
            ##  if (is.batch())
            ##  {
            ##      z$writefreq <-  z$writefreq * 2 ##compensation for it
            ##      ## running faster with no tcl/tk
            ##  }
            z$writefreq <- roundfreq(z$writefreq)
            z$writefreq <- 1
        }
        if ((z$nit <= 10) || (z$nit %% z$writefreq ==0))
        {
            DisplayIteration(z)
            if (is.batch())
            {
                val <- getProgressBar(z$pb)
                increment <- ifelse(z$nit <= 10, 1, z$writefreq)
                Report(paste('Phase ', z$Phase, ' Subphase ', subphase,
                             ' Iteration ', z$nit,' Progress: ',
                             round((increment + val) /
                                   z$pb$pbmax * 100),
                             '%\n', sep = ''))
                z$pb <- setProgressBar(z$pb, val + increment)
            }
            else
            {
                if (z$nit > 1)
                {
                    DisplayDeviations(z, fra)
                }
                if  (z$nit %% z$writefreq == 0)
                {
                    val <- getProgressBar(z$pb)
                    z$pb <-setProgressBar(z$pb, val + z$writefreq)
                }
            }
        }
        zsmall$nit <- z$nit

        zsmall$addChainToStore <- FALSE
        zsmall$needChangeContributions <- FALSE

        if (z$int == 1) ## not using a cluster
        {
            fra <- z$targets - z$targets
            for (i in 1:numberIterations)
            {
                zz <- x$FRAN(zsmall, xsmall, returnDeps=TRUE)
                ##  browser()
                fra0 <- colSums(zz$fra) - z$targets
                fra <- fra + fra0
                if (!zz$OK)
                {
                    z$OK <- zz$OK
                    break
                }
                z <- storeChainsAndFra(z, zz, i, fra0)
            }
            if (!z$OK)
            {
                break
            }
            fra <- fra / numberIterations
            z <- calculateLikelihoods(z)
        }
        else
        {
            zz <- clusterCall(z$cl, usesim, zsmall, xsmall)
            fra <- sapply(zz, function(x) colSums(x$fra) - z$targets)
            dim(fra) <- c(z$pp, z$int)
            fra <- rowMeans(fra)
            zz$OK <- sapply(zz, function(x) x$OK)
            if (!all(zz$OK))
            {
                z$OK <- FALSE
                break
            }
        }
        if (x$maxlike)
        {
            # z$phase2fras[subphase, ,z$nit] <- fra
            ##   z$rejectprops[subphase, , z$nit] <- zz$rejectprop
        }
        if (z$nit %% 2 == 1)
        {
            prev.fra <- fra
        }
        else
        {
            z$prod0 <- z$prod0 + fra * fra
            z$prod1 <- z$prod1 + fra * prev.fra
            ac <- ifelse (z$prod0 > 1e-12, z$prod1 / z$prod0, -2)
            z$maxacor <- max(-99, ac[!z$fixed]) ##note -2 > -99
            z$minacor <- min(1, ac[(!z$fixed) & ac > -1.0])
            z$ac <- ac
            if  (z$nit %% z$writefreq == 0)
            {
                DisplayThetaAutocor(z)
            }
        }
        z <- doChangeStep(z, x, fra)
        zsmall$theta <- z$theta

        if (x$maxlike && !is.null(x$moreUpdates) && x$moreUpdates > 0)
        {
            z <- doMoreUpdates(z, x, x$moreUpdates * subphase)
            zsmall$theta <- z$theta
        }
        ## importance sampling steps
        importSub <- 1
        repeat
        {
           ## browser()
            z <- predictOutcomes(z)
            #cat(z$varWeights,'\n')
            if (is.na(z$varWeights) ||
                z$varWeights > var(c(rep(0, numberIterations-1), 1))/10)
            {
                break
            }
            z <- doChangeStep(z, x, z$predictedStats)
            zsmall$theta <- z$theta
            importSub <- importSub + 1
            if (importSub > 25)
            {
                break
            }
        }
        if (zsmall$addChainToStore)
        {
            clearStoredChains()
        }
        ##check for user interrupt
        ##   browser()
        CheckBreaks()
        if (UserInterruptFlag() || UserRestartFlag() || EarlyEndPhase2Flag())
        {
            break
        }
        ## do we stop?
        if ( (z$nit >= z$n2min && z$maxacor < 1e-10)
            || (z$nit >= z$n2max) || (z$nit >= 50 && z$minacor < -0.8 &&
                z$repeatsubphase < z$maxrepeatsubphase))
        {
            break
        }
    }
    z
}

##@setUpPhase2Storage siena07 create stores and values in nonstandard Phase2
setUpPhase2Storage <- function(z, niter)
{
    nGroup <- z$f$nGroup
    z$chains <- lapply(1:nGroup, function(x)
                      lapply(1:z$groupPeriods[x], function(y)
                             vector("list", niter)))
    z$lik0 <- rep(0, niter * sum(z$groupPeriods - 1))
    z$iterFra <- matrix(0, ncol=z$pp, nrow=niter)
    z$thetaStore <- matrix(0, ncol=z$pp, nrow=z$n2max)
    z$sf <- matrix(0, ncol=z$pp, nrow=z$n2max)
    z
}

##@storeChainaAndFra siena07 update step for storage in non-standard phase 2
storeChainsAndFra <- function(z, zz, i, fra)
{
    for (ii in 1:z$nGroup)
    {
        for (jj in 1:z$groupPeriods[ii])
        {
            z$chains[[ii]][[jj]][[i]] <- zz$chain[[ii]][[jj]]
        }
    }
    z$iterFra[i, ] <- fra
    z$thetaStore[z$nit, ] <- z$theta
    z$sf[z$nit, ] <- fra
    z
}

##@calculateLikelihoods siena07 for use in non-standard phase 2
calculateLikelihoods <- function(z)
{
    storeSub <- 1
    for (ii in 1:z$nGroup)
    {
        for (jj in 1:z$groupPeriods[ii])
        {
            chains <- z$chains[[ii]][[jj]]
            for (chain in chains)
            {
                if (length(chain) == 0)
                {
                    stop("no events with this theta")
                }
                z$lik0[storeSub] <-
                    getLikelihoodPhase2(chain, z$nactors[[ii]],
                                  z$theta[z$rateParameterPosition[[ii]][[jj]]])
                storeSub <- storeSub + 1
            }
        }
    }
    z
}

##@getLikelihoodPhase2 siena07 for use in non-standard phase 2
getLikelihoodPhase2 <- function(chain, nactors, lambda)
{
    loglik <- 0
    ncvals <- sapply(chain, function(x)x[[3]])
    nc <- nactors
    nc[] <- 0
    ncvals <- table(ncvals)
    nc[names(ncvals)] <- ncvals
    logChoiceProb <- sapply(chain, function(x)x[[9]])
    logOptionSetProb <- sapply(chain, function(x)x[[8]])
    loglik <- sum(logChoiceProb) # + sum(logOptionSetProb)
    #print(sum(logOptionSetProb))
    loglik <- loglik - sum(nactors * lambda) + sum(nc * log(lambda))
    loglik
}

##@predictOutcomes siena07 for use in non-standard phase 2
predictOutcomes <- function(z)
{
    ## now the likelihood for the new theta
    ps <- getProbabilitiesPhase2(z$chain, z$theta, z$nactors,
                           z$rateParameterPosition, z$cl)$lik
    ps <- exp(unlist(ps) - z$lik0)
    ps <- ps / sum(ps)
    z$predictedStats <-  apply(z$iterFra, 2, function(x) sum(x * ps))
    z$varWeights <- var(ps)
    z
}
##@doChangeStep siena07 for use in non-standard phase 2, or outside siena07
doChangeStep <- function(z, x, fra)
{
    ## limit change.  Reporting is delayed to
    ## end of phase.
    ##   browser()
    if (x$diag)## !maxlike at present
    {
        maxrat<- max(ifelse(z$sd, abs(fra)/ z$sd, 1.0))#### check this
        if (maxrat > x$maxmaxrat)
        {
            maxrat <- x$maxmaxrat / maxrat
            z$truncated[z$nit] <- TRUE
        }
        else
            maxrat <- 1.0
        fchange<- z$gain * fra * maxrat / diag(z$dfra)
    }
    else
    {
        fchange <- as.vector(z$gain * fra %*% z$dinv)
    }
    ##   browser()
    fchange[z$fixed] <- 0.0
    ## check positivity restriction
    if (!is.null(z$nit))
    {
        z$positivized[z$nit, ] <- z$posj & (fchange > z$theta)
    }
    fchange <- ifelse(z$posj & (fchange > z$theta), z$theta * 0.5, fchange)

    z$theta <- z$theta - fchange
    if (!is.null(z$thav))
    {
      z$thav <- z$thav + z$theta
    z$thavn <- z$thavn + 1
    }
    z
}

##@clearStoredChains algorithms Clears storage used for chains in C.
clearStoredChains <- function()
{
    f <- FRANstore()
    .Call("clearStoredChains", PACKAGE=pkgname, f$pModel)
}

##@getProbabilitiesPhase2 algorithms Recalculates change contributions
getProbabilitiesPhase2 <- function(chain, theta, nactors, rateParameterPosition,
                              cl=NULL)
{
    f <-  FRANstore()#
    getScores <- FALSE
    ##print(theta)
    for (i in 1:length(chain)) # group
    {
        for (j in 1:length(chain[[i]])) # period
        {
            if (!is.null(cl) )
            {
                tmp <- parLapply(cl, chain[[i]][[j]],
                                 function(x, i, j, theta, getScores, k, n)
                             {
                                 f <- FRANstore()
                                 resp <- .Call("getChainProbabilitiesList",
                                               PACKAGE = pkgname, x,
                                               f$pData, f$pModel, as.integer(i),
                                               as.integer(j), f$myeffects,
                                               theta, getScores)
                                 lik <-
                                     getLikelihoodPhase2(resp[[1]], nactors[[i]],
                                                   theta[k])
                                 if (getScores)
                                 {
                                     sc <- resp[[2]]
                                 }
                                 else
                                 {
                                     sc <- NULL
                                 }
                                 list(lik=lik, sc=sc)
                             }, i=i, j=j, theta=theta, getScores=getScores,
                                 k=rateParameterPosition[[i]][[j]],
                                 n=nactors[[i]]
                                 )
                lik <- sapply(tmp, function(x)x[[1]])
                if (getScores)
                {
                    sc <- t(sapply(tmp, function(x)x[[2]]))
                }
                else
                {
                    sc <- NULL
                }
            }
            else
            {
                lik <- rep(0, length(chain[[i]][[j]]))
                sc <- matrix(0, nrow=length(chain[[i]][[j]]),
                             ncol=length(theta))
                for (k in 1:length(chain[[i]][[j]])) # one chain
                {
                    resp <- .Call("getChainProbabilitiesList",
                                  PACKAGE = pkgname,
                                  chain[[i]][[j]][[k]],
                                  f$pData, f$pModel, as.integer(i),
                                  as.integer(j), f$myeffects, theta,
                                  getScores)

                    lik[k] <-
                        getLikelihoodPhase2(resp[[1]], nactors[[i]],
                                      theta[rateParameterPosition[[i]][[j]]])
                    if (getScores)
                    {
                        sc[k,] <- resp[[2]]
                    }
                }
            }
        }
    }
    list(lik=lik, sc=sc)
}
##@initForAlgorithms siena07 stores values for use in nonstandard Phase2
initForAlgorithms <- function(z)
{
    if (z$cconditional)
    {
        return(z)
    }
    nGroup <- z$f$nGroup
    z$nDependentVariables <- length(z$f$depNames)
    atts <- attributes(z$f)
    z$groupPeriods <- atts$groupPeriods - 1
    netnames <- names(z$f[[1]]$depvars)
    z$nactors <-  lapply(1:nGroup, function(i, periods, data)
                   {
                       tmp <- sapply(data[[i]]$depvars, function(x)
                                          dim(x)[1])
                      tmp <- tmp[match(netnames, names(data[[i]]$depvars))]
                      tmp
                   }, periods=z$groupPeriods, data=z$f[1:nGroup]
                       )
    z$rateParameterPosition <-
        lapply(1:nGroup, function(i, periods, data)
           {
               lapply(1:periods[i], function(j)
                  {
                      rateEffects <-
                          z$effects[z$effects$basicRate &
                                    z$effects$period == j &
                                    z$effects$group == i,]
                      rateEffects <-
                          rateEffects[match(netnames,
                                            rateEffects$name), ]
                      tmp <- as.numeric(row.names(rateEffects))
                      names(tmp) <- netnames
                      tmp
                  }
                      )
           }, periods=z$groupPeriods, data=z$f[1:nGroup]
               )
    z$evalParameterPosition <-
        lapply(netnames, function(x)
           {
               as.numeric(row.names(z$effects[z$effects$name == x &
                                              z$effect$type =="eval", ]))
           }
               )
    names(z$evalParameterPosition) <- netnames
    z$endowParameterPosition <-
        lapply(netnames, function(x)
           {
               as.numeric(row.names(z$effects[z$effects$name == x &
                                              z$effect$type =="endow", ]))
           }
               )
    z$basicRate <- z$effects$basicRate
    z$nGroup <- nGroup
    z
}
