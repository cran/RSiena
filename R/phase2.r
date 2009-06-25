##/*****************************************************************************
## * SIENA: Simulation Investigation for Empirical Network Analysis
## *
## * Web: http://stat.gamma.rug.nl/siena.html
## *
## * File: phase2.r
## *
## * Description: This module contains the functions phase2.1, proc2subphase
## * and doIterations which together perform a robbins-monro stochastic
## * approximation algorithm.
## ****************************************************************************/
## args: z: internal control object
##       x: model object (readonly as not returned)
phase2.1<- function(z,x,...)
{
    #initialise phase2
    z$Phase <- 2
    z$writefreq<- 1
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
proc2subphase<- function(z, x, subphase, ...)
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
        z$prod0 <- rep(0, z$pp)
        z$prod1 <- rep(0, z$pp)
        ## ###############################################
        ## do the iterations for this repeat of this subphase
        ## ##############################################
        z <- doIterations(z, x, subphase, ...)
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
    ## browser()
    z$theta <- z$thav / (z$nit + 1)
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

doIterations<- function(z, x, subphase,...)
{
    z$nit <- 0
    ac <- 0
    repeat
    {
        z$n <- z$n+1
        z$nit <- z$nit + 1
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
          #     ## running faster with no tcl/tk
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
        zz <- x$FRAN(z, x, INIT = FALSE,...)
        if (!zz$OK)
        {
            z$OK <- zz$OK
            break
        }
        fra <- colSums(zz$fra) - z$targets
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
            if  (z$nit %% z$writefreq == 0)
            {
                z$ac <- ac
                DisplayThetaAutocor(z)
            }
        }
        ## limit change. not sure what to do here sd is not set up
        ## unless finite differences are used or ML and
        ## ML is specifically excluded here. Reporting is delayed to
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
            browser()
            fchange <- as.vector(z$gain * fra %*% z$dinv)
        }
        ##   browser()
        fchange[z$fixed] <- 0.0
        ## check positivity restriction
        z$positivized[z$nit, ] <- z$posj & (fchange > z$theta)
        fchange <- ifelse(z$posj & (fchange > z$theta), z$theta * 0.5, fchange)
        z$theta <- z$theta - fchange
        z$thav <- z$thav + z$theta
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
