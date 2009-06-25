#/******************************************************************************
# * SIENA: Simulation Investigation for Empirical Network Analysis
# *
# * Web: http://www.stats.ox.ac.uk/~snidjers/siena
# *
# * File: robmon.r
# *
# * Description: This module contains the function robmon which controls
# * the phases of the robbins-munro stochastic approximation algorithm.
# *****************************************************************************/
##args:x: model object - intended to be read only
##     z: model fitting object
## returns updated z
robmon <- function(z, x, useCluster, noClusters, initC, ...)
{
    z$FinDiff.method<- x$FinDiff.method
    z$n <- 0
    z$OK <-  TRUE
    z$error<- FALSE
    z$restarted <- FALSE
    z$DerivativeProblem <- FALSE
    z$ForceFinDifPhase1 <- FALSE
    z$Phase3Interrupt <- FALSE
    z$repeatsForEpsilon <- 0
    z$maxrepeatsubphase <- 4
    z$gain <- x$firstg
    z$haveDfra <- FALSE
    #######################################################
    ##do initial setup call of FRAN
    #######################################################
    z <- x$FRAN(z, x, INIT=TRUE, ...)
    ##
    ##if conditional, FRAN changes z$theta etc
    #######################################################
    z$newFixed <- rep(FALSE,z$pp)
    z$AllNowFixed <- FALSE
    z$dinv <- matrix(NA, nrow = z$pp, ncol = z$pp)
    z$scale <- rep(0.1, z$pp)
    Report('\n', outf)
    Report('\nStochastic approximation algorithm.\n', cf)
    if (x$firstg<=0)
    {
        Report(c('Initial value of the gain parameter is ', x$firstg,
                 '.\n'), outf)
        Report('This is not allowed; changed to 0.0001.\n', outf)
        z$gain <- 0.0001
    }
    Report(c('Initial value for gain parameter = ', format(z$gain),
             '.\nStart of the algorithm.\n'), cf, sep='')
    Report('Target function values are \n', cf)
    ftargets <- format(z$targets, width = 10, nsmall = 4)
    fnum<- format(1 : z$pp, width = 3)
    Report(c(paste(fnum, '. ', ftargets, sep = '')), cf, fill=80)
    z$epsilon<- pmin(0.1,z$scale)
    z$epsilon[z$posj]<- 0.1 * z$theta[z$posj]
    z$theta0<- z$theta ## store starting value without any conditioning variables
    z$anyposj <- any(z$posj)
    z$resist <- rep(1, z$pp)
    z$n1 <- 7 + 3 * z$pp
    if (any(!z$fixed))
    {
        z$AllUserFixed<- FALSE
    }
    else
    {
        z$AllUserFixed <- TRUE
    }
    ##browser()
    repeat  ##this is startagain:
    {
      z$epsilonProblem <- FALSE
      repeat ## this loop is simply to break out of! only intend to do it once
        {
            if (any(!z$fixed))
                z$AllNowFixed<- FALSE
            else
                z$AllNowFixed <- TRUE
            phase3Only <- FALSE
            if (!is.batch())
            {
                tkdelete(z$tkvars$subphase, 0, "end")
                tkconfigure(z$tkvars$earlyEndPhase2, state="disabled")
                tkconfigure(z$tkvars$subphase, state="disabled")
                tkconfigure(z$tkvars$subphaselabel, state="disabled")
            }
            if (z$AllNowFixed || x$nsub == 0)
            {
                if (z$AllNowFixed)
                {
                    Report('All parameters are fixed.\n', outf)
                    Report(c('All parameters fixed; no estimation;',
                             'only Phase 3.\n'), lf)
                }
                else
                {
                    Report('Number of subphases is specified as 0.\n', outf)
                    Report('0 subphases; no estimationl only phase 3.\n', lf)
                }
                Report(c('Therefore the estimation phase is skipped\n',
                         'and the program passes on immediately to phase 3\n',
                         'for checking the current parameter values and',
                         'calculating standard errors.\n'),
                       outf)
                phase3Only <- TRUE
            }
            NullChecks() ## reset user interrupt variables and flags
            z$Phase <- 0
            z <- AnnouncePhase(z, x)
            WriteOutTheta(z)
            z$fixed[z$newFixed] <- FALSE
            z$newFixed <- rep(FALSE,z$pp) ##not clear what we should do here
            z$restart <- FALSE
            z$OK <- TRUE
            if (!phase3Only)
            {
                if (!z$haveDfra)
                {
                    ##start phase1 and do 10 iterations,
                    z <- phase1.1(z, x, ...)
                    ##check epsilon
                    if (!z$OK || UserInterruptFlag() || UserRestartFlag())
                    {
                        break
                    }
                    if (z$epsilonProblem && z$repeatsForEpsilon<4)
                    {
                        z$repeatsForEpsilon<- z$repeatsForEpsilon+1
                        z$restart<- !z$restarted
                        break
                    }
                    z<- phase1.2(z,x,...)
                    ## browser()
                    if (!z$OK || z$DerivativeProblem ||
                        UserInterruptFlag() || UserRestartFlag())
                    {
                        break
                    }
                }
                if (x$nsub > 0)
                {
                    z <- phase2.1(z, x, ...)
                }
                if (!z$OK || UserInterruptFlag() || UserRestartFlag())
                {
                    if (!is.batch())
                    {
                        tkdelete(z$tkvars$subphase,0,'end')
                        tkconfigure(z$tkvars$earlyEndPhase2,state='disabled')
                        tkconfigure(z$tkvars$subphase,state='disabled')
                        tkconfigure(z$tkvars$subphaselabel,state='disabled')
                    }
                    break
                }
                if (x$nsub>1 && !EarlyEndPhase2Flag())
                {
                    z <- proc2subphase(z,x,2,...)
                }
                if (!z$OK || z$restart||
                    UserInterruptFlag() || UserRestartFlag() )
                {
                    if (!is.batch())
                    {
                        tkdelete(z$tkvars$subphase,0,'end')
                        tkconfigure(z$tkvars$earlyEndPhase2,state='disabled')
                        tkconfigure(z$tkvars$subphase,state='disabled')
                        tkconfigure(z$tkvars$subphaselabel,state='disabled')
                    }
                    break
                }
                if (x$nsub > 2 && !EarlyEndPhase2Flag())
                {
                    for (i in 3 : x$nsub)
                    {
                        z <- proc2subphase(z, x, i, ...)
                        if (!z$OK || UserInterruptFlag() ||
                            UserRestartFlag() || EarlyEndPhase2Flag())
                        {
                            if (!is.batch())
                            {
                                tkdelete(z$tkvars$subphase,0,'end')
                                tkconfigure(z$tkvars$earlyEndPhase2,state='disabled')
                                tkconfigure(z$tkvars$subphase,state='disabled')
                                tkconfigure(z$tkvars$subphaselabel,state='disabled')
                            }
                            break
                        }
                    }
                    if (!z$OK || UserInterruptFlag() ||
                        UserRestartFlag() )
                    {
                        if (!is.batch())
                        {
                            tkdelete(z$tkvars$subphase,0,'end')
                            tkconfigure(z$tkvars$subphase,state='disabled')
                            tkconfigure(z$tkvars$subphaselabel,state='disabled')
                            tkconfigure(z$tkvars$earlyEndPhase2,state='disabled')
                        }
                        break
                    }
                }
            }
            if (!is.batch())
            {
                tkdelete(z$tkvars$subphase,0,'end')
                tkconfigure(z$tkvars$subphase,state='disabled')
                tkconfigure(z$tkvars$subphaselabel,state='disabled')
                tkconfigure(z$tkvars$earlyEndPhase2,state='disabled')
            }
            z<- phase3(z, x, useCluster, noClusters, initC, ...)
            break
        }
        ##stop if not OK or user has asked to
        if (!z$OK || UserInterruptFlag())
        {
            break
        }
        ##break unless we want to start at phase 1 again.
        if (!z$DerivativeProblem && !z$restart && !UserRestartFlag())
        {
            break
        }
    }
    if (!z$OK || UserInterruptFlag())
    {
        if (z$OK)
        {
            if (!z$Phase3Interrupt)
                z$termination <- 'UserInterrupt'
            else
                z$termination <- 'OK'
        }
        else
        {
            z$termination <- 'Error'
        }
        if (!z$OK || !z$Phase3Interrupt)
            return(z)
    }
    ## #####################################################
    ## do final call of FRAN
    ## #####################################################
    z <- x$FRAN(z, x, TERM=TRUE,...)
    ## #####################################################
    ## call to FRAN changes covariance matrix for conditional estimation
    z$diver<- (z$fixed | z$diver | diag(z$covtheta) < 1e-9) & (!z$AllUserFixed)
    z$covtheta[z$diver, ] <- Root(diag(z$covtheta)) * 33
    ##not sure this does not use very small vals
    z$covtheta[, z$diver] <- Root(diag(z$covtheta)) * 33
    diag(z$covtheta)[z$diver] <- 999
    z$termination <- 'OK'
    z
}

