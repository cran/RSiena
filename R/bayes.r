##/*****************************************************************************
## * SIENA: Simulation Investigation for Empirical Network Analysis
## *
## * Web: http://www.stats.ox.ac.uk/~snidjers/siena
## *
## * File: bayes.r
## *
## * Description: This file contains the code to run Bayesian simulation
## *
## ****************************************************************************/
##@bayes algorithm  fit a Bayesian model
bayes <- function(data, effects, model, nwarm=100, nmain=100, nrunMHBatches=20,
                 nrunMH=100, save=TRUE)
{
    createStores <- function()
    {
        z$posteriorTot <- rep(0, npar)
        z$posteriorMII <- matrix(0, nrow=npar, ncol=npar)
        z$lambdadist <- matrix(NA, nrow=nmain * nrunMHBatches,
                               ncol=sum(basicRate))
        z$lambdas  <- matrix(NA, nrow=nmain * nrunMHBatches,
                             ncol=sum(basicRate))
        z$betas <- matrix(NA, nrow=nmain * nrunMHBatches, ncol=sum(!basicRate))
        z$candidates <- matrix(NA, nrow=nmain * nrunMHBatches,
                               ncol=sum(!basicRate))
        z$acceptances <- rep(NA, nmain * nrunMHBatches)
        z$MHacceptances <- matrix(NA, nrow=nmain * nrunMHBatches, ncol=6)
        z$MHrejections <- matrix(NA, nrow=nmain * nrunMHBatches , ncol=6)
        z$MHproportions <- matrix(NA, nrow=nmain *  nrunMHBatches, ncol=6)
        z
    }
    storeData <- function()
    {
        start <- z$sub + 1
        nrun <- nrow(z$parameters)
        end <- start + nrun - 1
        z$accepts <- as.logical(z$accepts)
        z$acceptances[start:end] <- z$accepts
        z$lambdas[start:end, ] <- z$parameters[, basicRate]
        z$lambdadist[start:end, ] <- z$shapes[, basicRate]
        z$candidates[start:end, ] <- z$parameters[, !basicRate]
        z$parameters[!z$accepts, !basicRate] <- NA
        z$parameters[, !basicRate] <-
            carryForward(z$parameters[, !basicRate], z$betas[1:end,])
        z$betas[start:end, ] <- z$parameters[, !basicRate]

        z$posteriorTot <- z$posteriorTot + colSums(z$parameters)
        for (i in npar)
        {
            z$posteriorMII <- z$posteriorMII +
                outer(z$parameters[i, ], z$parameters[i, ])
        }
        z$MHacceptances[start:end, ] <- z$MHaccepts
        z$MHrejections[start:end, ] <- z$MHrejects
        z$MHproportions[start:end, ] <-z$MHaccepts/ (z$MHaccepts + z$MHrejects)
        z$sub <- z$sub + nrun
        z
    }

    carryForward <- function(parameters, betas)
    {
        npar <- nrow(parameters)
        nbeta <- nrow(betas)
        if (npar < nbeta)
        {
            parameters <- rbind(betas[(nbeta - npar), ], parameters)
        }
        parameters <-
            apply(parameters, 2, function(x)
              {
                  x.pos <- which(!is.na(x))
                  if (length(x.pos) == 0 || x.pos[1] != 1)
                  {
                      x.pos <- c(1, x.pos)
                  }
                  x[rep(x.pos, c(diff(x.pos),
                                 length(x) - x.pos[length(x.pos)] + 1))]
              }
                  )
        if (npar < nbeta)
        {
            parameters[-1, ]
        }
        else
        {
            parameters
        }
    }

    improveMH <- function(z, x, tiny=1.0e-15, desired=40, maxiter=100,
                          tolerance=15)
    {
        rescaleCGD <- function(iter)
        {
            if (actual > desired)
            {
                u <-  2 - ((iter - actual) / (iter - desired))
            }
            else
            {
                u <- 1 / (2 - (actual / desired))
            }
            if (tol <- abs(actual - desired) <= tolerance)
            {
                number <<- number + 1
                if (number == 2) success <<- TRUE
            }
            else
            {
                number <<- 0
            }
            u
        }
        iter <- 0
        number <- 0
        success <- FALSE
        repeat
        {
            iter <- iter + 1
            z <- MCMCcycle(z, nrunMH=1, nrunMHBatches=100)
            actual <- z$BayesAcceptances ## acceptances
            ans <- rescaleCGD(100)
            z$scaleFactor <- z$scaleFactor * ans
            if (success | iter == maxiter)
            {
                break
            }
            if (z$scaleFactor < tiny)
            {
                cat('scalefactor < tiny\n')
                browser()
            }
        }
        cat('fine tuning took ', iter, ' iterations. Scalefactor:',
            z$scaleFactor, '\n')
        z
    }

    ## initialise
    Report(openfiles=TRUE, type="n") #initialise with no file
    if (!save)
    {
        require(lattice)
        dev.new()
        thetaplot = dev.cur()
        dev.new()
        acceptsplot = dev.cur()
    }
    z  <-  NULL
    z$FinDiff.method <- FALSE
    z$maxlike <- TRUE
    model$maxlike <- TRUE
    model$FRANname <- "maxlikec"
    z$int <- 1
    z$int2 <- 1
    model$cconditional <-  FALSE
    if (!is.null(model$randomSeed))
    {
        set.seed(model$randomSeed)
    }
    model$FRAN <- getFromNamespace(model$FRANname, pos=grep("RSiena",
                                                   search())[1])
    z <- model$FRAN(z, model, INIT=TRUE, data=data, effects=effects,
                    returnDeps=TRUE)
    ##if (useCluster)
    ##   cl <- makeCluster(rep("localhost", 2), type = "SOCK")

    is.batch(TRUE)

    WriteOutTheta(z)

    npar <- length(z$theta)
    basicRate <- z$effects$basicRate
    iter <- 0
    z$numm <- 20
    z$scaleFactor <- 1
    z <- createStores()
    z$sub <- 0
    z <- improveMH(z)

    for (ii in 1:nwarm)
    {
        z <- MCMCcycle(z, nrunMH=4, nrunMHBatches=20)
        numm <- z$numm
    }

    for (ii in 1:nmain)
    {
        z <- MCMCcycle(z, nrunMH=100, nrunMHBatches=20)
        z <- storeData()

        numm <- z$numm
        if (ii %% 10 == 0 && !save) ## do some plots
        {
            cat('main after ii',ii,numm, '\n')
            dev.set(thetaplot)
            thetadf <- data.frame(z$lambdas, z$betas)
            acceptsdf <- data.frame(z$MHproportions, z$acceptances)
            lambdaNames <- paste(z$effects$name[basicRate],
                                 z$effects$shortName[basicRate],
                                 z$effects$period[basicRate],
                                 z$effects$group[basicRate], sep=".")
            betaNames <- paste(z$effects$name[!basicRate],
                               z$effects$shortName[!basicRate], sep=".")
            names(thetadf) <- c(lambdaNames, betaNames)
            names(acceptsdf) <- c("InsDiag", "CancDiag", "Permute", "InsPerm",
                                  "CancPerm", "Missing", "BayesAccepts")
            varnames <- paste(names(thetadf), sep="", collapse= " + ")
            varcall <- paste("~ ", varnames,  sep="", collapse="")
            print(histogram(as.formula(varcall), data=thetadf, scales="free",
                            outer=TRUE, breaks=NULL))
            dev.set(acceptsplot)
            varnames <- paste(names(acceptsdf), sep="", collapse= " + ")
            varcall <- paste("~ ", varnames,  sep="", collapse="")
            print(histogram(as.formula(varcall), data=acceptsdf,
                            scales=list(x="same", y="free"),
                            outer=TRUE, breaks=NULL, type="density"))
        }
    }
    z
}

MCMCcycle <- function(z, nrunMH, nrunMHBatches)
{
    ## this function assumes only one period at the moment, but would deal with
    ## multiple periods locally, calling MCMCcycle for each. Not yet available
    ## in C, due to need to keep multiple chains in existence. I have not
    ## decided how to do it yet!
    period <- 1
    group <- 1
    f <- FRANstore()
    ans <- .Call("MCMCcycle", PACKAGE=pkgname, f$pData, f$pModel,
                 f$pMLSimulation, f$myeffects, as.integer(period),
                 as.integer(group),
                 z$scaleFactor, nrunMH, nrunMHBatches)
    ## process the return values
    z$BayesAcceptances <- sum(ans[[1]])
    z$accepts <- ans[[1]]
    z$parameters <- ans[[2]]
    z$shapes <- ans[[3]]
    z$MHaccepts <- ans[[4]]
    z$MHrejects <- ans[[5]]
    z
}



