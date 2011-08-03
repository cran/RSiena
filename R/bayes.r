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
##@bayes Bayesian fit a Bayesian model
bayes <- function(data, effects, model, nwarm=100, nmain=100, nrunMHBatches=20,
                  nrunMH=100, save=TRUE, nbrNodes=1, seed=1, dfra=NULL,
                  priorSigma=NULL, prevAns=NULL, getDocumentation=FALSE)
{
    ##@createStores internal bayes Bayesian set up stores
    createStores <- function(z)
    {
        npar <- length(z$theta)
        numberRows <- nmain * nrunMHBatches
        z$posteriorTot <- matrix(0, nrow=z$nGroup, ncol=npar)
        z$posteriorMII <- array(0, dim=c(z$nGroup, npar, npar))
        z$candidates <- array(NA, dim=c(numberRows, z$nGroup, npar))
        z$acceptances <- matrix(NA, nrow=z$nGroup, ncol=numberRows)
        z$MHacceptances <- array(NA, dim=c(numberRows, z$nGroup, 7))
        z$MHrejections <- array(NA, dim=c(numberRows, z$nGroup, 7))
        z$MHproportions <- array(NA, dim=c(numberRows, z$nGroup, 7))
        z
    }
    ##@storeData internal bayes Bayesian put data in stores
    storeData <- function()
    {
        start <- z$sub + 1
        nrun <- nrow(z$parameters)
        npar <- length(z$theta)
        end <- start + nrun - 1
        z$acceptances[, start:end] <- z$accepts
        z$candidates[start:end,, ] <- z$parameters
        z$posteriorTot <- z$posteriorTot + colSums(z$parameters)
        for (group in 1:z$nGroup)
        {
            for (i in npar)
            {
                z$posteriorMII[group, , ] <- z$posteriorMII[group, ,] +
                    outer(z$parameters[i, group, ], z$parameters[i, group,])
            }
        }
        z$MHacceptances[start:end, , ] <- z$MHaccepts
        z$MHrejections[start:end, , ] <- z$MHrejects
        z$MHproportions[start:end, , ] <- z$MHaccepts /
            (z$MHaccepts + z$MHrejects)
        z$sub <- z$sub + nrun
        z
    }
	##@improveMH internal bayes Bayesian find scale factors
    improveMH <- function(z, x, tiny=1.0e-15, desired=40, maxiter=100,
                          tolerance=15, getDocumentation=FALSE)
    {
		##@rescaleCGD internal improveMH Bayesian
        rescaleCGD <- function(iter)
        {
			u <- ifelse (actual > desired,
                         2 - ((iter - actual) / (iter - desired)),
                         1 / (2 - (actual / desired)))
            number <<- ifelse(abs(actual - desired) <= tolerance,
                               number + 1, 0 )
            success <<- number >= 2
            u
        }
		if (getDocumentation)
		{
			tt <- getInternals()
			return(tt)
		}
        iter <- 0
        number <- rep(0, z$nGroup)
        success <- rep(FALSE, z$nGroup)
        repeat
        {
            iter <- iter + 1
            z <- MCMCcycle(z, nrunMH=1, nrunMHBatches=100, change=FALSE)
            actual <- z$BayesAcceptances ## acceptances
            ans <- rescaleCGD(100)
            update <- number < 3
            z$scaleFactors[update] <- z$scaleFactors[update] * ans[update]
            cat(actual, ans, z$scaleFactors, '\n')
            if (all(success) | iter == maxiter)
            {
                break
            }
            if (any(z$scaleFactors < tiny))
            {
                cat('scalefactor < tiny\n')
                browser()
            }
        }
        cat('fine tuning took ', iter, ' iterations. Scalefactor:',
            z$scaleFactors, '\n')
        z
    }
    ## ################################
    ## start of function proper
    ## ################################
	if (getDocumentation != FALSE)
	{
		if (getDocumentation == TRUE)
		{
			tt <- getInternals()
			return(tt)
		}
		else ## need to run getInternals on the argument value
		{
			targs <- formals(getDocumentation[1])
			targs[1:length(targs)] <- 1
			targs['getDocumentation'] <- TRUE
			if (length(getDocumentation) > 1)
			{
				targs['getDocumentation'] <- getDocumentation[-1]
			}
			return(do.call(getDocumentation[1], targs))
		}
	}

    z <- initializeBayes(data, effects, model, nbrNodes, seed, priorSigma,
                         prevAns=prevAns)
    z <- createStores(z)

    z$sub <- 0

    if (is.null(z$dfra) && is.null(dfra))
    {
        z <- getDFRA(z, 10)
    }
    else
    {
        if (!is.null(dfra))
        {
            z$dfra <- dfra
        }
        lambda <- z$theta[z$basicRate]
        dfra <- z$dfra
        z$dfra[z$basicRate, ] <- z$dfra[z$basicRate,] * lambda
        z$dfra[, z$basicRate] <- z$dfra[, z$basicRate] * lambda
        diag(z$dfra)[z$basicRate] <- lambda *
            diag(dfra)[z$basicRate] + lambda * lambda *
                colMeans(z$sf)[z$basicRate]
        print(z$dfra)
        z$dfra <- solve(z$dfra)
        print(z$dfra)


    }
    z <- improveMH(z)

    if (!save)
    {
        require(lattice)
        dev.new()
        thetaplot = dev.cur()
        dev.new()
        ratesplot = dev.cur()
        dev.new()
        tseriesplot = dev.cur()
        dev.new()
        tseriesratesplot = dev.cur()
   }

    for (ii in 1:nwarm)
    {
        z <- MCMCcycle(z, nrunMH=4, nrunMHBatches=20)
    }

    for (ii in 1:nmain)
    {
        z <- MCMCcycle(z, nrunMH=nrunMH, nrunMHBatches=nrunMHBatches)
        z <- storeData()

        if (ii %% 10 == 0 && !save) ## do some plots
        {
            cat('main after ii', ii, '\n')
            dev.set(thetaplot)
            thetadf <-
                lapply(1:z$nGroup, function(i)
                   {
                       data.frame(Group=rep(i, ii * nrunMHBatches),
                                  z$candidates[1:(ii * nrunMHBatches), i, ])
                   }
                       )
            thetadf <- do.call(rbind, thetadf)
            basicRate <- z$basicRate
            ##thetadf <- data.frame(z$candidates)
            acceptsdf <- data.frame(z$MHproportions,
                                    z$acceptances)
            ratesdf <- thetadf[, -1, drop=FALSE][, z$basicRate, drop=FALSE]
            thetadf <- cbind(Group=thetadf[, 1, drop=FALSE],
							 thetadf[, -1, drop=FALSE][,
										   !z$basicRate, drop=FALSE])
            thetaNames<- paste(z$effects$name[!z$basicRate],
                               z$effects$shortName[!z$basicRate], sep=".")
            rateNames <- paste(z$effects$name[basicRate],
                                           z$effects$shortName[basicRate],
                                           z$effects$period[basicRate],
                                           z$effects$group[basicRate], sep=".")
            names(ratesdf) <- rateNames
            ratesdf <- cbind(Group=thetadf[, 1, drop=FALSE], ratesdf)
            names(thetadf)[-1] <- make.names(thetaNames, unique=TRUE)
            names(acceptsdf) <- c("InsDiag", "CancDiag", "Permute", "InsPerm",
                                  "DelPerm", "InsMissing", "DelMissing",
                                  "BayesAccepts")
            varnames <- paste(names(thetadf)[-1], sep="", collapse= " + ")
            varcall <- paste("~ ", varnames,  " | Group", sep="", collapse="")
            print(histogram(as.formula(varcall), data=thetadf, scales="free",
                            outer=TRUE, breaks=NULL, type="density",
                            panel=function(x, ...)
                        {
                            panel.histogram(x, ...)
                            panel.densityplot(x, darg=list(na.rm=TRUE), ...)
                        }
                            ))
            dev.set(ratesplot)
            varnames <- paste(names(ratesdf)[-1], sep="", collapse= " + ")
            varcall <- paste("~ ", varnames, sep="", collapse="")
            print(histogram(as.formula(varcall), data=ratesdf, scales="free",
                            outer=TRUE, breaks=NULL, type="density",
                            panel=function(x, ...)
                        {
                            panel.histogram(x, ...)
                            panel.densityplot(x, darg=list(na.rm=TRUE), ...)
                        }
                            ))
            varnames <- paste(names(thetadf)[-1], sep="", collapse= " + ")
            varcall <- paste(varnames,  "~ 1:", ii * nrunMHBatches * z$nGroup,
                             " | Group", sep="", collapse="")
            dev.set(tseriesplot)
            print(xyplot(as.formula(varcall), data=thetadf, scales="free",
                         outer=TRUE))
            varnames <- paste(names(ratesdf)[-1], sep="", collapse= " + ")
            varcall <- paste(varnames,  "~ 1:", ii * nrunMHBatches * z$nGroup,
                             sep="", collapse="")
            dev.set(tseriesratesplot)
            print(xyplot(as.formula(varcall), data=ratesdf, scales="free",
                         outer=TRUE))
            ## dev.set(acceptsplot)
            ## varnames <- paste(names(acceptsdf), sep="", collapse= " + ")
            ## varcall <- paste("~ ", varnames,  sep="", collapse="")
            ## print(histogram(as.formula(varcall), data=acceptsdf,
            ##                 scales=list(x="same", y="free"),
            ##                 outer=TRUE, breaks=NULL, type="density",
            ##                 panel=function(x, ...)
            ##             {
            ##                 panel.histogram(x, ...)
            ##                 panel.densityplot(x, darg=list(na.rm=TRUE), ...)
            ##             }))
        }
    }
    z$FRAN <- NULL
    z
}
##@MCMCcycle algorithms do some loops of (MH steps and sample parameters)
MCMCcycle <- function(z, nrunMH, nrunMHBatches, change=TRUE)
{
    z$accepts <- matrix(NA, nrow=z$nGroup, nrunMHBatches)
    z$parameters <- array(NA, dim=c(nrunMHBatches, z$nGroup, z$pp))
    z$MHaccepts <- array(NA, dim=c(nrunMHBatches, z$nGroup, 7))
    z$MHrejects <- array(NA, dim=c(nrunMHBatches, z$nGroup, 7))
    z$MHaborts <- array(NA, dim=c(nrunMHBatches, z$nGroup, 7))
    storeNrunMH <- z$nrunMH
    z$nrunMH <- nrunMH
    for (i in 1:nrunMHBatches)
    {
        ans <- z$FRAN(z, returnChains=TRUE, byGroup=TRUE)
        z$chain <- ans$chain
        z <- sampleParameters(z, change)
        z$accepts[, i] <- z$accept
        z$parameters[i, , ] <- z$thetaMat
        z$MHaccepts[i, ,] <-
            t(do.call(cbind,
                      tapply(ans$accepts, factor(z$callGrid[, 1]),
                             function(x)Reduce("+", x))))
        z$MHrejects[i, , ] <-
            t(do.call(cbind, tapply(ans$rejects, factor(z$callGrid[, 1]),
                                    function(x)Reduce("+", x))))
        z$MHaborts[i, , ] <- t(do.call(cbind,
                                       tapply(ans$aborts,
                                              factor(z$callGrid[, 1]),
                                              function(x)Reduce("+", x))))
    }
    z$BayesAcceptances <- rowSums(z$accepts)
    z$nrunMH <- storeNrunMH
    z
}
##@sampleParameters algorithms propose new parameters and accept them or not
sampleParameters <- function(z, change=TRUE)
{
    ## get a multivariate normal with covariance matrix dfra multiplied by a
    ## scale factor which varies between groups
    require(MASS)
    thetaChanges <- t(sapply(1:z$nGroup, function(i)
                         {
                             tmp <- z$thetaMat[i, ]
                             use <- !is.na(z$thetaMat[i, ])
                             tmp[use] <-
                                 mvrnorm(1, mu=rep(0, sum(use)),
                                         Sigma=z$scaleFactors[i] *
                                         z$dfra[use, use])
                             tmp
                         }
                             ))

    thetaOld <- z$thetaMat
    thetaOld[, z$basicRate] <- log(thetaOld[, z$basicRate])
    thetaNew <- thetaOld + thetaChanges

    priorOld <- sapply(1:z$nGroup, function(i)
                   {
                       tmp <- thetaOld[i, ]
                       use <- !is.na(tmp)
                       dmvnorm(tmp[use],  mean=rep(0, sum(use)),
                               sigma=z$priorSigma[use, use])
                   }
                       )
    priorNew <- sapply(1:z$nGroup, function(i)
                   {
                       tmp <- thetaNew[i, ]
                       use <- !is.na(tmp)
                       dmvnorm(tmp[use],  mean=rep(0, sum(use)),
                               sigma=z$priorSigma[use, use])
                   }
                       )
    logpOld <- getProbs(z, z$chain)

    storeNrunMH <- z$nrunMH
    z$nrunMH <- 0
    thetaNew[, z$basicRate] <- exp(thetaNew[, z$basicRate])
    z$thetaMat <- thetaNew
    resp <- z$FRAN(z, returnChains=TRUE, byGroup=TRUE)
    logpNew <- getProbs(z, resp$chain)

    proposalProbability <- priorNew - priorOld + logpNew - logpOld

    z$accept <- log(runif(length(proposalProbability))) < proposalProbability
    thetaOld[, z$basicRate] <- exp(thetaOld[, z$basicRate])
    if (!change)
    {
        z$thetaMat <- thetaOld
    }
    else
    {
        ##print(z$thetaMat)
        z$thetaMat[!z$accept, ] <- thetaOld[!z$accept, ]
    }
    z$nrunMH <- storeNrunMH
    ## cat(thetaNew, priorNew, priorOld, logpNew, logpOld,
    ##    exp(proposalProbability),
    ##    z$accept,'\n')
    z
}

##@initializeBayes algorithms do set up for Bayesian model
initializeBayes <- function(data, effects, model, nbrNodes, seed, priorSigma,
                            prevAns)
{
    ## initialise
    set.seed(seed)
    Report(openfiles=TRUE, type="n") #initialise with no file
    z  <-  NULL
    z$Phase <- 1
    z$Deriv <- FALSE
    z$FinDiff.method <- FALSE
    z$maxlike <- TRUE
    model$maxlike <- TRUE
    model$FRANname <- "maxlikec"
    z$print <- FALSE
    z$int <- 1
    z$int2 <- nbrNodes
    model$cconditional <-  FALSE
    if (!is.null(model$randomSeed))
    {
        set.seed(model$randomSeed)
    }
    z$FRAN <- getFromNamespace(model$FRANname, pos=grep("RSiena",
                                               search())[1])
    z <- z$FRAN(z, model, INIT=TRUE, data=data, effects=effects,
                prevAns=prevAns)
    is.batch(TRUE)

    WriteOutTheta(z)

    if (nbrNodes > 1)
    {
        require(snow)
        require(rlecuyer)
        clusterString <- rep("localhost", nbrNodes)
        z$cl <- makeCluster(clusterString, type = "SOCK",
                            outfile = "cluster.out")
        clusterCall(z$cl, library, pkgname, character.only = TRUE)
        clusterCall(z$cl, storeinFRANstore,  FRANstore())
        clusterCall(z$cl, FRANstore)
        clusterCall(z$cl, initializeFRAN, z, model,
                    initC = TRUE, profileData=FALSE, returnDeps=FALSE)
        clusterSetupRNG(z$cl,
                        seed = as.integer(runif(6, max=.Machine$integer.max)))
    }

    z$scaleFactors <- rep(1, z$nGroup)
    ## z$returnDataFrame <- TRUE # chains come back as data frames not lists
    z$returnChains <- TRUE
    if (is.null(priorSigma))
    {
        z$priorSigma <- diag(z$pp) * 10000
    }
    else
    {
        z$priorSigma <- priorSigma
    }
    z$ratePositions <- lapply(z$rateParameterPosition, unlist)
    for (i in 1:z$nGroup)
    {
        use <- rep(FALSE, z$pp)
        use[z$ratePositions[[i]]] <- TRUE
        use[!z$basicRate] <- TRUE
        z$thetaMat[i, !use] <- NA
    }
    z$callGrid <- cbind(rep(1:z$nGroup, z$groupPeriods),
                        as.vector(unlist(sapply(z$groupPeriods,
                                                function(x) 1:x))))
    z
}
##@getDFRA algorithms do a few ML iterations and calculate a derivative matrix
getDFRA <- function(z, n)
{
    ## do n MLmodelsteps with the initial thetas and get
    ## derivs
    z$sdf <- array(0, dim=c(n, z$pp, z$pp))
    z$ssc <- matrix(0, nrow=n, ncol=z$pp)
    z$Deriv <- TRUE
    for (i in 1:n)
    {
        ans <- z$FRAN(z)
        z$sdf[i, , ] <- ans$dff
        z$ssc[i,  ] <- colSums(ans$fra)
    }
    dfra <- t(apply(z$sdf, c(2, 3), mean))
    ssc <- colMeans(z$ssc)
    z$dfra <- dfra
    lambda <- z$theta[z$basicRate]
    z$dfra[z$basicRate, ] <- z$dfra[z$basicRate,] * lambda
    z$dfra[, z$basicRate] <- z$dfra[, z$basicRate] * lambda
    diag(z$dfra)[z$basicRate] <- lambda *
        diag(dfra)[z$basicRate] + lambda * lambda * ssc[z$basicRate]
    ##print(z$dfra)
    z$dfra <- solve(z$dfra)
    ##print(z$dfra)
    ##$eS <- eigen(z$dfra, symmetric=TRUE, EISPACK=TRUE)
    ##z$ev <- z$eS$values
    ##z$covFactor <- z$eS$vectors %*% diag(sqrt(pmax(z$ev, 0)), z$pp)
    z
}
##@getLikelihood algorithms calculated likelihood for one chain
getLikelihood <- function(chain, nactors, lambda, simpleRates)
{
    loglik <- 0
    ncvals <- sapply(chain, function(x)x[[3]])
    nc <- nactors
    nc[] <- 0
    ncvals <- table(ncvals)
    nc[names(ncvals)] <- ncvals
    logChoiceProb <- sapply(chain, function(x)x[[9]])
    logOptionSetProb <- sapply(chain, function(x)x[[8]])
    loglik <- sum(logChoiceProb)  + sum(logOptionSetProb)
	##print(sum(logOptionSetProb))
    if (simpleRates)
    {
        loglik <- loglik - sum(nactors * lambda) + sum(nc * log(lambda))# -
		## sum(lfactorial(nc)) don't need factorial in bayes!
    }
    else
    {
        mu <- attr(chain, "mu")
        sigma <- sqrt(attr(chain, "sigma2"))
        finalReciprocalRate <- attr(chain, "finalReciprocalRate")
        loglik <- loglik + dnorm(1, mu, sigma, log=TRUE) +
            log(finalReciprocalRate)
    }
    loglik
}
##@getProbs algorithms calculates likelihood sum over nested list of chains
getProbs <- function(z, chain)
{
    sapply(1: length(chain), function(i)
       {
           groupChain <- chain[[i]]
           sum(sapply(1:length(groupChain), function(j)
                  {
                      periodChain <- groupChain[[j]]
                      theta1 <- z$thetaMat[i, ]
                      k <- z$rateParameterPosition[[i]][[j]]
                      getLikelihood(periodChain, z$nactors[[i]], theta1[k],
                                    z$simpleRates)
                  }
                      ))
       }
           )
}
##@flattenChains algorithms converts a nested list of chains to a single list
flattenChains <- function(zz)
{
        for (i in 1:length(zz)) ##group
        {
            for (j in 1:length(zz[[i]])) ## period
            {
                attr(zz[[i]][[j]], "group") <- i
                attr(zz[[i]][[j]], "period") <- j
            }
        }
    zz <- do.call(c, zz)
    zz
}
##@dmvnorm algorithms calculated multivariate normal density:
##inefficient: should not call mahalanobis and eigen with same sigma repeatedly
dmvnorm <- function(x, mean , sigma)
{
    if (is.vector(x))
    {
        x <- matrix(x, ncol=length(x))
    }
    distval <- mahalanobis(x, center=mean, cov=sigma)
    logdet <- sum(log(eigen(sigma, symmetric=TRUE, only.values=TRUE)$values))
    -(ncol(x) * log(2 * pi) + logdet + distval) / 2
}
