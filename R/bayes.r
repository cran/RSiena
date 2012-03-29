##/*****************************************************************************
## * SIENA: Simulation Investigation for Empirical Network Analysis
## *
## * Web: http://www.stats.ox.ac.uk/~snidjers/siena
## *
## * File: bayes.r
## *
## * Description: This file contains the code to run Bayesian simulation.
## * Many functions are defined within others to reduce copying of objects.
## *
## ****************************************************************************/
##@bayes Bayesian fit a Bayesian model
bayes <- function(data, effects, model, nwarm=100, nmain=100, nrunMHBatches=20,
                  plotit=FALSE, nbrNodes=1, dfra=NULL, n=10,
                  priorSigma=NULL, prevAns=NULL, clusterType=c("PSOCK", "FORK"),
				  getDocumentation=FALSE)
{
    ##@createStores internal bayes Bayesian set up stores
    createStores <- function()
    {
        npar <- length(z$theta)
        numberRows <- nmain * nrunMHBatches
        z$posteriorTot <<- matrix(0, nrow=z$nGroup, ncol=npar)
        z$posteriorMII <<- array(0, dim=c(z$nGroup, npar, npar))
        z$candidates <<- array(NA, dim=c(numberRows, z$nGroup, npar))
        z$acceptances <<- matrix(NA, ncol=z$nGroup, nrow=numberRows)
        z$MHacceptances <<- array(NA, dim=c(numberRows, z$nGroup,
									 z$nDependentVariables, 9))
        z$MHrejections <<- array(NA, dim=c(numberRows, z$nGroup,
									 z$nDependentVariables, 9))
        z$MHproportions <<- array(NA, dim=c(numberRows, z$nGroup,
									 z$nDependentVariables, 9))
    }
    ##@storeData internal bayes Bayesian put data in stores
    storeData <- function()
    {
        start <- z$sub + 1
        nrun <- nrow(z$parameters)
        end <- start + nrun - 1
        z$acceptances[start:end, ] <<- z$accepts
        z$candidates[start:end,, ] <<- z$parameters
        z$posteriorTot <<- z$posteriorTot + colSums(z$parameters)
        for (group in 1:z$nGroup)
        {
            for (i in dim(z$parameters)[1])
            {
                z$posteriorMII[group, , ] <<- z$posteriorMII[group, ,] +
                    outer(z$parameters[i, group, ], z$parameters[i, group, ])
            }
        }
        z$MHacceptances[start:end, , , ] <<- z$MHaccepts
        z$MHrejections[start:end, , , ] <<- z$MHrejects
        z$MHproportions[start:end, , , ] <<- z$MHaccepts /
            (z$MHaccepts + z$MHrejects)
        z$sub <<- z$sub + nrun

    }
	##@improveMH internal bayes Bayesian find scale factors
    improveMH <- function(tiny=1.0e-15, desired=40, maxiter=100,
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
            MCMCcycle(nrunMH=1, nrunMHBatches=100, change=FALSE)
            actual <- z$BayesAcceptances ## acceptances
            ans <- rescaleCGD(100)
            update <- number < 3
            z$scaleFactors[update] <<- z$scaleFactors[update] * ans[update]
            cat(actual, ans, z$scaleFactors, '\n')
            if (all(success) || iter == maxiter)
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
    }
	##@MCMCcycle internal Bayes do some loops of (MH steps and sample parameters)
	MCMCcycle <- function(nrunMH, nrunMHBatches, change=TRUE)
	{
		z$accepts <<- matrix(NA, nrow=z$nGroup, nrunMHBatches)
		z$parameters <<- array(NA, dim=c(nrunMHBatches, z$nGroup, z$pp))
		z$MHaccepts <<- array(NA, dim=c(nrunMHBatches, z$nGroup,
								  z$nDependentVariables, 9))
		z$MHrejects <<- array(NA, dim=c(nrunMHBatches, z$nGroup,
								  z$nDependentVariables, 9))
		z$MHaborts <<- array(NA, dim=c(nrunMHBatches, z$nGroup,
								 z$nDependentVariables, 9))
		storeNrunMH <- z$nrunMH
		z$nrunMH <<- nrunMH
		for (i in 1:nrunMHBatches)
		{
		#	cc <- proc.time()[1]

			ans <- z$FRAN(z, byGroup=TRUE, returnLoglik=TRUE, onlyLoglik=TRUE)
										#	c2 <- proc.time()[1]
										#	cat ('fran',c2-cc,'\n')
			z$loglik <<- ans$loglik
										#	cc <- proc.time()[1]
			sampleParameters(change)
										#	cc1 <- proc.time()[1]
										#	cat('samp',cc1-cc, '\n')
			z$accepts[, i] <<- z$accept
			z$parameters[i, , ] <<- z$thetaMat
			z$MHaccepts[i, , , ] <<-
				t(do.call(cbind,
						  tapply(ans$accepts, factor(z$callGrid[, 1]),
								 function(x)Reduce("+", x))))
			z$MHrejects[i, , , ] <<-
				t(do.call(cbind, tapply(ans$rejects, factor(z$callGrid[, 1]),
										function(x)Reduce("+", x))))
			z$MHaborts[i, , , ] <<- t(do.call(cbind,
											  tapply(ans$aborts,
													 factor(z$callGrid[, 1]),
													 function(x)Reduce("+", x))))
		}
		z$BayesAcceptances <<- rowSums(z$accepts)
		z$nrunMH <<- storeNrunMH
	}
	##@sampleParameters algorithms propose new parameters and accept them or not
	sampleParameters <- function(change=TRUE)
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
		logpOld <- z$loglik

		thetaNew[, z$basicRate] <- exp(thetaNew[, z$basicRate])
		z$thetaMat <<- thetaNew
		logpNew <- getProbabilitiesFromC(z)[[1]]
		proposalProbability <- priorNew - priorOld + logpNew - logpOld
		##cat(proposalProbability, priorNew, priorOld, logpNew, logpOld, '\n')
		z$accept <<- log(runif(length(proposalProbability))) <
			proposalProbability
		thetaOld[, z$basicRate] <- exp(thetaOld[, z$basicRate])
		if (!change)
		{
			z$thetaMat <<- thetaOld
		}
		else
		{
			##print(z$thetaMat)
			z$thetaMat[!z$accept, ] <<- thetaOld[!z$accept, ]
		}
##		print(thetaNew)
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
	ctime <- proc.time()[1]

	z <- initializeBayes(data, effects, model, nbrNodes, priorSigma,
                         prevAns=prevAns, clusterType=clusterType)
    createStores()

    z$sub <- 0

    if (is.null(z$dfra) && is.null(dfra))
    {
        z <- getDFRA(z, n)
    }
    else
    {
        if (!is.null(dfra))
        {
            z$dfra <- dfra
        }
		else
		{
			if (is.null(z$sf))
			{
				stop("need some scores to scale dfra")
			}
			z$dfra <- scaleDfra(z)
		}

    }
	ctime1 <- proc.time()[1]
	cat(ctime1-ctime,'\n')
	improveMH()
	ctime2<- proc.time()[1]

	cat('improvMh', ctime2-ctime1,'\n')

    if (plotit)
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
        MCMCcycle(nrunMH=4, nrunMHBatches=20)
    }
	print('endof warm')
	ctime3<- proc.time()[1]

 	cat('warm', ctime3-ctime2,'\n')
    for (ii in 1:nmain)
    {
		MCMCcycle(nrunMH=z$nrunMH, nrunMHBatches=nrunMHBatches)
		storeData()
		ctime4<- proc.time()[1]
		cat('main', ii, ctime4-ctime3,'\n')
		ctime3 <- ctime4

        if (ii %% 10 == 0 && plotit) ## do some plots
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
			if (z$nGroup > 1)
			{
				varcall <- paste("~ ", varnames,  " | Group", sep="",
								 collapse="")
			}
			else
			{
				varcall <- paste("~ ", varnames,  sep="", collapse="")

			}
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
			if (z$nGroup > 1)
			{
				varcall <- paste(varnames,  "~ 1:", ii *
								 nrunMHBatches * z$nGroup,
                             " | Group", sep="", collapse="")
			}
			else
			{
				varcall <- paste(varnames,  "~ 1:", ii *
								 nrunMHBatches * z$nGroup,
								 sep="", collapse="")
			}
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

##@initializeBayes algorithms do set up for Bayesian model
initializeBayes <- function(data, effects, model, nbrNodes, priorSigma,
                            prevAns, clusterType=c("PSOCK", "FORK"))
{
    ## initialise
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
		##seed <- model$randomSeed
    }
	else
	{
		if (exists(".Random.seed"))
		{
			rm(.Random.seed, pos=1)
		}
		newseed <- trunc(runif(1) * 1000000)
		set.seed(newseed)  ## get R to create a random number seed for me.
		##seed <- NULL
	}
   	z$FRAN <- getFromNamespace(model$FRANname, pkgname)
    z <- initializeFRAN(z, model, data=data, effects=effects,
                prevAns=prevAns, initC=FALSE, onlyLoglik=TRUE)
	z$basicRate <- z$effects$basicRate
    z$nGroup <- z$f$nGroup
	is.batch(TRUE)

    WriteOutTheta(z)

    if (nbrNodes > 1 && z$observations > 1)
    {
        require(parallel)
		clusterType <- match.arg(clusterType)
		if (clusterType == "PSOCK")
		{
        clusterString <- rep("localhost", nbrNodes)
        z$cl <- makeCluster(clusterString, type = "PSOCK",
                            outfile = "cluster.out")
		}
		else
		{
			z$cl <- makeCluster(nbrNodes, type = "FORK",
								outfile = "cluster.out")
		}
        clusterCall(z$cl, library, pkgname, character.only = TRUE)
        clusterCall(z$cl, storeinFRANstore,  FRANstore())
        clusterCall(z$cl, FRANstore)
        clusterCall(z$cl, initializeFRAN, z, model,
                    initC = TRUE, profileData=FALSE, returnDeps=FALSE)
		clusterSetRNGStream(z$cl, iseed = as.integer(runif(1,
								max=.Machine$integer.max)))
    }

    z$scaleFactors <- rep(1, z$nGroup)
    ## z$returnDataFrame <- TRUE # chains come back as data frames not lists
    z$returnChains <- FALSE
    if (is.null(priorSigma))
    {
        z$priorSigma <- diag(z$pp) * 10000
    }
    else
    {
        z$priorSigma <- priorSigma
    }
  	groupPeriods <- attr(z$f, "groupPeriods")
    netnames <- z$f$depNames
	z$rateParameterPosition <-
        lapply(1:z$nGroup, function(i, periods, data)
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
           }, periods=groupPeriods - 1, data=z$f[1:z$nGroup]
               )
	z$ratePositions <- lapply(z$rateParameterPosition, unlist)
    for (i in 1:z$nGroup)
    {
        use <- rep(FALSE, z$pp)
        use[z$ratePositions[[i]]] <- TRUE
        use[!z$basicRate] <- TRUE
        z$thetaMat[i, !use] <- NA
    }
    z
}
##@getDFRA algorithms do a few ML iterations and calculate a derivative matrix
getDFRA <- function(z, n)
{
    ## do n MLmodelsteps with the initial thetas and get
    ## derivs
    z$sdf <- vector("list", n)
    z$sf <- matrix(0, nrow=n, ncol=z$pp)
    z$Deriv <- TRUE
    for (i in 1:n)
    {
        ans <- z$FRAN(z)
        z$sdf[[i]] <- ans$dff
        z$sf[i,  ] <- colSums(ans$fra)
    }
	dfra <- t(as.matrix(Reduce("+", z$sdf) / length(z$sdf)))
    z$dfra <- dfra
	z$dfra <- scaleDfra(z)
	z$Deriv <- FALSE
    z
}
scaleDfra <- function(z)
{
    lambda <- z$theta[z$basicRate]
	dfra <- z$dfra
    z$dfra[z$basicRate, ] <- z$dfra[z$basicRate,] * lambda
    z$dfra[, z$basicRate] <- z$dfra[, z$basicRate] * lambda
    diag(z$dfra)[z$basicRate] <- lambda *
		diag(dfra)[z$basicRate] + lambda * lambda *
			colMeans(z$sf)[z$basicRate]
	chol2inv(chol(z$dfra))
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

##@getProbabilitiesFromC bayes gets loglik from chains in C
getProbabilitiesFromC <- function(z, index=1, getScores=FALSE)
{
	## expects maximum likelihood parallelisations
    f <- FRANstore()

	callGrid <- z$callGrid
    ## z$int2 is the number of processors if iterating by period, so 1 means
    ## we are not. Can only parallelize by period1
    if (nrow(callGrid) == 1)
    {
		theta <- z$thetaMat[1,]
        ans <- .Call("getChainProbabilities", PACKAGE = pkgname, f$pData,
					 f$pModel, as.integer(1), as.integer(1),
					 as.integer(index), f$myeffects, theta, getScores)
        anss <- list(ans)
	}
    else
    {
        if (z$int2 == 1 )
        {
            anss <- apply(callGrid, 1,
						  doGetProbabilitiesFromC, z$thetaMat, index, getScores)
        }
        else
        {
            use <- 1:(min(nrow(callGrid), z$int2))
            anss <- parRapply(z$cl[use], callGrid,
							  doGetProbabilitiesFromC, z$thetaMat, index,
							  getScores)
        }
    }
	ans <- list()
	ans[[1]] <- sum(sapply(anss, "[[", 1))
	if (getScores)
	{
		ans[[2]] <- rowSums(sapply(anss, "[[", 2))
	}
	ans[[3]] <- sapply(anss, "[[", 3)
	ans
}

##@doGetProbabilitiesFromC Maximum likelihood
doGetProbabilitiesFromC <- function(x, thetaMat, index, getScores)
{
    f <- FRANstore()
	theta <- thetaMat[x[1], ]
    .Call("getChainProbabilities", PACKAGE = pkgname, f$pData,
		  f$pModel, as.integer(x[1]), as.integer(x[2]),
		  as.integer(index), f$myeffects, theta, getScores)

}
