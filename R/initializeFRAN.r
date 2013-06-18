#/******************************************************************************
# * SIENA: Simulation Investigation for Empirical Network Analysis
# *
# * Web: http://www.stats.ox.ac.uk/~snijders/siena
# *
# * File: initializeFRAN.r
# *
# * Description: This module contains the code for setting up the data in C++
# * and for ML making the initial chain.
# *****************************************************************************/
##@initializeFRAN siena07 reformat data and send to C. get targets.
initializeFRAN <- function(z, x, data, effects, prevAns=NULL, initC,
						   profileData=FALSE, returnDeps=FALSE,
						   returnChains=FALSE, byGroup=FALSE,
						   returnDataFrame=FALSE, byWave=FALSE,
						   returnLoglik=FALSE, onlyLoglik=FALSE)
{
	z$effectsName <- deparse(substitute(effects))
	## fix up the interface so can call from outside robmon framework
    if (is.null(z$FinDiff.method))
    {
        z$FinDiff.method <- FALSE
    }
    if (is.null(z$int))
    {
        z$int <- 1
    }
    if (is.null(z$int2))
    {
        z$int2 <- 1
    }
    if (!initC) ## ie first time round
    {
        if (!inherits(data,"siena"))
        {
            stop("not valid siena data object")
        }
        ## check the effects object
        defaultEffects <- getEffects(data)
        if (is.null(effects))
        {
            effects <- defaultEffects
        }
        else
        {
            ## todo check that the effects match the data dependent variables
            userlist <- apply(effects[effects$include,], 1, function(x)
                              paste(x[c("name", "effectName",
                                        "type", "groupName")],
                                    collapse="|"))
            deflist <- apply(defaultEffects, 1, function(x)
                             paste(x[c("name", "effectName",
                                       "type", "groupName")],
                                   collapse="|"))
            if (!all(userlist %in% deflist))
            {
                bad <- which(!(userlist %in% deflist))
                print(userlist[bad])
                stop("invalid effect requested: see above ")
            }
        }
        if (!inherits(effects, "data.frame"))
        {
            stop("effects is not a data.frame")
        }
        if (x$useStdInits)
        {
            if (any(effects$effectName != defaultEffects$effectName))
            {
                stop("Cannot use standard initialisation with a ",
                     "different effect list")
            }
            effects$initialValue <- defaultEffects$initialValue
        }
        ## get data object into group format to save coping with two
        ## different formats
        if (inherits(data, "sienaGroup"))
        {
            nGroup <- length(data)
        }
        else
        {
            nGroup <- 1
            data <- sienaGroupCreate(list(data), singleOK=TRUE)
        }
		## add any effects needed for time dummies
        tmp <- sienaTimeFix(effects, data)
        data <- tmp$data
        effects <- tmp$effects
		if (!x$useStdInits)
        {
			if (!is.null(prevAns) && inherits(prevAns, "sienaFit"))
			{
				effects <- updateTheta(effects, prevAns)
			}
        }
		## add any effects needed for settings model
# this now is replaced by adding the settings in getEffects,
# which is the more logical place.
# If all works, this can be deleted,
# and also the function addSettingsEffects can be deleted.
# I used this function as a template for the change to getEffects.
# I wonder why the next 8 lines cannot be dropped;
# gives error message "cannot find setting col".
		if (!is.null(x$settings))
		{
			effects <- addSettingsEffects(effects, x)
		}
		else
		{
			effects$setting <- rep("", nrow(effects))
		}
        ## find any effects not included which are needed for interactions
        tmpEffects <- effects[effects$include, ]
        interactionNos <- unique(c(tmpEffects$effect1, tmpEffects$effect2,
                                   tmpEffects$effect3))
        interactionNos <- interactionNos[interactionNos > 0]
        interactions <- effects$effectNumber %in%
                                          interactionNos
        effects$requested <- effects$include
        requestedEffects <- effects[effects$include, ]

        effects$include[interactions] <- TRUE
        effects <- effects[effects$include, ]

        ## split and rejoin both versions before continuing
        depvarnames <- names(data[[1]]$depvars)

        effects1 <- split(requestedEffects, requestedEffects$name)
        effects1order <- match(depvarnames, names(effects1))
        requestedEffects <- do.call(rbind, effects1[effects1order])
        row.names(requestedEffects) <- 1:nrow(requestedEffects)

        effects1 <- split(effects, effects$name)
        effects1order <- match(depvarnames, names(effects1))
        effects <- do.call(rbind, effects1[effects1order])
        row.names(effects) <- 1:nrow(effects)

        ## now set up z provisionally
        z$theta <- requestedEffects$initialValue
        z$fixed <- requestedEffects$fix
        z$test <- requestedEffects$test
        z$pp <- length(z$test)
        z$posj <- rep(FALSE, z$pp)
        z$posj[requestedEffects$basicRate] <- TRUE
        z$BasicRateFunction <- z$posj

        ## sort out names of user specified interaction effects
        effects <- fixUpEffectNames(effects)
        ## copy interaction names to the requested effects
        requestedEffects$effectName <- effects[effects$requested,
                                               "effectName"]
        requestedEffects$functionName <- effects[effects$requested,
                                                 "functionName"]

		if (attr(data, "compositionChange"))
		{
			if (x$maxlike)
			{
				stop("Not possible to use maximum likelihood estimation ",
					 "with composition change")
			}
		}

        ## if not specified whether conditional or nor, set to conditional
        ## iff there is only one dependent variable (therefore number 1)
        ## and not maxlike
        if (is.na(x$cconditional))
        {
            x$cconditional <- !x$maxlike && (length(depvarnames) == 1)
            if (x$cconditional)
            {
                x$condvarno <- 1
            }
        }
        types <- sapply(data[[1]]$depvars, function(x) attr(x, "type"))
        ## now check if conditional estimation is OK and copy to z if so
        z$cconditional <- FALSE
        if ((x$cconditional) & (!(attr(data, "compositionChange"))))
        {
            if (x$maxlike)
            {
                stop("Conditional estimation is not possible with",
                     "maximum likelihood method")
            }
            ##  if (nets == 1) not sure if this is necessary
            ##  {
            z$cconditional <- TRUE
            ## find the conditioning variable
            observations <- attr(data, "observations")
            ## this is actual number of waves to process
            if (x$condname != "")
            {
                z$condvarno <- match(x$condname, attr(data, "netnames"))
                z$condname <- x$condname
            }
            else
            {
                z$condvarno <- x$condvarno
                z$condname <- attr(data, "netnames")[x$condvarno]
            }
            z$condtype <- attr(data, "types")[z$condvarno]
            if (z$condtype == "oneMode")
                z$symmetric  <-  attr(data, "symmetric")[[z$condvarno]]
            else
                z$symmetric <- FALSE
            ## find the positions of basic rate effects for this network
            z$condvar <-
                (1:nrow(requestedEffects))[requestedEffects$name==
                                           z$condname][1:observations]
            z$theta<- z$theta[-z$condvar]
            z$fixed<- z$fixed[-z$condvar]
            z$test<- z$test[-z$condvar]
            z$pp<- z$pp - length(z$condvar)
            z$scale<- z$scale[-z$condvar]
            z$BasicRateFunction <- z$posj[-z$condvar]
            z$posj <- z$posj[-z$condvar]
            z$theta[z$posj] <-
                z$theta[z$posj] / requestedEffects$initialValue[z$condvar]
        }

        ## unpack data and put onto f anything we may need next time round.
        f <- lapply(data, function(xx, x) unpackData(xx, x), x=x)

        attr(f, "netnames") <- attr(data, "netnames")
        attr(f, "symmetric") <- attr(data, "symmetric")
        attr(f, "allUpOnly") <- attr(data, "allUpOnly")
        attr(f, "allDownOnly") <- attr(data, "allDownOnly")
        attr(f, "allHigher") <- attr(data, "allHigher")
        attr(f, "allDisjoint") <- attr(data, "allDisjoint")
        attr(f, "allAtLeastOne") <- attr(data, "allAtLeastOne")
        attr(f, "anyUpOnly") <- attr(data, "anyUpOnly")
        attr(f, "anyDownOnly") <- attr(data, "anyDownOnly")
        attr(f, "anyHigher") <- attr(data, "anyHigher")
        attr(f, "anyDisjoint") <- attr(data, "anyDisjoint")
        attr(f, "anyAtLeastOne") <- attr(data, "anyAtLeastOne")
        attr(f, "types") <- attr(data, "types")
        attr(f, "observations") <- attr(data, "observations")
        attr(f, "compositionChange") <- attr(data, "compositionChange")
        attr(f, "exooptions") <- attr(data, "exooptions")
        attr(f, "groupPeriods") <- attr(data, "groupPeriods")
        attr(f, "periodNos") <- attr(data, "periodNos")
        attr(f, "numberNonMissingNetwork") <-
            attr(data, "numberNonMissingNetwork")
        attr(f, "numberMissingNetwork") <- attr(data, "numberMissingNetwork")
        attr(f, "numberNonMissingBehavior") <-
            attr(data, "numberNonMissingBehavior")
        attr(f, "numberMissingBehavior") <- attr(data, "numberMissingBehavior")

		##  attr(f, "totalMissings") <- attr(data, "totalMissings")

        if (x$maxlike && x$FinDiff.method)
        {
            stop("Finite difference method for derivatives not available",
                 "with Maximum likelihood method")
        }
        ## if any networks symmetric must use finite differences and not maxlike
        ## scores are now available (30.10.10) but still not maxlike.
        ## check model type: default for symmetric is type 2 (forcing model).
        syms <- attr(data,"symmetric")
        z$FinDiffBecauseSymmetric <- FALSE
        z$modelType <- x$modelType
        if (any(!is.na(syms) & syms))
        {
            ##     z$FinDiff.method <- TRUE
            ##     z$FinDiffBecauseSymmetric <- TRUE
            if (x$maxlike)
            {
                stop("Maximum likelihood method not implemented",
                     "for symmetric networks")
            }
            if (x$modelType == 1)
            {
                z$modelType <- 2
            }
        }
        if (z$cconditional)
        {
            attr(f, "change") <-
                sapply(f, function(xx)attr(xx$depvars[[z$condname]],
                                           "distance"))
            attr(f,"condEffects") <- requestedEffects[z$condvar,]
            effcondvar <-
                (1:nrow(effects))[effects$name==
                                  z$condname][1:observations]
            effects <- effects[-effcondvar, ]
            requestedEffects <- requestedEffects[-z$condvar,]
        }
        ## use previous dfra only if everything matches including method
        if (!is.null(prevAns) && inherits(prevAns, "sienaFit") &&
			prevAns$maxlike == z$maxlike)
        {
            if (!is.null(prevAns$dfra) &&
				 nrow(prevAns$dfra) == nrow(requestedEffects) &&
                 all(rownames(prevAns$dfra) == requestedEffects$shortName) &&
				 all(prevAns$fix == requestedEffects$fix) &&
				 !is.null(prevAns$sf))
            {
                z$haveDfra <- TRUE
                z$dfra <- prevAns$dfra
                z$dinv <- prevAns$dinv
				# z$dinvv must not be taken from prevAns,
				# because the value of diagonalize
				# is defined in x and may have changed.
				# Therefore here we copy the corresponding lines
				# from phase1.r.
				if (!x$diagg)
				{
					# Partial diagonalization of derivative matrix
					# for use if 0 < x$diagonalize < 1.
					temp <- (1-x$diagonalize)*z$dfra +
							x$diagonalize*diag(diag(z$dfra))
					temp[z$fixed, ] <- 0.0
					temp[, z$fixed] <- 0.0
					diag(temp)[z$fixed] <- 1.0
					# Invert this matrix
					z$dinvv <- solve(temp)
				}
				z$sf <- prevAns$sf
				# check for backward compatibility with pre-1.1-220 versions:
				if (is.null(prevAns$regrCoef))
				{
					z$regrCoef <- rep(0, z$pp)
					z$regrCor <- rep(0, z$pp)
				}
				else
				{
					z$regrCoef <- prevAns$regrCoef
					z$regrCor <- prevAns$regrCor
				}
            }
        }
        z$effects <- effects
        z$requestedEffects <- requestedEffects
    }
    else ## initC, i.e just send already set up data into new processes
    {
        f <- FRANstore()
        ## Would like f to be just the data objects plus the attributes
        ## but need the effects later. Also a few other things,
        ## which probably could be attributes but are not!
        ## They will be automatically removed: if needed they must be readded
        ff <- f
        nGroup <- f$nGroup
        f[(nGroup + 1): length(f)] <- NULL
    }
    pData <- .Call("setupData", PACKAGE=pkgname,
                   lapply(f, function(x)(as.integer(x$observations))),
                   lapply(f, function(x)(x$nodeSets)))
	ans <- .Call("OneMode", PACKAGE=pkgname,
                 pData, lapply(f, function(x)x$nets))
    ans <- .Call("Bipartite", PACKAGE=pkgname,
                 pData, lapply(f, function(x)x$bipartites))
    ans <- .Call("Behavior", PACKAGE=pkgname,
                 pData, lapply(f, function(x)x$behavs))
    ans <-.Call("ConstantCovariates", PACKAGE=pkgname,
                pData, lapply(f, function(x)x$cCovars))
    ans <-.Call("ChangingCovariates", PACKAGE=pkgname,
                pData, lapply(f, function(x)x$vCovars))
	ans <-.Call("DyadicCovariates", PACKAGE=pkgname,
                pData, lapply(f, function(x)x$dycCovars))
	ans <-.Call("ChangingDyadicCovariates", PACKAGE=pkgname,
                pData, lapply(f, function(x)x$dyvCovars))
	ans <-.Call("ExogEvent", PACKAGE=pkgname,
                pData, lapply(f, function(x)x$exog))
   ## split the names of the constraints
    higher <- attr(f, "allHigher")
    disjoint <- attr(f, "allDisjoint")
    atLeastOne <- attr(f, "allAtLeastOne")
    froms <- sapply(strsplit(names(higher), ","), function(x)x[1])
    tos <- sapply(strsplit(names(higher), ","), function(x)x[2])
    ans <- .Call("Constraints", PACKAGE=pkgname,
                 pData, froms[higher], tos[higher],
                 froms[disjoint], tos[disjoint],
                 froms[atLeastOne], tos[atLeastOne])

    ##store the address
    f$pData <- pData
    ## register a finalizer
    ans <- reg.finalizer(f$pData, clearData, onexit = FALSE)

    if (!initC)
    {
        storage.mode(effects$parm) <- "integer"
        storage.mode(effects$group) <- "integer"
        storage.mode(effects$period) <- "integer"
        effects$effectPtr <- rep(NA, nrow(effects))
        splitFactor <- factor(effects$name, levels=attr(f, "netnames"))
        if (!all(attr(f,"netnames") %in% effects$name))
        {
            stop("Must have at least one effect for each dependent variable")
        }
        myeffects <- split(effects, splitFactor)
        myCompleteEffects <- myeffects
        ## remove interaction effects and save till later
        basicEffects <-
            lapply(myeffects, function(x)
               {
                   x[!x$shortName %in% c("unspInt", "behUnspInt"), ]
               }
                   )
         basicEffectsl <-
            lapply(myeffects, function(x)
               {
                   !x$shortName %in% c("unspInt", "behUnspInt")
               }
                   )

        interactionEffects <-
            lapply(myeffects, function(x)
               {
                   x[x$shortName %in% c("unspInt", "behUnspInt"), ]
               }
                   )
         interactionEffectsl <-
            lapply(myeffects, function(x)
               {
                   x$shortName %in% c("unspInt", "behUnspInt")
               }
                   )
       ## store effects objects as we may need to recreate them
        f$interactionEffects <- interactionEffects
        f$basicEffects <- basicEffects
        f$interactionEffectsl <- interactionEffectsl
        f$basicEffectsl <- basicEffectsl
    }
    else
    {
        myCompleteEffects <- ff$myCompleteEffects
        basicEffects <- ff$basicEffects
        interactionEffects <- ff$interactionEffects
        basicEffectsl <- ff$basicEffectsl
        interactionEffectsl <- ff$interactionEffectsl
        types <- ff$types
    }
	ans <- .Call("effects", PACKAGE=pkgname, pData, basicEffects)
    pModel <- ans[[1]][[1]]
    for (i in seq(along=(ans[[2]]))) ## ans[[2]] is a list of lists of
        ## pointers to effects. Each list corresponds to one
        ## dependent variable
    {
        effectPtr <- ans[[2]][[i]]
        basicEffects[[i]]$effectPtr <- effectPtr
        interactionEffects[[i]]$effect1 <-
            basicEffects[[i]]$effectPtr[match(interactionEffects[[i]]$effect1,
                                              basicEffects[[i]]$effectNumber)]
        interactionEffects[[i]]$effect2 <-
            basicEffects[[i]]$effectPtr[match(interactionEffects[[i]]$effect2,
                                              basicEffects[[i]]$effectNumber)]
        interactionEffects[[i]]$effect3 <-
            basicEffects[[i]]$effectPtr[match(interactionEffects[[i]]$effect3,
                                              basicEffects[[i]]$effectNumber)]
    }
    ans <- .Call("interactionEffects", PACKAGE=pkgname,
                 pModel, interactionEffects)
    ## copy these pointers to the interaction effects and then insert in
    ## effects object in the same rows for later use
    for (i in 1:length(ans[[1]])) ## ans is a list of lists of
        ## pointers to effects. Each list corresponds to one
        ## dependent variable
    {
        if (nrow(interactionEffects[[i]]) > 0)
        {
            effectPtr <- ans[[1]][[i]]
            interactionEffects[[i]]$effectPtr <- effectPtr
        }
        myCompleteEffects[[i]][basicEffectsl[[i]], ] <- basicEffects[[i]]
        myCompleteEffects[[i]][interactionEffectsl[[i]],] <-
            interactionEffects[[i]]
        ##myeffects[[i]] <- myeffects[[i]][order(myeffects[[i]]$effectNumber),]
    }
    ## remove the effects only created as underlying effects
    ## for interaction effects. first store the original for use next time
    myeffects <- lapply(myCompleteEffects, function(x)
                    {
                        x[x$requested, ]
                    }
                        )
    if (!initC)
    {
        ans <- .Call("getTargets", PACKAGE=pkgname, pData, pModel, myeffects,
					 z$parallelTesting)
		##stop("done")
		## create a grid of periods with group names in case want to
		## parallelize using this or to access chains easily
		groupPeriods <- attr(f, "groupPeriods")
		z$callGrid <- cbind(rep(1:nGroup, groupPeriods - 1),
							as.vector(unlist(sapply(groupPeriods - 1,
													function(x) 1:x))))
        if (!x$maxlike)
        {
            z$targets <- rowSums(ans)
            z$targets2 <- ans
        }
        else
        {
            z$targets <- rep(0, z$pp)
            z$targets2 <- ans
            z$targets2[] <- 0
            z$maxlikeTargets <- rowSums(ans)
            z$maxlikeTargets2 <- ans
            z$mult <- x$mult
            z$nrunMH <- floor(z$mult *
					colSums(z$maxlikeTargets2[z$requestedEffects$basicRate,
												   , drop=FALSE ]))
			## make the number pretty
			z$nrunMH <- ifelse (z$nrunMH > 100,
								 round(z$nrunMH / 100) * 100, z$nrunMH)
            ##thetaMat is to allow different thetas for each group in Bayes
            z$thetaMat <- matrix(z$theta, nrow=nGroup, ncol=z$pp, byrow=TRUE)
        }
    }

    ##store address of model
    f$pModel <- pModel
    ans <- reg.finalizer(f$pModel, clearModel, onexit = FALSE)
    if (x$MaxDegree == 0 || is.null(x$MaxDegree))
    {
        MAXDEGREE <-  NULL
    }
    else
    {
        MAXDEGREE <- x$MaxDegree
        storage.mode(MAXDEGREE) <- "integer"
    }
    if (z$cconditional)
    {
        CONDVAR <- z$condname
        CONDTARGET <- attr(f, "change")
        ##   cat(CONDTARGET, "\n")
    }
    else
    {
        CONDVAR <- NULL
        CONDTARGET <- NULL
    }

	simpleRates <- TRUE
	if (any(!z$effects$basicRate & z$effects$type =="rate"))
	{
		## browser()
		simpleRates <- FALSE
	}
	z$simpleRates <- simpleRates
	ans <- .Call("setupModelOptions", PACKAGE=pkgname,
                 pData, pModel, MAXDEGREE, CONDVAR, CONDTARGET,
                 profileData, z$parallelTesting, x$modelType, z$simpleRates)
    if (x$maxlike)
    {
        if (!initC)
        {
            ## set up chains and do initial steps

            types <- attr(f, "types")
            nbrNonMissNet <- attr(f, "numberNonMissingNetwork")
            nbrMissNet <- attr(f, "numberMissingNetwork")
            nbrNonMissBeh <- attr(f, "numberNonMissingBehavior")
            nbrMissBeh <- attr(f, "numberMissingBehavior")

            if (sum(nbrMissNet + nbrNonMissNet) > 0)
            {
                z$prmin <- nbrMissNet/ (nbrMissNet + nbrNonMissNet)
            }
            else
            {
                z$prmin <- rep(0, length(nbrMissNet))
            }
            if (sum(nbrMissBeh + nbrNonMissBeh) > 0)
            {
                z$prmib <-   nbrMissBeh/ (nbrMissBeh + nbrNonMissBeh)
            }
            else
            {
                z$prmib <- rep(0, length(nbrMissBeh))
            }
            z$probs <- c(x$pridg, x$prcdg, x$prper, x$pripr, x$prdpr, x$prirms,
                         x$prdrms)
            ans <- .Call("mlMakeChains", PACKAGE=pkgname, pData, pModel,
                         z$probs, z$prmin, z$prmib,
                         x$minimumPermutationLength,
                         x$maximumPermutationLength,
                         x$initialPermutationLength)
            f$minimalChain <- ans[[1]]
            f$chain <- ans[[2]]
        }
        else ## set up the initial chains in the sub processes
        {
			ans <- .Call("mlInitializeSubProcesses",
                         PACKAGE=pkgname, pData, pModel,
                         z$probs, z$prmin, z$prmib,
                         x$minimumPermutationLength,
                         x$maximumPermutationLength,
                         x$initialPermutationLength, ff$chain)
            f$chain <- ff$chain
       }
    }
	f$simpleRates <- simpleRates
    f$myeffects <- myeffects
    f$myCompleteEffects <- myCompleteEffects
    if (!initC)
    {
        if (is.null(z$print) || z$print)
        {
            DataReport(z, x, f)
        }
        f$randomseed2 <- z$randomseed2
    }
    else
    {
        f$randomseed2 <- ff$randomseed2
    }
    f$observations <- attr(f, "observations") + 1
    f$depNames <- names(f[[1]]$depvars)
    f$groupNames <- names(f)[1:nGroup]
    f$nGroup <- nGroup
    f$types <- types
    if (!initC)
    {
        z$f <- f
        ##z <- initForAlgorithms(z)
        z$periodNos <- attr(data, "periodNos")
        z$f$myeffects <- NULL
        z$f$myCompleteEffects <- NULL
		if (!returnDeps)
		{
			z$f[1:nGroup] <- NULL
		}
    }
    if (initC || (z$int == 1 && z$int2 == 1 &&
                  (is.null(z$nbrNodes) || z$nbrNodes == 1)))
    {
        f[1:nGroup] <- NULL
    }
    FRANstore(f) ## store f in FRANstore
    z$returnDeps <- returnDeps
    z$returnDepsStored <- returnDeps
    z$observations <- f$observations
	z$returnChains <- returnChains
	z$byGroup <- byGroup
	z$byWave <- byWave
	z$returnDataFrame <- returnDataFrame
    z$nDependentVariables <- length(z$f$depNames)
	if (initC)
	{
		NULL
	}
	else
	{
		z
	}
}
##@createEdgeLists siena07 Reformat data for C++
createEdgeLists<- function(mat, matorig, bipartite)
{
    ## mat1 is basic values, with missings and structurals replaced
    tmp <- lapply(1 : nrow(mat), function(x, y)
              {
                  mymat <- matrix(0, nrow = sum(y[x, ] > 0), ncol = 3)
                  mymat[, 1] <- x
                  mymat[, 2] <- which(y[x, ] != 0)
                  mymat[, 3] <- y[x, mymat[, 2]]
                  mymat
              }, y = mat)
    mat1 <- do.call(rbind, tmp)
    ## mat2 reverts to matorig to get the missing values
    tmp <- lapply(1 : nrow(matorig), function(x, y)
              {
                  mymat <- matrix(0, nrow = sum(is.na(y[x, ])), ncol = 3)
                  mymat[, 1] <- x
                  mymat[, 2] <- which(is.na(y[x, ]))
                  mymat[, 3] <- 1
                  mymat
              }, y = matorig)
    mat2 <- do.call(rbind, tmp)
    ## remove the diagonal if not bipartite
    if (!bipartite)
    {
        mat2 <- mat2[mat2[, 1] != mat2[, 2], , drop=FALSE]
    }
    ## mat3 structurals
    struct <- mat1[,3] %in% c(10, 11)
    mat1[struct, 3] <- mat1[struct,3] - 10
    mat3 <- mat1[struct, , drop=FALSE]
    mat3[, 3] <- 1
    mat1 <- mat1[!mat1[,3] == 0, , drop=FALSE] ##remove any zeros just created
    ##fix up storage mode to be integer
    storage.mode(mat1) <- "integer"
    storage.mode(mat2) <- "integer"
    storage.mode(mat3) <- "integer"
    ## add attribute of size
    if (bipartite)
    {
        attr(mat1, "nActors") <- c(nrow(mat), ncol(mat))
        attr(mat2, "nActors") <- c(nrow(mat), ncol(mat))
        attr(mat3, "nActors") <- c(nrow(mat), ncol(mat))
    }
    else
    {
        attr(mat1, "nActors") <- nrow(mat)
        attr(mat2, "nActors") <- nrow(mat)
        attr(mat3, "nActors") <- nrow(mat)
    }

    list(mat1 = t(mat1), mat2 = t(mat2), mat3 = t(mat3))
}
##@createCovarEdgeLists siena07 Reformat data for C++
createCovarEdgeList<- function(mat, matorig)
{
    tmp <- lapply(1 : nrow(mat), function(x, y)
              {
                  mymat <- matrix(0, nrow = sum(y[x, ] != 0), ncol = 3)
                  mymat[, 1] <- x
                  mymat[, 2] <- which(y[x, ] != 0)
                  mymat[, 3] <- y[x, mymat[, 2]]
                  mymat
              }, y = mat)
    mat1 <- do.call(rbind, tmp)
    ##mat2 reverts to matorig to get the missing values
    tmp <- lapply(1 : nrow(matorig), function(x, y)
              {
                  mymat <- matrix(0, nrow = sum(is.na(y[x, ])), ncol = 3)
                  mymat[, 1] <- x
                  mymat[, 2] <- which(is.na(y[x, ]))
                  mymat[, 3] <- 1
                  mymat
              }, y = matorig)
    mat2 <- do.call(rbind, tmp)
    ## add attribute of size
    attr(mat1, "nActors1") <- nrow(mat)
    attr(mat1, "nActors2") <- ncol(mat)
    list(mat1=t(mat1), mat2=t(mat2))
}
##@unpackOneMode siena07 Reformat data for C++
unpackOneMode <- function(depvar, observations, compositionChange)
{
    edgeLists <- vector("list", observations)
    networks <- vector("list", observations)
    actorSet <- attr(depvar, "nodeSet")
    compActorSets <- sapply(compositionChange, function(x)attr(x, "nodeSet"))
    thisComp <- match(actorSet, compActorSets)
    compChange <- !is.na(thisComp)
    if (compChange)
    {
        action <- attr(compositionChange[[thisComp]], "action")
        ccOption <- attr(compositionChange[[thisComp]], "ccOption")
    }
    else
    {
        ccOption <- 0
        action <- matrix(0, nrow=attr(depvar, "netdims")[1], ncol=observations)
    }
    ## sort out composition change
    ##      convertToStructuralZeros()?
    sparse <- attr(depvar, "sparse")
	allowOnly <- attr(depvar, "allowOnly")
    if (sparse)
    {
        ## require(Matrix)
        ## have a list of sparse matrices in triplet format
        ## with missings and structurals embedded and 0 based indices!
        netmiss <- vector("list", observations)
        for (i in 1:observations)
        {
            ## extract this matrix
            networks[[i]] <- depvar[[i]]
            nActors <- nrow(depvar[[i]])
            ## stop if any duplicates
            netmat <- cbind(networks[[i]]@i+1, networks[[i]]@j+1,
                            networks[[i]]@x)
            if (any(duplicated(netmat[, 1:2])))
            {
                stop("duplicate entries in sparse matrix")
            }
            ## extract missing entries
            netmiss[[i]] <- netmat[is.na(netmat[,3]), , drop = FALSE]
            netmiss[[i]] <-
                netmiss[[i]][netmiss[[i]][, 1] != netmiss[[i]][, 2], ,
                             drop=FALSE]
            ## carry forward missing values if any
            if (i == 1)
            {
                netmat <- netmat[!is.na(netmat[,3]), ]
                networks[[i]] <- spMatrix(nActors, nActors, netmat[, 1],
                                          netmat[, 2], netmat[,3])
            }
            else
            {
                netmiss1 <- netmiss[[i]][, 1:2]
                storage.mode(netmiss1) <- "integer"
                networks[[i]][netmiss1[, 1:2]] <-
                    networks[[i-1]][netmiss1[, 1:2]]
            }
        }
        for (i in 1:observations)
        {
            mat1 <- networks[[i]]
            ## drop the diagonal, if present
            diag(mat1) <- 0
            mat1 <- cbind(mat1@i + 1, mat1@j + 1, mat1@x)
            ##missing edgelist
            mat2 <- netmiss[[i]]
            mat2[, 3] <- 1
            ## rows of mat1 with structural values
            struct <- mat1[, 3] %in% c(10, 11)
            ## reset real data
            mat1[struct, 3] <- mat1[struct, 3] - 10
            ## copy reset data to structural edgelist
            mat3 <- mat1[struct, , drop = FALSE]
            mat3[, 3] <- 1
            ## now remove the zeros from reset data
            mat1 <- mat1[!mat1[, 3] == 0, ]
            ## do comp change
            if (compChange)
            {
                ## revert to sparse matrices temporarily
                mat1 <- spMatrix(nrow=nActors, ncol=nActors, i = mat1[, 1],
                                 j=mat1[, 2], x=mat1[, 3])
                mat2 <- spMatrix(nrow=nActors, ncol=nActors, i = mat2[, 1],
                                 j=mat2[, 2], x=mat2[, 3])
                mat3 <- spMatrix(nrow=nActors, ncol=nActors, i = mat3[, 1],
                                 j=mat3[, 2], x=mat3[, 3])
                ones <- which(action[, i] == 1)
                twos <- which(action[, i] == 2)
                threes <- which(action[, i] == 3)
                for (j in ones) ## False data is not preceded by anything real
                {
                    if (ccOption %in% c(1, 2))
                    {
                        ## find missing values for this actor
                        use <- mat2[j, ] > 0
                        ## remove from real data (i.e. zero)
                        mat1[j, use] <- 0
                        mat1[use, j] <- 0
                        ## remove from missing data
                        mat2[j, use] <- 0
                        mat2[use, j] <- 0
                        ## remove from raw data for distances later
                        depvar[[i]][j, use] <- 0 ## zero
                        depvar[[i]][use, j] <- 0
                        depvar[[i]][j, j] <- NA
                    }
                    else if (ccOption == 3)
                    {
                        ## add the row and column to the missing data
                        mat2[j, ] <- 1
                        mat2[, j] <- 1
                        mat2[j, j] <- 0
                        ## set to missing in raw data for distances later
                        depvar[[i]][j, ] <- NA
                        depvar[[i]][, j] <- NA
                    }
                }
                for (j in threes) ## False data is preceded and followed by real
                {
                    if (ccOption %in% c(1, 2))
                    {
                        ## find missing values for this actor
                        use <- mat2[j, ] > 0
                        ## remove these from mat2, the missing data
                        mat2[j, use] <- 0
                        mat2[use, j] <- 0
                        ## carry forward
                        if (i == 1)
                        {
                            ## 0 any matches from mat1, the real data
                            mat1[j, use] <- 0
                            mat1[use, j] <- 0
                        }
                        else
                        {
                            mat1[j, use] <- networks[[i-1]][j, use]
                            mat1[use, j] <- networks[[i-1]][use, j]
                        }
                        depvar[[i]][j, use] <- 0 ##  not missing
                        depvar[[i]][use, j] <- 0
                        depvar[[i]][j, j] <- NA
                    }
                    else if (ccOption == 3)
                    {
                        ## add the row and column to the missing data
                        mat2[j, ] <- 1
                        mat2[, j] <- 1
                        mat2[j, j] <- 0
                        depvar[[i]][j, ] <- NA
                        depvar[[i]][, j] <- NA
                    }
                }
                for (j in twos) ## False data is not followed by anything real
                {
                    if (ccOption == 1)
                    {
                        ## find missing values for this actor
                        use <- mat2[j, ] > 0
                        ## remove these from mat2, the missing data
                        mat2[j, use] <- 0
                        mat2[use, j] <- 0
                        depvar[[i]][j, use] <- 0 ##  not missing
                        depvar[[i]][use, j] <- 0
                        depvar[[i]][j, j] <- NA
                        ## carry forward
                        if (i == 1)
                        {
                            ## 0 any matches from mat1, the real data
                            mat1[j, use] <- 0
                            mat1[use, j] <- 0
                        }
                        else
                        {
                            mat1[j, use] <- networks[[i-1]][j , use]
                            mat1[use, j] <- networks[[i-1]][use, j]
                        }
                    }
                    else if (ccOption %in% c(2, 3))
                    {
                        ## add the row and column to the missing data
                        mat2[j, ] <- 1
                        mat2[, j] <- 1
                        mat2[j, j] <- 0
                        depvar[[i]][j, ] <- NA
                        depvar[[i]][, j] <- NA
                    }
                }
                ## now revert to triplet matrices, after updating networks
                networks[[i]] <- mat1
                mat1 <- cbind(mat1@i + 1, mat1@j + 1, mat1@x)
                mat2 <- cbind(mat2@i + 1, mat2@j + 1, mat2@x)
                mat3 <- cbind(mat3@i + 1, mat3@j + 1, mat3@x)
                if (any (mat1[, 3] == 0) || any (mat2[, 3] == 0) ||
                    any (mat3[, 3] == 0))
                {
                    stop("zero values in sparse matrices")
                }
                if (any (duplicated(mat1[, -3])) ||
                    any (duplicated(mat2[, -3])) ||
                    any (duplicated(mat3[, -3])))
                {
                    stop("duplicate values in sparse matrices")
                }
                if (any (mat1[, 1] == mat1[, 2]) ||
                    any (mat2[, 1] == mat2[, 2]) ||
                    any (mat3[, 1] == mat3[, 2]))
                {
                    stop("loop values in sparse matrices")
                }
            }
            ##fix up storage mode to be integer
            storage.mode(mat1) <- "integer"
            storage.mode(mat2) <- "integer"
            storage.mode(mat3) <- "integer"
            ## add attribute of size
            attr(mat1,"nActors") <- nActors
            attr(mat2,"nActors") <- nActors
            attr(mat3,"nActors") <- nActors
            if (i < observations)
            {
                ## recreate the distance etc
                mymat1 <- depvar[[i]]
                mymat2 <- depvar[[i + 1]]
                ##remove structural values
                x1 <- mymat1@x
                x2 <- mymat2@x
                x1[x1 %in% c(10, 11)] <- NA
                x2[x2 %in% c(10, 11)] <- NA
                mymat1@x <- x1
                mymat2@x <- x2
                diag(mymat1) <- 0
                diag(mymat2) <- 0
                mydiff <- mymat2 - mymat1
                attr(depvar, "distance")[i] <- sum(mydiff != 0,
                                                   na.rm = TRUE)
				if (allowOnly)
				{
					if (all(mydiff@x >= 0, na.rm=TRUE))
					{
						attr(depvar, "uponly")[i] <- TRUE
					}
					if (all(mydiff@x <= 0, na.rm=TRUE))
					{
						attr(depvar, "downonly")[i] <- TRUE
					}
				}
            }
            edgeLists[[i]] <- list(mat1 = t(mat1), mat2 = t(mat2),
                                   mat3 = t(mat3))
        }
    }
    else
    {
        for (i in 1:observations) ## carry missings forward  if exist
        {
            networks[[i]] <- depvar[, , i]
            if (i == 1)
                networks[[i]][is.na(depvar[, , i])] <-0
            else ##carry missing forward!
                networks[[i]][is.na(depvar[, , i])] <-
                    networks[[i-1]][is.na(depvar[, , i])]
        }
        for (i in 1:observations)
        {
            ones <- which(action[, i] == 1)
            twos <- which(action[, i] == 2)
            threes <- which(action[, i] == 3)
            for (j in ones) ## False data is not preceded by anything real
            {
                if (ccOption %in% c(1, 2))
                {
                    use <- is.na(depvar[j, , i])
                    depvar[j, use, i] <- 0 ## not missing
                    depvar[use, j, i] <- 0
                    depvar[j, j, i] <- NA
                    networks[[i]][j, use] <- 0 ## zero
                    networks[[i]][use, j] <- 0
                }
                else if (ccOption == 3)
                {
                    depvar[j, , i] <- NA ## missing
                    depvar[, j, i] <- NA
                }
            }
            for (j in threes) ## False data is preceded and followed by real
            {

                if (ccOption %in% c(1, 2))
                {
                    use <- is.na(depvar[j, , i])
                    depvar[j, use, i] <- 0 ##  not missing
                    depvar[use, j, i] <- 0
                    depvar[j, j, i] <- NA
                    ## carry forward already done
                    if (i == 1)
                    {
                        networks[[i]][j, use] <- 0
                        networks[[i]][use, j] <- 0
                    }
                    else
                    {
                        networks[[i]][j, use] <- networks[[i-1]][j, use]
                        networks[[i]][use, j] <- networks[[i-1]][use, j]
                    }
                }
                else if (ccOption == 3)
                {
                    depvar[j, , i] <- NA ## missing
                    depvar[, j, i] <- NA
                }
            }
            for (j in twos) ## False data is not followed by anything real
            {
                if (ccOption == 1)
                {
                    use <- is.na(depvar[j, , i])
                    depvar[j, use, i] <- 0 ##  not missing
                    depvar[use, j, i] <- 0
                    depvar[j, j, i] <- NA
                    ## carry forward already done
                    if (i == 1)
                    {
                        networks[[i]][j, use] <- 0
                        networks[[i]][use, j] <- 0
                    }
                    else
                    {
                        networks[[i]][j, use] <- networks[[i-1]][j, use]
                        networks[[i]][use, j] <- networks[[i-1]][use, j]

                    }
                }
                else if (ccOption %in% c(2, 3))
                {
                    depvar[j, , i] <- NA ## missing
                    depvar[, j, i] <- NA
                }
            }
        }
        for (i in 1:observations)
        {
            if (i < observations)
            {
                ## recreate distances, as we have none in c++. (no longer true)
                mymat1 <- depvar[, ,i, drop=FALSE]
                mymat2 <- depvar[, ,i + 1, drop=FALSE]
                ##remove structural values
                mymat1[mymat1 %in% c(10, 11)] <- NA
                mymat2[mymat2 %in% c(10, 11)] <- NA
                ## and the diagonal
                diag(mymat1[, ,1]) <- 0
                diag(mymat2[, ,1]) <- 0
                mydiff <- mymat2 - mymat1
                attr(depvar, "distance")[i] <- sum(mydiff != 0,
                                                   na.rm = TRUE)
				if (allowOnly)
				{
					if (all(mydiff >= 0, na.rm=TRUE))
					{
						attr(depvar, "uponly")[i] <- TRUE
					}
					if (all(mydiff <= 0, na.rm=TRUE))
					{
						attr(depvar, "downonly")[i] <- TRUE
					}
				}
            }
            diag(networks[[i]]) <- 0
            edgeLists[[i]] <- createEdgeLists(networks[[i]], depvar[, , i],
											  FALSE)
        }
    }
    ## add attribute of nodeset
    attr(edgeLists, "nodeSet") <- attr(depvar, "nodeSet")
    ## add attribute of name
    attr(edgeLists, "name") <- attr(depvar, "name")
    ## add attribute of distance
    attr(edgeLists, "distance") <- attr(depvar, "distance")
    ## attr uponly and downonly
    attr(edgeLists, "uponly") <- attr(depvar, "uponly")
    attr(edgeLists, "downonly") <- attr(depvar, "downonly")
    ## attr symmetric
    attr(edgeLists, "symmetric") <- attr(depvar, "symmetric")
    ## attr balmean
    attr(edgeLists, "balmean") <- attr(depvar, "balmean")
    attr(edgeLists, "structmean") <- attr(depvar, "structmean")
    attr(edgeLists, "averageInDegree") <- attr(depvar, "averageInDegree")
    attr(edgeLists, "averageOutDegree") <- attr(depvar, "averageOutDegree")
	attr(edgeLists, "settings") <- attr(depvar, "settings")
    return(edgeLists = edgeLists)
}
##@unpackBipartite siena07 Reformat data for C++
unpackBipartite <- function(depvar, observations, compositionChange)
{
    edgeLists <- vector("list", observations)
    networks <- vector("list", observations)
    actorSet <- attr(depvar, "nodeSet")
    compActorSets <- sapply(compositionChange, function(x)attr(x, "nodeSet"))
    thisComp <- match(actorSet, compActorSets)
    compChange <- any(!is.na(thisComp))
    if (compChange)
    {
      #  stop("Composition change is not yet implemented for bipartite",
      #       "networks")
        action <- attr(compositionChange[[thisComp]], "action")
        ccOption <- attr(compositionChange[[thisComp]], "ccOption")
    }
    else
    {
        ccOption <- 0
        action <- matrix(0, nrow=attr(depvar, "netdims")[1], ncol=observations)
    }
    sparse <- attr(depvar, "sparse")
	allowOnly <- attr(depvar, "allowOnly")
    if (sparse)
    {
        ## require(Matrix)
        ## have a list of sparse matrices in triplet format
        ## with missings and structurals embedded and 0 based indices!
        netmiss <- vector("list", observations)
        for (i in 1:observations)
        {
            ## extract this matrix
            networks[[i]] <- depvar[[i]]
            nActors <- nrow(depvar[[i]])
            nReceivers <- ncol(depvar[[i]])
            ## stop if any duplicates
            netmat <- cbind(networks[[i]]@i+1, networks[[i]]@j+1,
                            networks[[i]]@x)
            if (any(duplicated(netmat[, 1:2])))
            {
                stop("duplicate entries in sparse matrix")
            }
            ## extract missing entries
            netmiss[[i]] <- netmat[is.na(netmat[, 3]), , drop = FALSE]
            ## carry forward missing values if any
            if (i == 1) # set missings to zero
            {
                netmat <- netmat[!is.na(netmat[,3]), ]
                networks[[i]] <- spMatrix(nActors, nReceivers, netmat[, 1],
                                          netmat[, 2], netmat[,3])
            }
            else
            {
                netmiss1 <- netmiss[[i]][, 1:2]
                storage.mode(netmiss1) <- "integer"
                networks[[i]][netmiss1[, 1:2]] <-
                    networks[[i-1]][netmiss1[, 1:2]]
            }
        }
        for (i in 1:observations)
        {
            mat1 <- networks[[i]]
            mat1 <- cbind(mat1@i + 1, mat1@j + 1, mat1@x)
            ##missing edgelist
            mat2 <- netmiss[[i]]
            mat2[, 3] <- 1
            ## rows of mat1 with structural values
            struct <- mat1[, 3] %in% c(10, 11)
            ## reset real data
            mat1[struct, 3] <- mat1[struct, 3] - 10
            ## copy reset data to structural edgelist
            mat3 <- mat1[struct, , drop = FALSE]
            ## now remove the zeros from reset data
            mat1 <- mat1[!mat1[, 3] == 0, ]
            ## do comp change
            if (compChange)
            {
                ## revert to sparse matrices temporarily
                mat1 <- spMatrix(nrow=nActors, ncol=nReceivers, i = mat1[, 1],
                                 j=mat1[, 2], x=mat1[, 3])
                mat2 <- spMatrix(nrow=nActors, ncol=nReceivers, i = mat2[, 1],
                                 j=mat2[, 2], x=mat2[, 3])
                mat3 <- spMatrix(nrow=nActors, ncol=nReceivers, i = mat3[, 1],
                                 j=mat3[, 2], x=mat3[, 3])
                ones <- which(action[, i] == 1)
                twos <- which(action[, i] == 2)
                threes <- which(action[, i] == 3)
                for (j in ones) ## False data is not preceded by anything real
                {
                    if (ccOption %in% c(1, 2))
                    {
                        ## find missing values for this actor
                        use <- mat2[j, ] > 0
                        ## remove from real data (i.e. zero)
                        mat1[j, use] <- 0
                        ## remove from missing data
                        mat2[j, use] <- 0
                        ## remove from raw data for distances later
                        depvar[[i]][j, use] <- 0 ## zero
                    }
                    else if (ccOption == 3)
                    {
                        ## add the row  to the missing data
                        mat2[j, ] <- 1
                        ## set to missing in raw data for distances later
                        depvar[[i]][j, ] <- NA
                    }
                }
                for (j in threes) ## False data is preceded and followed by real
                {
                    if (ccOption %in% c(1, 2))
                    {
                        ## find missing values for this actor
                        use <- mat2[j, ] > 0
                        ## remove these from mat2, the missing data
                        mat2[j, use] <- 0
                        ## carry forward
                        if (i == 1)
                        {
                            ## 0 any matches from mat1, the real data
                            mat1[j, use] <- 0
                        }
                        else
                        {
                            mat1[j, use] <- networks[[i-1]][j, use]
                        }
                        depvar[[i]][j, use] <- 0 ##  not missing
                    }
                    else if (ccOption == 3)
                    {
                        ## add the row to the missing data
                        mat2[j, ] <- 1
                        depvar[[i]][j, ] <- NA
                    }
                }
                for (j in twos) ## False data is not followed by anything real
                {
                    if (ccOption == 1)
                    {
                        ## find missing values for this actor
                        use <- mat2[j, ] > 0
                        ## remove these from mat2, the missing data
                        mat2[j, use] <- 0
                        depvar[[i]][j, use] <- 0 ##  not missing
                        ## carry forward
                        if (i == 1)
                        {
                            ## 0 any matches from mat1, the real data
                            mat1[j, use] <- 0
                        }
                        else
                        {
                            mat1[j, use] <- networks[[i-1]][j , use]
                        }
                    }
                    else if (ccOption %in% c(2, 3))
                    {
                        ## add the row  to the missing data
                        mat2[j, ] <- 1
                        depvar[[i]][j, ] <- NA
                    }
                }

                ## now revert to triplet matrices, after updating networks
                networks[[i]] <- mat1
                mat1 <- cbind(mat1@i + 1, mat1@j + 1, mat1@x)
                mat2 <- cbind(mat2@i + 1, mat2@j + 1, mat2@x)
                mat3 <- cbind(mat3@i + 1, mat3@j + 1, mat3@x)
                if (any (mat1[, 3] == 0) || any (mat2[, 3] == 0) ||
                    any (mat3[, 3] == 0))
                {
                    stop("zero values in sparse matrices")
                }
                if (any (duplicated(mat1[, -3])) ||
                    any (duplicated(mat2[, -3])) ||
                    any (duplicated(mat3[, -3])))
                {
                    stop("duplicate values in sparse matrices")
                }
            }
            ##fix up storage mode to be integer
            storage.mode(mat1) <- "integer"
            storage.mode(mat2) <- "integer"
            storage.mode(mat3) <- "integer"
            ## add attribute of size
            attr(mat1,"nActors") <- c(nActors, nReceivers)
            attr(mat2,"nActors") <- c(nActors, nReceivers)
            attr(mat3,"nActors") <- c(nActors, nReceivers)
            if (i < observations)
            {
                ## recreate the distance etc
                mymat1 <- depvar[[i]]
                mymat2 <- depvar[[i + 1]]
                ##remove structural values
                x1 <- mymat1@x
                x2 <- mymat2@x
                x1[x1 %in% c(10, 11)] <- NA
                x2[x2 %in% c(10, 11)] <- NA
                mymat1@x <- x1
                mymat2@x <- x2
                mydiff <- mymat2 - mymat1
                attr(depvar, "distance")[i] <- sum(mydiff != 0,
                                                   na.rm = TRUE)
				if (allowOnly)
				{
					if (all(mydiff@x >= 0, na.rm=TRUE))
					{
						attr(depvar, "uponly")[i] <- TRUE
					}
					if (all(mydiff@x <= 0, na.rm=TRUE))
					{
						attr(depvar, "downonly")[i] <- TRUE
					}
				}
            }
            edgeLists[[i]] <- list(mat1 = t(mat1), mat2 = t(mat2),
                                   mat3 = t(mat3))
        }
    }
    else
    {
        for (i in 1:observations) ## carry missings forward  if exist
        {
            networks[[i]] <- depvar[, , i]
            if (i == 1)
                networks[[i]][is.na(depvar[, , i])] <-0
            else ##carry missing forward!
                networks[[i]][is.na(depvar[, , i])] <-
                    networks[[i-1]][is.na(depvar[, , i])]
        }
        for (i in 1:observations)
        {
            ones <- which(action[, i] == 1)
            twos <- which(action[, i] == 2)
            threes <- which(action[, i] == 3)
            for (j in ones) ## False data is not preceded by anything real
            {
                if (ccOption %in% c(1, 2))
                {
                    use <- is.na(depvar[j, , i])
                    depvar[j, use, i] <- 0 ## not missing
                    networks[[i]][j, use] <- 0 ## zero
                }
                else if (ccOption == 3)
                {
                    depvar[j, , i] <- NA ## missing
                }
            }
            for (j in threes) ## False data is preceded and followed by real
            {

                if (ccOption %in% c(1, 2))
                {
                    use <- is.na(depvar[j, , i])
                    depvar[j, use, i] <- 0 ##  not missing
                    ## carry forward already done
                    if (i == 1)
                    {
                        networks[[i]][j, use] <- 0
                    }
                    else
                    {
                        networks[[i]][j, use] <- networks[[i-1]][j, use]
                    }
                }
                else if (ccOption == 3)
                {
                    depvar[j, , i] <- NA ## missing
                }
            }
            for (j in twos) ## False data is not followed by anything real
            {
                if (ccOption == 1)
                {
                    use <- is.na(depvar[j, , i])
                    depvar[j, use, i] <- 0 ##  not missing
                    ## carry forward already done
                    if (i == 1)
                    {
                        networks[[i]][j, use] <- 0
                    }
                    else
                    {
                        networks[[i]][j, use] <- networks[[i-1]][j, use]

                    }
                }
                else if (ccOption %in% c(2, 3))
                {
                    depvar[j, , i] <- NA ## missing
                }
            }
        }
        for (i in 1:observations)
        {
            if (i < observations)
            {
                ## recreate distances, as we have none in c++. (no longer true)
                mymat1 <- depvar[,,i, drop=FALSE]
                mymat2 <- depvar[,,i + 1,drop=FALSE]
                ##remove structural values
                mymat1[mymat1 %in% c(10,11)] <- NA
                mymat2[mymat2 %in% c(10,11)] <- NA
                mydiff <- mymat2 - mymat1
                attr(depvar, "distance")[i] <- sum(mydiff != 0,
                                                         na.rm = TRUE)
				if (allowOnly)
				{
					if (all(mydiff >= 0, na.rm=TRUE))
					{
						attr(depvar, "uponly")[i] <- TRUE
					}
					if (all(mydiff <= 0, na.rm=TRUE))
					{
						attr(depvar, "downonly")[i] <- TRUE
					}
				}
            }

            edgeLists[[i]] <- createEdgeLists(networks[[i]], depvar[, , i], TRUE)
        }
    }
    ## add attribute of nodeset
    attr(edgeLists, "nodeSet") <- attr(depvar, "nodeSet")
    ## add attribute of name
    attr(edgeLists, "name") <- attr(depvar, "name")
    ## add attribute of distance
    attr(edgeLists, "distance") <- attr(depvar, "distance")
    ## attr uponly and downonly
    attr(edgeLists, "uponly") <- attr(depvar, "uponly")
    attr(edgeLists, "downonly") <- attr(depvar, "downonly")
    ## attr symmetric
    attr(edgeLists, "symmetric") <- attr(depvar, "symmetric")
    ## attr balmean
    attr(edgeLists, "balmean") <- attr(depvar, "balmean")
    ## attr structmean
    attr(edgeLists, "structmean") <- attr(depvar, "structmean")
    attr(edgeLists, "averageOutDegree") <- attr(depvar, "averageOutDegree")
    return(edgeLists = edgeLists)
}
##@unpackBehavior siena07 Reformat data for C++
unpackBehavior<- function(depvar, observations)
{
    beh <- depvar[, 1, ]
    behmiss <- is.na(beh)
    allna <- apply(beh, 1, function(x)all(is.na(x)))
    modes <- attr(depvar, "modes")
    ## carry forward missings ### nb otherwise use the mode
    ## allNAs: use modes
    beh[allna, ] <- rep(modes, each=sum(allna))
    for (i in 2:observations)
    {
            beh[is.na(beh[, i]), i] <-  beh[is.na(beh[, i]), i - 1]
    }
    for (i in (observations-1):1)
    {
            beh[is.na(beh[, i]), i] <-  beh[is.na(beh[, i]), i +1]
    }
    ## need a better definition of structural for behavior variables
    ## struct <- beh %in% c(10, 11)
    ## beh[struct] <- beh[struct] - 10
    ## behstruct <- beh
    ## behstruct[!struct] <- 0

    ## add attribute of nodeset
    attr(beh, "nodeSet") <- attr(depvar, "nodeSet")
    ## add attribute of name
    attr(beh, "name") <- attr(depvar, "name")
    ## attr uponly and downonly
    attr(beh, "uponly") <- attr(depvar, "uponly")
    attr(beh, "downonly") <- attr(depvar, "downonly")
    ## attr symmetric
    attr(beh, "symmetric") <- attr(depvar, "symmetric")
    ## attr distance
    attr(beh, "distance") <- attr(depvar, "distance")
    ## attr simMean
    attr(beh, "simMean") <- attr(depvar, "simMean")
    ## attr simMeans
    attr(beh, "simMeans") <- attr(depvar, "simMeans")
    storage.mode(beh) <- "integer"
    list(beh=beh, behmiss=behmiss)
}
##@convertToStructuralZeros Miscellaneous To be implemented
convertToStructuralZeros <- function()
{
}

##@unpackCDyad siena07 Reformat data for C++
unpackCDyad<- function(dycCovar)
{
    sparse <- attr(dycCovar, "sparse")
	bipartite <- attr(dycCovar, "type") == "bipartite"
    if (sparse)
    {
        ## have a list containing 1 sparse matrix in triplet format
        ## with missings embedded
        ## with 0 based indices!
        varmat <- cbind(dycCovar[[1]]@i+1, dycCovar[[1]]@j+1, dycCovar[[1]]@x)
        if (any(duplicated(varmat[, 1:2])))
        {
            stop("duplicate entries in sparse matrix dyadic covariate")
        }
        ##drop the diagonal, if present - not for bipartite
        if (!bipartite)
        {
            varmat <- varmat[varmat[,1] != varmat[, 2],]
        }
        mat1 <- varmat
        mat1[is.na(varmat[, 3]), 3] <- attr(dycCovar, "mean")
        mat1 <- mat1[!mat1[, 3] == 0, ]
        ## add attribute of dim
        attr(mat1, "nActors1") <- nrow(dycCovar[[1]])
        attr(mat1, "nActors2") <- ncol(dycCovar[[1]])
        mat2 <- varmat[is.na(varmat[, 3]), , drop=FALSE]
        mat2[, 3] <- 1
        ## add attribute of dim
        attr(mat2,"nActors1") <- nrow(dycCovar[[1]])
        attr(mat2,"nActors2") <- ncol(dycCovar[[1]])
        edgeLists <-  list(t(mat1), t(mat2))
    }
    else
    {
        if (!bipartite)
        {
            diag(dycCovar) <- 0
        }
        dycCovar1 <- dycCovar
        dycCovar1[is.na(dycCovar1)] <- attr(dycCovar, "mean")
        edgeLists <- createCovarEdgeList(dycCovar1, dycCovar)
    }
    ## add attribute of nodesets
    attr(edgeLists, "nodeSet") <- attr(dycCovar, "nodeSet")
    ## add attribute of type
    attr(edgeLists, "type") <- attr(dycCovar, "type")
    ## add attribute of name
    attr(edgeLists, "name") <- attr(dycCovar, "name")
    ## add attribute of mean
    attr(edgeLists, "mean") <- attr(dycCovar, "mean")
    return(edgeLists = edgeLists)
}


##@unpackVDyad siena07 Reformat data for C++
unpackVDyad<- function(dyvCovar, observations)
{
    edgeLists <- vector("list", observations)
    sparse <- attr(dyvCovar, "sparse")
    means <- attr(dyvCovar, "meanp")
	bipartite <- attr(dyvCovar, "type") == "bipartite"
    if (sparse)
    {
        ## have a list of sparse matrices in triplet format
        ## with 0 based indices!
        for (i in 1:(observations - 1))
        {
            thisvar <- dyvCovar[[i]]
            varmat <- cbind(thisvar@i+1, thisvar@j+1, thisvar@x)
            ## drop the diagonal, if present no - bipartite?
            if (!bipartite)
            {
                varmat <- varmat[varmat[,1] != varmat[, 2],]
            }
            mat1 <- varmat
            mat1[is.na(varmat[, 3]), 3] <- means[i]
            mat1 <- mat1[!mat1[, 3] == 0, ]
            mat2 <- varmat[is.na(varmat[, 3]),, drop=FALSE ]
            mat2[, 3] <- 1
            ## add attribute of size
            attr(mat1, "nActors1") <- nrow(dyvCovar[[i]])
            attr(mat1, "nActors2") <- ncol(dyvCovar[[i]])
            attr(mat2, "nActors1") <- nrow(dyvCovar[[i]])
            attr(mat2, "nActors2") <- ncol(dyvCovar[[i]])
            edgeLists[[i]] <- list(t(mat1), t(mat2))
        }
    }
    else
    {
        for (i in 1:(observations - 1))
        {
            if (!bipartite)
            {
                diag(dyvCovar[, , i]) <- 0
            }
            thisvar <- dyvCovar[, , i]
            thisvar[is.na(thisvar) ] <- means[i]
            edgeLists[[i]] <- createCovarEdgeList(thisvar, dyvCovar[, , i])
        }
    }
    ## add attribute of nodeset
    attr(edgeLists, "nodeSet") <- attr(dyvCovar, "nodeSet")
    ## add attribute of type
    attr(edgeLists, "type") <- attr(dyvCovar, "type")
    ## add attribute of name
    attr(edgeLists, "name") <- attr(dyvCovar, "name")
    ## add attribute of mean
    attr(edgeLists, "mean") <- attr(dyvCovar, "mean")
    return(edgeLists = edgeLists)
}

##@unpackData siena07 Reformat data for C++
unpackData <- function(data, x)
{
    f <- NULL
    observations<- data$observations
    types <- sapply(data$depvars, function(x) attr(x, "type"))
    f$nDepvars <- length(data$depvars)
    oneModes <- data$depvars[types == "oneMode"]
    Behaviors <- data$depvars[types == "behavior"]
    bipartites <- data$depvars[types == "bipartite"]
	## add the settings
	oneModes <- lapply(oneModes, function(depvar)
	   {
		   name <- attr(depvar, "name")
		   if (name %in% names(x$settings))
		   {
			   attr(depvar, "settings") <- c("universal", "primary",
											 x$settings[[name]])
		   }
		   depvar
	   }
		   )
    f$nets <- lapply(oneModes, function(x, n, comp)
					 unpackOneMode(x, n, comp),
                     n = observations, comp=data$compositionChange)
    names(f$nets) <- names(oneModes)
    f$bipartites <- lapply(bipartites, function(x, n, comp)
                           unpackBipartite(x, n, comp),
                     n = observations, comp=data$compositionChange)
    names(f$bipartites) <- names(bipartites)
    f$behavs <-  lapply(Behaviors, function(x, n) unpackBehavior(x, n),
                        n = observations)
    names(f$behavs) <- names(Behaviors)
    f$observations <- observations
    f$seed<- vector("list", observations - 1)
    f$depvars <- data$depvars
    f$nodeSets <- data$nodeSets
    f$oneModes <- oneModes
    f$Behaviors <- Behaviors
    f$oneModeUpOnly <- sapply(oneModes, function(x) attr(x, "uponly"))
    f$oneModeDownOnly <- sapply(oneModes, function(x) attr(x, "downonly"))
    f$behaviorUpOnly <- sapply(Behaviors, function(x) attr(x, "uponly"))
    f$behaviorDownOnly <- sapply(Behaviors, function(x) attr(x,
                                                             "downonly"))
    f$distances <- sapply(data$depvars, function(x) attr(x, "distance"))
    f$cCovars <- data$cCovars
    f$vCovars <- data$vCovars
    ## dyadic covars need to be edgelists
    f$dycCovars <- lapply(data$dycCovars, function(x) unpackCDyad(x))
    f$dyvCovars <- lapply(data$dyvCovars, function(x,n) unpackVDyad(x,n),
                          n=observations)
    ## create the composition change event lists
    f$exog <- lapply(data$compositionChange, function(x)
                     unpackCompositionChange(x))
    f
}

##@unpackCompositionChange siena07 Reformat data for C++
unpackCompositionChange <- function(compositionChange)
{
    atts <- attributes(compositionChange)
    events <- atts$events
    activeStart <- atts$activeStart
    observations <- ncol(activeStart)
    ## check that there is someone there always
    for (i in 1:(observations - 1))
    {
        activeAll <- sum(activeStart[, i] & activeStart[, i + 1])
        if (activeAll < 2)
        {
            active <- sum(activeStart[, i])
            if (active == 1)
                stop("Only one active actor at start of period", i)
            else if (active == 0)
                stop("No active actors at start of period", i)
            perEvents <- events[events$period == i,]
            perEvents <- perEvents[order(perEvents$time),]
            changes <- c(1, -1)[as.numeric(as.character(perEvents$event))]
            active <- active + cumsum(changes)
            if (any(active < 2))
            {
                stop("No/only one active actor(s) left.")
            }
        }
    }
    events <- events[events$time > 1e-10,]
    exog <- list(events=events, activeStart=activeStart)
    attr(exog, "nodeSet") <- attr(compositionChange, "nodeSet")
    exog
}
##@fixUpEffectNames siena07 Replace # and construct interaction names
fixUpEffectNames <- function(effects)
{
    ## replace # by the parm value in function and effect names
    effects$effectName <-
        sapply(1:nrow(effects), function(x, y)
           {
               y <- y[x, ]
               gsub("#", y$parm, y$effectName)
           }, effects)
    effects$functionName <-
        sapply(1:nrow(effects), function(x, y)
           {
               y <- y[x, ]
               gsub("#", y$parm, y$functionName)
           }, y=effects)

    ##validate user-specified network interactions
    interactions <- effects[effects$shortName == "unspInt" & effects$include &
                            effects$effect1 > 0, ]
    if (nrow(interactions) > 0)
    {
        unspIntNames <-
            sapply(1:nrow(interactions), function(x, y, z)
               {
                   y <- y[x, ] ## get the interaction effect
                   twoway <- y$effect3 == 0
                   ## now get the rows which are to interact
                   inter1 <- z[z$effectNumber == y$effect1, ]
                   if (nrow(inter1) != 1 )
                   {
                       stop("invalid network interaction specification: ",
                            "effect number 1")
                   }
                   inter2 <- z[z$effectNumber == y$effect2, ]
                   if (nrow(inter2) != 1 )
                   {
                       stop("invalid network interaction specification: ",
                            "effect number 2")
                   }
                   if (!twoway)
                   {
                       inter3 <- z[z$effectNumber == y$effect3, ]
                       if (nrow(inter3) != 1)
                       {
                           stop("invalid network interaction specification: ",
                                "effect number 3")
                       }
                   }
                   else
                   {
                       inter3 <- z[is.na(z$effectNumber), ]
                       ## should be empty row
                   }
                   if (twoway)
                   {
                       if (inter1$name != inter2$name)
                       {
                           stop("invalid network interaction specification: ",
                                "must all be same network")
                       }
                       if (inter1$type != inter2$type)
                       {
                           stop("invalid network interaction specification: ",
                                "must all be same type: ",
                                "evaluation, endowment or creation")
                       }
                   }
                   else
                   {
                       if (inter1$name != inter2$name ||
                           inter1$name != inter3$name)
                       {
                           stop("invalid network interaction specification: ",
                                "must all be same network")
                       }
                       if (inter1$type != inter2$type ||
                           inter1$type != inter3$type)
                       {
                           stop("invalid network interaction specification: ",
                                "must all be ",
                                "same type: evaluation, endowment or creation ")
                       }
                   }
                   ## check types
                   inters <- rbind(inter1, inter2, inter3)
                   egos <- which(inters$interactionType == "ego")
                   egoCount <- length(egos)
                   dyads <- which(inters$interactionType == "dyadic")
                   dyadCount <- length(dyads)
                   if (twoway)
                   {
                       if (egoCount < 1 && dyadCount != 2)
                       {
                           stop("invalid network interaction specification: ",
                                "must be at least one ego or both dyadic ",
                                "effects")
                       }
                   }
                   else
                   {
                       if (egoCount < 2 && dyadCount != 3)
                       {
                           stop("invalid network interaction specification: ",
                                "must be at least two ego or all dyadic ",
                                "effects")
                       }
                   }
                   ## construct a name
                   ## make sure the egos are at the front of inters
                   if (egoCount > 0)
                   {
                       inters <- rbind(inters[egos, ], inters[-egos, ])
                   }
 				   tmpnames <- inters$effectName
				   tmpnames[-1] <- sub(paste(inters$name[1], ": ",
											 sep=""), "", tmpnames[-1])
				   tmpname <- paste(tmpnames, collapse = " x ")
                   if (twoway && nchar(tmpname) < 38)
                   {
                       tmpname <- paste("int. ", tmpname)
                   }
                   if (!twoway)
                   {
                       tmpname <- paste("i3.", tmpname)
                   }
                   tmpname
               }, y=interactions, z=effects)
        effects[effects$shortName == "unspInt" & effects$include &
                !is.na(effects$effect1), c("effectName", "functionName")] <-
                    unspIntNames
    }
    ##validate user-specified behavior interactions
    interactions <- effects[effects$shortName == "behUnspInt" &
                            effects$include &
                            effects$effect1 > 0, ]
    if (nrow(interactions) > 0)
    {
        unspIntNames <-
            sapply(1:nrow(interactions), function(x, y, z)
               {
                   y <- y[x, ] ## get the interaction effect
                   twoway <- y$effect3 == 0
                   ## now get the rows which are to interact
                   inter1 <- z[z$effectNumber == y$effect1, ]
                   if (nrow(inter1) != 1 )
                   {
                       stop("invalid behavior interaction specification: ",
                            "effect number 1")
                   }
                   inter2 <- z[z$effectNumber == y$effect2, ]
                   if (nrow(inter2) != 1 )
                   {
                       stop("invalid behavior interaction specification: ",
                            "effect number 2")
                   }
                   if (!twoway)
                   {
                       inter3 <- z[z$effectNumber == y$effect3, ]
                       if (nrow(inter3) != 1)
                       {
                           stop("invalid behavior interaction specification: ",
                                "effect number 3")
                       }
                   }
                   else
                   {
                       inter3 <- z[is.na(z$effectNumber), ]
                       ## should be empty row
                   }
                   if (twoway)
                   {
                       if (inter1$name != inter2$name)
                       {
                           stop("invalid behavior interaction specification: ",
                                "must all be same behavior variable")
                       }
                       if (inter1$type != inter2$type)
                       {
                           stop("invalid behavior interaction specification: ",
                                "must be same type: evaluation, endowment ",
                                "or creation")
                       }
                   }
                   else
                   {
                       if (inter1$name != inter2$name ||
                           inter1$name != inter3$name)
                       {
                           stop("invalid behavior interaction specification: ",
                                "must all be same behavior variable")
                       }
                       if (inter1$type != inter2$type ||
                           inter1$type != inter3$type)
                       {
                           stop("invalid behavior interaction specification: ",
                                "must all be ",
                                "same type: evaluation, endowment or creation")
                       }
                   }
                   ## check types - at most one should be not OK here
                   inters <- rbind(inter1, inter2, inter3)
                   if (length(which(inters$interactionType != "OK")) > 1)
				   {
				   	   stop("invalid behavior interaction specification: ",
				   			"at most one effect with interactionType ",
				   			"not OK is allowed")
                   }
                   ##if (any(inters$interactionType != "OK"))
                   ##{
                   ##    stop("invalid behavior interaction specification: ",
                   ##         "only effects with interactionType OK are allowed")
                   ##}
                   ## construct a name
				   tmpnames <- inters$effectName
				   tmpnames[-1] <- sub(paste("behavior ", inters$name[1], " ",
											 sep=""), "", tmpnames[-1])
				   tmpname <- paste(tmpnames, collapse = " x ")
                   if (twoway && nchar(tmpname) < 38)
                   {
                       tmpname <- paste("int. ", tmpname)
                   }
                   if (!twoway)
                   {
                       tmpname <- paste("i3.", tmpname)
                   }
                   tmpname
               }, y=interactions, z=effects)
        effects[effects$shortName == "behUnspInt" & effects$include &
                !is.na(effects$effect1), c("effectName", "functionName")] <-
                    unspIntNames
    }
    effects
}

##@updateTheta siena07 Copy theta values from previous fit
updateTheta <- function(effects, prevAns)
{
	if (!inherits(effects, "data.frame"))
	{
		stop("effects is not a data.frame")
	}
	if (!inherits(prevAns, "sienaFit"))
	{
		stop("prevAns is not an RSiena fit object")
	}
	prevEffects <- prevAns$requestedEffects
	prevEffects$initialValue <- prevAns$theta
	if (prevAns$cconditional)
	{
		condEffects <- attr(prevAns$f, "condEffects")
		condEffects$initialValue <- prevAns$rate
		prevEffects <- rbind(prevEffects, condEffects)
	}
	oldlist <- apply(prevEffects, 1, function(x)
					 paste(x[c("name", "shortName",
							   "type", "groupName",
							   "interaction1", "interaction2",
							   "period")],
						   collapse="|"))
	efflist <- apply(effects, 1, function(x)
					 paste(x[c("name", "shortName",
							   "type", "groupName",
							   "interaction1", "interaction2",
							   "period")],
						   collapse="|"))
	use <- efflist %in% oldlist
	effects$initialValue[use] <-
		prevEffects$initialValue[match(efflist, oldlist)][use]
	effects
}

##@addSettingseffects siena07 add extra rate effects for settings model
addSettingsEffects <- function(effects, x)
{
	## need a list of settings (constant dyadic covariate name) for some or all
	## dependent networks. If symmetric this is equivalent to a forcing model.
	## The covariate actor sets should match the network actor sets.
	## ?? would it ever make sense for bipartites? Allow it here for now and see!
	## maybe better to add the settings pars to the data object but for now
	## they are on the model with maxdegree. TODO write some validation
	varNames <- names(x$settings)
	nbrSettings <- sapply(x$settings, length)
	tmp <-
		lapply(1:length(x$settings), function(i)
		   {
			   ## find effects for this variable
			   varEffects <- effects[effects$name == varNames[i], ]
			   ## find relevant rate effects
			   dupl <- varEffects[varEffects$basicRate, ]
			   ## make extra copies
			   newEffects <- dupl[rep(1:nrow(dupl),
									  each = nbrSettings[i] + 2), ]
			   newEffects <- split(newEffects, list(newEffects$group,
								   newEffects$period))
			   newEffects <- lapply(newEffects, function(dd)
				  {
					  dd$setting <- c("universal", "primary",
											  x$settings[[i]])
					  i1 <- regexpr("rate", dd$effectName)
					  dd$effectName <-
						  paste(substr(dd$effectName, 1, i1 - 2),
											 dd$setting,
								substring(dd$effectName, i1))
					  dd
				  }
					  )
			   newEffects <- do.call(rbind, newEffects)
			   ## add extra column to non basic rate effects
			   varEffects$setting <- rep("", nrow(varEffects))
			   ## combine them
			   rbind(newEffects, varEffects[!varEffects$basicRate, ])
		   })
	## join them together again
	do.call(rbind, tmp)
}

