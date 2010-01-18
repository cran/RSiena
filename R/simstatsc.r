#/******************************************************************************
# * SIENA: Simulation Investigation for Empirical Network Analysis
# *
# * Web: http://www.stats.ox.ac.uk/~snidjers/siena
# *
# * File: simstatsc.r
# *
# * Description: This module contains the code for simulating the process,
# * communicating with C++.
# *****************************************************************************/
##@simstats0c siena07 Simulation Module
simstats0c <-function(z, x, INIT=FALSE, TERM=FALSE, initC=FALSE, data=NULL,
                        effects=NULL, fromFiniteDiff=FALSE,
                      profileData=FALSE, prevAns=NULL, returnDeps=FALSE)
{
    if (INIT || initC)  ## initC is to initialise multiple C processes in phase3
    {
        if (!initC)
        {
            if (!inherits(data,'siena'))
                stop('not valid siena data object')
            defaultEffects <- getEffects(data)
            if (is.null(effects))
                effects <- defaultEffects
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
                    stop("invalid effect requested: ", userlist[bad])
                }
            }
            if (!inherits(effects, 'data.frame'))
                stop('effects is not a data.frame')
            if (x$useStdInits)
            {
                if (any(effects$effectName != defaultEffects$effectName))
                {
                    stop('Cannot use standard initialisation with a ',
                          'different effect list')
                }
                effects$initialValue <- defaultEffects$initialValue
            }
            ## find any effects not included which are needed for interactions
            interactionNos <- unique(c(effects$effect1, effects$effect2,
                                       effects$effect3))
            interactionNos <- interactionNos[interactionNos > 0]
            interactionMainEffects <- effects[interactionNos, ]
            effects$requested <- effects$include
            requestedEffects <- effects[effects$include, ]

            effects$include[interactionNos] <- TRUE
            effects <- effects[effects$include,]

            ## split and rejoin both versions before continuing
            effects1 <- split(requestedEffects, requestedEffects$name)
            if (inherits(data, "sienaGroup"))
                depvarnames <- names(data[[1]]$depvars)
            else
                depvarnames <- names(data$depvars)
            effects1order <- match(depvarnames, names(effects1))
            requestedEffects <- do.call(rbind, effects1[effects1order])
            row.names(requestedEffects) <- 1:nrow(requestedEffects)
            effects1 <- split(effects, effects$name)
            effects1order <- match(depvarnames, names(effects1))
            effects <- do.call(rbind, effects1[effects1order])
            row.names(effects) <- 1:nrow(effects)
            z$theta <- requestedEffects$initialValue
            z$fixed <- requestedEffects$fix
            z$test <- requestedEffects$test
            z$pp <- length(z$test)
            z$posj <- rep(FALSE,z$pp)
            z$posj[requestedEffects$basicRate] <- TRUE
            z$BasicRateFunction <- z$posj
            effects <- fixUpEffectNames(effects)

            ## copy interaction names to the requested effects
            requestedEffects$effectName <- effects[effects$requested,
                                                   "effectName"]
            requestedEffects$functionName <- effects[effects$requested,
                                                   "functionName"]

            if (inherits(data, 'sienaGroup'))
            {
                nGroup <- length(data)
            }
            else
            {
                nGroup <- 1
                data <- sienaGroupCreate(list(data), singleOK=TRUE)
            }
            if (is.na(x$cconditional))
            {
                x$cconditional <- length(depvarnames) == 1
                if (x$cconditional)
                {
                    x$condvarno <- 1
                }
            }
            z$cconditional <- FALSE
            if (x$cconditional)
            {
                types <- sapply(data[[1]]$depvars, function(x) attr(x, 'type'))
                nets <- sum(types != "behavior")
                if (nets == 1)
                {
                    z$cconditional <- TRUE
                    ## find the conditioning variable
                    observations <- attr(data, 'observations')
                    if (x$condname != '')
                    {
                        z$condvarno <- match(x$condname, attr(data, "netnames"))
                        z$condname <- x$condname
                    }
                    else
                    {
                        z$condvarno <- x$condvarno
                        z$condname <- attr(data, 'netnames')[x$condvarno]
                    }
                    z$condtype <- attr(data, "types")[z$condvarno]
                    if (z$condtype == 'oneMode')
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
                    z$pp<- z$pp-length(z$condvar)
                    z$scale<- z$scale[-z$condvar]
                    z$BasicRateFunction <- z$posj[-z$condvar]
                    z$posj <- z$posj[-z$condvar]
                    z$theta[z$posj] <-
                        z$theta[z$posj] /
                            requestedEffects$initialValue[z$condvar]
                    z$ntim<- matrix(NA, nrow=x$n3, ncol=observations)
                }
            }
            ## unpack data and put onto f anything we may need next time round.
            f <- lapply(data, function(x) unpackData(x))
            attr(f, "netnames") <- attr(data, "netnames")
            attr(f, "symmetric") <- attr(data, "symmetric")
            attr(f, "allUpOnly") <- attr(data, "allUpOnly")
            attr(f, "allDownOnly") <- attr(data, "allDownOnly")
            attr(f, "anyUpOnly") <- attr(data, "anyUpOnly")
            attr(f, "anyDownOnly") <- attr(data, "anyDownOnly")
            attr(f, "types") <- attr(data, "types")
            attr(f, "observations") <- attr(data, "observations")
            attr(f, "compositionChange") <- attr(data, "compositionChange")
            attr(f, "exooptions") <- attr(data, "exooptions")
           ## if any networks symmetric must use finite differences
            syms <- attr(data,"symmetric")
            z$FinDiffBecauseSymmetric <- FALSE
            if (any(!is.na(syms) & syms))
            {
                z$FinDiff.method <- TRUE
                z$FinDiffBecauseSymmetric <- TRUE
            }
        	if (z$cconditional)
            {
                attr(f, "change") <-
                    sapply(f, function(xx)attr(xx$depvars[[z$condname]],
                                               'distance'))
                attr(f,"condEffects") <- requestedEffects[z$condvar,]
                effcondvar <-
                    (1:nrow(effects))[effects$name==
                                      z$condname][1:observations]
                effects <- effects[-effcondvar, ]
                requestedEffects <- requestedEffects[-z$condvar,]
            }
            ## see if we can use the original dfra
            if (!is.null(prevAns) && inherits(prevAns, "sienaFit"))
            {
                if (all(rownames(prevAns$dfra) == requestedEffects$shortName)
                    && !is.null(prevAns$sf))
                {
                    z$haveDfra <- TRUE
                    z$dfra <- prevAns$dfra
                    z$sf <- prevAns$sf
                    ## use thetas too, unless use standard values
                    if (!x$useStdInits)
                    {
                        requestedEffects$initialValue <- prevAns$theta
                        if (!is.null(prevAns$condvar))
                        {
                            ## z$condvar has the subscripts of included
                            ## parameters
                            ## that correspond to the conditional variable
                            ## need to scale the other rates again
                            requestedEffects$initialValue[z$posj] <-
                                requestedEffects$initialValue[z$posj] /
                                    prevAns$rate
                        }
                        z$theta <- requestedEffects$initialValue
                    }
                }
            }
            z$effects <- effects
            z$requestedEffects <- requestedEffects
        }
        else
        {
            f <- FRANstore()
            ## Would like f to be just the data objects plus the attributes
            ## but need the effects later. Also returnDeps flag
            ff <- f
            f$pData <- NULL
            f$pModel <-  NULL
            f$myeffects <-  NULL
            f$observations <-  NULL
            f$randomseed2 <- NULL
            f$seeds <- NULL
            f$returnDeps <- NULL
            f$depNames <- NULL
            f$groupNames <- NULL
            f$nGroup <- NULL
       }
        ##browser()
        pData <- .Call('setupData', PACKAGE="RSiena",
                       lapply(f, function(x)(as.integer(x$observations))),
                       lapply(f, function(x)(x$nodeSets)))
        ans <- .Call('OneMode', PACKAGE="RSiena",
                    pData, lapply(f, function(x)x$nets))
        ans <- .Call('Bipartite', PACKAGE="RSiena",
                    pData, lapply(f, function(x)x$bipartites))
        ans <- .Call('Behavior', PACKAGE="RSiena",
                     pData, lapply(f, function(x)x$behavs))
        ans <-.Call('ConstantCovariates', PACKAGE="RSiena",
                   pData, lapply(f, function(x)x$cCovars))
        ans <-.Call('ChangingCovariates', PACKAGE="RSiena",
                   pData, lapply(f, function(x)x$vCovars))
        ans <-.Call('DyadicCovariates', PACKAGE="RSiena",
                   pData, lapply(f, function(x)x$dycCovars))
        ans <-.Call('ChangingDyadicCovariates', PACKAGE="RSiena",
                   pData, lapply(f, function(x)x$dyvCovars))
        ans <-.Call('ExogEvent', PACKAGE="RSiena",
                   pData, lapply(f, function(x)x$exog))
        ##store the address
        f$pData <- pData
        ## register a finalizer
        ans <- reg.finalizer(f$pData, clearData, onexit = FALSE)
        if (!initC)
        {
            storage.mode(effects$parm) <- 'integer'
            storage.mode(effects$group) <- 'integer'
            storage.mode(effects$period) <- 'integer'
            effects$effectPtr <- NA
            splitFactor <- factor(effects$name, levels=attr(f, "netnames"))
            myeffects <- split(effects, splitFactor)
        }
        else
        {
            myeffects <- ff$myeffects
            returnDeps <- ff$returnDeps
            nGroup <- ff$nGroup
        }
        ## remove interaction effects and save till later
        basicEffects <- lapply(myeffects, function(x)
                        {
                            x[!x$shortName %in% c("unspInt", "behUnspInt"), ]
                        }
                            )
        interactionEffects <- lapply(myeffects, function(x)
                        {
                            x[x$shortName %in% c("unspInt", "behUnspInt"), ]
                        }
                            )
        ans <- .Call('effects', PACKAGE="RSiena",
                    pData, basicEffects)
        pModel <- ans[[1]][[1]]
       ## browser()
        for (i in 1:length(ans[[2]])) ## ans[[2]] is a list of lists of
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
        ans <- .Call('interactionEffects', PACKAGE="RSiena",
                     pData, pModel, interactionEffects)
        ## copy these pointers to the interaction effects and then rejoin
        for (i in 1:length(ans[[1]])) ## ans is a list of lists of
            ## pointers to effects. Each list corresponds to one
            ## dependent variable
        {
            if (nrow(interactionEffects[[i]]) > 0)
            {
                effectPtr <- ans[[1]][[i]]
                interactionEffects[[i]]$effectPtr <- effectPtr
            }
            myeffects[[i]] <- rbind(basicEffects[[i]], interactionEffects[[i]])
        }
        ## remove the effects only created as underlying effects
        ## for interaction effects
        myeffects <- lapply(myeffects, function(x)
                        {
                            x[x$requested, ]
                        }
                            )
        if (!initC)
        {
            ans <- .Call('getTargets', PACKAGE="RSiena",
                         pData, pModel, myeffects)
            z$targets <- rowSums(ans)
            z$targets2 <- ans
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
        }
        if (z$cconditional)
        {
            CONDVAR <- z$condname
            CONDTARGET <- attr(f, "change")
         ##   cat(CONDTARGET, '\n')
        }
        else
        {
            CONDVAR <- NULL
            CONDTARGET <- NULL
        }
        ans <- .Call("setupModelOptions", PACKAGE="RSiena",
                     pData, pModel, MAXDEGREE, CONDVAR, CONDTARGET,
                     profileData, z$parallelTesting)
        f$myeffects <- myeffects
        if (!initC)
        {
            DataReport(z, x, f)
            f$randomseed2 <- z$randomseed2
        }
        else
        {
            f$randomseed2 <- ff$randomseed2
        }
        f$observations <- attr(f, "observations") + 1
        f$returnDeps <- returnDeps
        f$depNames <- names(f[[1]]$depvars)
        f$groupNames <- names(f)[1:nGroup]
        f$nGroup <- nGroup
        if (!initC)
        {
            z$f <- f
        }
        if (initC || z$int == 1)
        {
            f[1:nGroup] <- NULL
        }
        FRANstore(f) ## store f in FRANstore
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
        if (z$cconditional)
        {
            z$rate<- colMeans(z$ntim, na.rm=TRUE)
            z$vrate <- apply(z$ntim, 2, sd, na.rm=TRUE)
            z$theta[z$posj] <- z$theta[z$posj] * z$rate
            z$covtheta[z$posj, ] <- z$covtheta[z$posj, ] * z$rate
            z$covtheta[, z$posj] <- z$covtheta[,z$posj ] * z$rate
        }
        f <- FRANstore()
        f$pModel <- NULL
        f$pData <- NULL
        FRANstore(NULL) ## clear the stored object
        PrintReport(z, x)
        if (sum(z$test))
        {
            z$fra <- colMeans(z$sf, na.rm=TRUE)
            ans <- ScoreTest(z$pp, z$dfra, z$msf, z$fra, z$test, x$maxlike)
            z <- c(z, ans)
            TestOutput(z, x)
        }
        dimnames(z$dfra)[[1]] <- as.list(z$requestedEffects$shortName)
        return(z)
    }
    ## iteration entry point
    f <- FRANstore()
   # browser()
   # cat(f$randomseed2, f$storedseed, '\n')
    if (fromFiniteDiff || z$Phase == 2)
    {
        returnDeps <- FALSE
    }
    else
    {
        returnDeps <- f$returnDeps
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
    ans <- .Call('model', PACKAGE="RSiena",
                 z$Deriv, f$pData, seeds,
                 fromFiniteDiff, f$pModel, f$myeffects, z$theta,
                 randomseed2, returnDeps, z$FinDiff.method)
   #  browser()
   if (!fromFiniteDiff)
    {
        if (z$FinDiff.method)
            f$seeds <- ans[[3]]
    }
    if (z$Deriv)
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
        sims <- ans[[6]]
    else
        sims <- NULL
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
 #   cat('fra', fra, '\n')
 #    cat(f$randomseed2, f$storedseed, '\n')
   list(sc = sc, fra = fra, ntim0 = ntim, feasible = TRUE, OK = TRUE,
         sims=sims)
}
##@clearData siena07 Finalizer to clear Data object in C++
clearData <- function(pData)
{
    ans <- .Call('deleteData', PACKAGE="RSiena",
                 pData)
}
##@clearModel siena07 Finalizer to clear Model object in C++
clearModel <- function(pModel)
{
    ans <- .Call('deleteModel', PACKAGE="RSiena",
                 pModel)
}
##@createEdgeLists siena07 Reformat data for C++
createEdgeLists<- function(mat, matorig)
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
     },y = matorig)
    mat2 <- do.call(rbind, tmp)
    ## remove the diagonal
    mat2 <- mat2[mat2[, 1] != mat2[, 2], , drop=FALSE]
    ## mat3 structurals
    struct <- mat1[,3] %in% c(10, 11)
    mat1[struct, 3] <- mat1[struct,3] - 10
    mat3 <- mat1[struct, , drop=FALSE]
    mat3[, 3] <- 1
    mat1 <- mat1[!mat1[,3] == 0, , drop=FALSE] ##remove any zeros just created
    ##fix up storage mode to be integer
    storage.mode(mat1) <- 'integer'
    storage.mode(mat2) <- 'integer'
    storage.mode(mat3) <- 'integer'
    ## add attribute of size
    attr(mat1,'nActors') <- nrow(mat)
    attr(mat2,'nActors') <- nrow(mat)
    attr(mat3,'nActors') <- nrow(mat)

    list(mat1 = t(mat1), mat2 = t(mat2), mat3 = t(mat3))
}
##@createCovarEdgeLists siena07 Reformat data for C++
createCovarEdgeList<- function(mat)
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
    ##drop the diagonal : no, in case bipartite
   ## mat1 <- mat1[mat1[,1] != mat1[, 2],]
    ## add attribute of size
    attr(mat1,'nActors1') <- nrow(mat)
    attr(mat1,'nActors2') <- ncol(mat)
    t(mat1)
}
##@unpackOneMode siena07 Reformat data for C++
unpackOneMode <- function(depvar, observations, compositionChange)
{
    edgeLists <- vector('list', observations)
    networks <- vector('list', observations)
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
    sparse <- attr(depvar, 'sparse')
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
            for (j in seq(along=netmiss[[i]][,1]))
            {
                if (i == 1) # set missings to zero
                {
                    networks[[i]][netmiss[[i]][j, 1],
                                  netmiss[[i]][j, 2]] <- 0
                }
                else
                {
                    networks[[i]][netmiss[[i]][j, 1], netmiss[[i]][j, 2]] <-
                        networks[[i-1]][netmiss[[i]][j, 1], netmiss[[i]][j, 2]]
                }
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
            storage.mode(mat1) <- 'integer'
            storage.mode(mat2) <- 'integer'
            storage.mode(mat3) <- 'integer'
            ## add attribute of size
            attr(mat1,'nActors') <- nActors
            attr(mat2,'nActors') <- nActors
            attr(mat3,'nActors') <- nActors
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
                attr(depvar, 'distance')[i] <- sum(mydiff != 0,
                                                   na.rm = TRUE)
                if (all(mydiff@x >= 0, na.rm=TRUE))
                    attr(depvar, 'uponly')[i] <- TRUE
                if (all(mydiff@x <= 0, na.rm=TRUE))
                    attr(depvar, 'downonly')[i] <- TRUE
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
            if (i < observations)
            {
                ## recreate distances, as we have none in c++. (no longer true)
                mymat1 <- depvar[,,i, drop=FALSE]
                mymat2 <- depvar[,,i + 1,drop=FALSE]
                ##remove structural values
                mymat1[mymat1 %in% c(10,11)] <- NA
                mymat2[mymat2 %in% c(10,11)] <- NA
                ## and the diagonal
                diag(mymat1[,,1]) <- 0
                diag(mymat2[,,1]) <- 0
                mydiff <- mymat2 - mymat1
                attr(depvar, 'distance')[i] <- sum(mydiff != 0,
                                                         na.rm = TRUE)
                if (all(mydiff >= 0, na.rm=TRUE))
                    attr(depvar, 'uponly')[i] <- TRUE
                if (all(mydiff <= 0, na.rm=TRUE))
                    attr(depvar, 'downonly')[i] <- TRUE
            }

            diag(networks[[i]]) <- 0
            edgeLists[[i]] <- createEdgeLists(networks[[i]], depvar[, , i])
        }
    }
    ## add attribute of nodeset
    attr(edgeLists, 'nodeSet') <- attr(depvar, 'nodeSet')
    ## add attribute of name
    attr(edgeLists, 'name') <- attr(depvar, 'name')
    ## add attribute of distance
    attr(edgeLists, 'distance') <- attr(depvar, 'distance')
    ## attr uponly and downonly
    attr(edgeLists, 'uponly') <- attr(depvar, 'uponly')
    attr(edgeLists, 'downonly') <- attr(depvar, 'downonly')
    ## attr symmetric
    attr(edgeLists, 'symmetric') <- attr(depvar, 'symmetric')
    ## attr balmean
    attr(edgeLists, 'balmean') <- attr(depvar, 'balmean')
    return(edgeLists = edgeLists)
}
##@unpackBipartite siena07 Reformat data for C++
unpackBipartite <- function(depvar, observations, compositionChange)
{
    edgeLists <- vector('list', observations)
    networks <- vector('list', observations)
    actorSet <- attr(depvar, "nodeSet")
    compActorSets <- sapply(compositionChange, function(x)attr(x, "nodeSet"))
    sparse <- attr(depvar, 'sparse')
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
            for (j in seq(along=netmiss[[i]][,1]))
            {
                if (i == 1) # set missings to zero
                {
                    networks[[i]][netmiss[[i]][j, 1],
                                  netmiss[[i]][j, 2]] <- 0
                }
                else
                {
                    networks[[i]][netmiss[[i]][j, 1], netmiss[[i]][j, 2]] <-
                        networks[[i-1]][netmiss[[i]][j, 1], netmiss[[i]][j, 2]]
                }
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
            ##fix up storage mode to be integer
            storage.mode(mat1) <- 'integer'
            storage.mode(mat2) <- 'integer'
            storage.mode(mat3) <- 'integer'
            ## add attribute of size
            attr(mat1,'nActors') <- nActors
            attr(mat2,'nActors') <- nActors
            attr(mat3,'nActors') <- nActors
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
                attr(depvar, 'distance')[i] <- sum(mydiff != 0,
                                                   na.rm = TRUE)
                if (all(mydiff@x >= 0, na.rm=TRUE))
                    attr(depvar, 'uponly')[i] <- TRUE
                if (all(mydiff@x <= 0, na.rm=TRUE))
                    attr(depvar, 'downonly')[i] <- TRUE
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
            if (i < observations)
            {
                ## recreate distances, as we have none in c++. (no longer true)
                mymat1 <- depvar[,,i, drop=FALSE]
                mymat2 <- depvar[,,i + 1,drop=FALSE]
                ##remove structural values
                mymat1[mymat1 %in% c(10,11)] <- NA
                mymat2[mymat2 %in% c(10,11)] <- NA
                ## and the diagonal
                diag(mymat1[,,1]) <- 0
                diag(mymat2[,,1]) <- 0
                mydiff <- mymat2 - mymat1
                attr(depvar, 'distance')[i] <- sum(mydiff != 0,
                                                         na.rm = TRUE)
                if (all(mydiff >= 0, na.rm=TRUE))
                    attr(depvar, 'uponly')[i] <- TRUE
                if (all(mydiff <= 0, na.rm=TRUE))
                    attr(depvar, 'downonly')[i] <- TRUE
            }

            diag(networks[[i]]) <- 0
            edgeLists[[i]] <- createEdgeLists(networks[[i]], depvar[, , i])
        }
    }
    ## add attribute of nodeset
    attr(edgeLists, 'nodeSet') <- attr(depvar, 'nodeSet')
    ## add attribute of name
    attr(edgeLists, 'name') <- attr(depvar, 'name')
    ## add attribute of distance
    attr(edgeLists, 'distance') <- attr(depvar, 'distance')
    ## attr uponly and downonly
    attr(edgeLists, 'uponly') <- attr(depvar, 'uponly')
    attr(edgeLists, 'downonly') <- attr(depvar, 'downonly')
    ## attr symmetric
    attr(edgeLists, 'symmetric') <- attr(depvar, 'symmetric')
    ## attr balmean
    attr(edgeLists, 'balmean') <- attr(depvar, 'balmean')
    return(edgeLists = edgeLists)
}
##@unpackBehavior siena07 Reformat data for C++
unpackBehavior<- function(depvar, observations)
{
    beh <- depvar[, 1, ]
    behmiss <- is.na(beh)
    ## carry forward missings ### nb otherwise use the mode
    for (i in 1:observations)
    {
        if (i == 1)
            beh[is.na(beh[, i]), i] <- 0
        else ##carry missing forward!
            beh[is.na(beh[, i]), i] <-
                beh[is.na(beh[, i]), i - 1]
    }
    struct <- beh[beh %in% c(10,11)]
    beh[struct] <- beh[struct] - 10
    behstruct <- beh
    behstruct[!struct] <- 0
    ## add attribute of nodeset
    attr(beh, 'nodeSet') <- attr(depvar, 'nodeSet')
    ## add attribute of name
    attr(beh, 'name') <- attr(depvar, 'name')
    ## attr uponly and downonly
    attr(beh, 'uponly') <- attr(depvar, 'uponly')
    attr(beh, 'downonly') <- attr(depvar, 'downonly')
    ## attr symmetric
    attr(beh, 'symmetric') <- attr(depvar, 'symmetric')
    ## attr distance
    attr(beh, 'distance') <- attr(depvar, 'distance')
    ## attr simMean
    attr(beh, 'simMean') <- attr(depvar, 'simMean')
    storage.mode(beh) <- 'integer'
    list(beh=beh, behmiss=behmiss)
}
##@convertToStructuralZeros Miscellaneous To be implemented
convertToStructuralZeros <- function()
{
}

##@unpackCDyad siena07 Reformat data for C++
unpackCDyad<- function(dycCovar)
{
    sparse <- attr(dycCovar, 'sparse')
    if (sparse)
    {
        ## have a sparse matrix in triplet format
        ## with missings embedded - not allowed!
        ## with 0 based indices!
            varmat <- cbind(dycCovar@i+1, dycCovar@j+1, dycCovar@x)
            ##drop the diagonal, if present - not for bipartite
           ## varmat <- varmat[varmat[,1] != varmat[, 2],]
            mat1 <- varmat[!is.na(varmat[, 3]), ]
            mat1 <- mat1[!mat1[, 3] == 0, ]
            ## add attribute of dim
            attr(mat1,'nActors1') <- nrow(dycCovar)
            attr(mat1,'nActors2') <- ncol(dycCovar)
            edgeLists <-  t(mat1)
    }
    else
    {
        edgeLists <- createCovarEdgeList(dycCovar)
    }
    ## add attribute of nodesets
    attr(edgeLists, 'nodeSet') <- attr(dycCovar, 'nodeSet')
    ## add attribute of name
    attr(edgeLists, 'name') <- attr(dycCovar, 'name')
    ## add attribute of mean
    attr(edgeLists, 'mean') <- attr(dycCovar, 'mean')
    return(edgeLists = edgeLists)
}


##@unpackVDyad siena07 Reformat data for C++
unpackVDyad<- function(dyvCovar, observations)
{
    edgeLists <- vector('list', observations)
    varmats <- vector('list', observations)
    sparse <- attr(dyvCovar, 'sparse')
    if (sparse)
    {
        ## have a list of sparse matrices in triplet format
        ## with 0 based indices!
        for (i in 1:observations)
        {
            thisvar <- dyvCovar[[i]]
            varmat <- cbind(var@i+1, var@j+1, var@x)
            ##drop the diagonal, if present no - bipartite?
           ## varmat <- varmat[varmat[,1] != varmat[, 2],]
            mat1 <- varmat[!is.na(varmat[, 3]), ]
            mat1 <- mat1[!mat1[, 3] == 0, ]
            ## add attribute of size
            attr(mat1, 'nActors1') <- nrow(dyvCovar[[i]])
            attr(mat1, 'nActors2') <- ncol(dyvCovar[[i]])
            edgeLists[[i]] <- t(mat1)
        }
    }
    else
    {
        for (i in 1:(observations - 1))
        {
            thisvar <- dyvCovar[, , i]
            thisvar[is.na(thisvar) ]<- 0
            edgeLists[[i]] <- createCovarEdgeList(thisvar)
        }
    }
    ## add attribute of nodeset
    attr(edgeLists, 'nodeSet') <- attr(dyvCovar, 'nodeSet')
    ## add attribute of name
    attr(edgeLists, 'name') <- attr(dyvCovar, 'name')
     ## add attribute of mean
    attr(edgeLists, 'mean') <- attr(dyvCovar, 'mean')
    return(edgeLists = edgeLists)
}

##@unpackData siena07 Reformat data for C++
unpackData <- function(data)
{
    f <- NULL
    observations<- data$observations
    types <- sapply(data$depvars, function(x) attr(x, 'type'))
    f$nDepvars <- length(data$depvars)
    oneModes <- data$depvars[types == 'oneMode']
    Behaviors <- data$depvars[types == 'behavior']
    bipartites <- data$depvars[types == 'bipartite']
    f$nets <- lapply(oneModes, function(x, n, comp) unpackOneMode(x, n, comp),
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
    f$seed<- vector('list', observations - 1)
    f$depvars <- data$depvars
    f$nodeSets <- data$nodeSets
    f$oneModes <- oneModes
    f$Behaviors <- Behaviors
    f$oneModeUpOnly <- sapply(oneModes, function(x) attr(x, 'uponly'))
    f$oneModeDownOnly <- sapply(oneModes, function(x) attr(x, 'downonly'))
    f$behaviorUpOnly <- sapply(Behaviors, function(x) attr(x, 'uponly'))
    f$behaviorDownOnly <- sapply(Behaviors, function(x) attr(x,
                                                             'downonly'))
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
    nActors <- nrow(activeStart)
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
    interactions <- effects[effects$shortName == "unspInt" &
                            effects$effect1 > 0, ]
    if (nrow(interactions) > 0)
    {
        unspIntNames <- sapply(1:nrow(interactions), function(x, y, z)
           {
               y <- y[x, ] ## get the interaction effect
               twoway <- y$effect3 == 0
               ## now get the rows which are to interact
               inter1 <- z[z$effectNumber == y$effect1, ]
               if (nrow(inter1) != 1 )
               {
                   stop("invalid interaction specification effect number 1")
               }
               inter2 <- z[z$effectNumber == y$effect2, ]
               if (nrow(inter2) != 1 )
               {
                   stop("invalid interaction specification effect number 2")
               }
               if (!twoway)
               {
                   inter3 <- z[z$effectNumber == y$effect3, ]
                   if (nrow(inter3) != 1)
                   {
                       stop("invalid interaction specification effect number 3")
                   }
               }
               else
               {
                   inter3 <- z[is.na(z$effectNumber),] ## should be empty row
               }
               if (twoway)
               {
                   if (inter1$name != inter2$name)
                   {
                       stop("invalid interaction specification: ",
                            "must be same network")
                   }
                   if (inter1$type != inter2$type)
                   {
                       stop("invalid interaction specification: ",
                            "must be same type: evaluation or endowment")
                   }
               }
               else
               {
                   if (inter1$name != inter2$name ||
                       inter1$name != inter3$name)
                   {
                       stop("invalid interaction specification:",
                            "must all be same network")
                   }
                   if (inter1$type != inter2$type ||
                       inter1$type != inter3$type)
                   {
                       stop("invalid interaction specification:",
                            "must all be same type: evaluation or endowment")
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
                       stop("invalid interaction specification:",
                            "must be at least one ego or both dyadic effects")
                   }
               }
               else
               {
                   if (egoCount < 2 && dyadCount != 3)
                   {
                       stop("invalid interaction specification:",
                            "must be at least two ego or all dyadic effects")
                   }
               }
               ## construct a name
               ### make sure the egos are at the front of inters
               inters <- rbind(inters[egos, ], inters[-egos, ])
               tmpname <- paste(inters$effectName, collapse = " x ")
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
        effects[effects$shortName == "unspInt" &
                !is.na(effects$effect1), c("effectName", "functionName")] <-
                    unspIntNames
    }
    effects
}
