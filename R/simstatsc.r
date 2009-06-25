#/******************************************************************************
# * SIENA: Simulation Investigation for Empirical Network Analysis
# *
# * Web: http://www.stats.ox.ac.uk/~snidjers/siena
# *
# * File: siena01.r
# *
# * Description: This module contains the code for simulating the process,
# * communicating with C++.
# *****************************************************************************/
simstats0c <-function(z, x, INIT=FALSE, TERM=FALSE, initC=FALSE, data=NULL,
                        effects=NULL, fromFiniteDiff=FALSE,
                      profileData=FALSE, prevAns=NULL)
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
            effects <- effects[effects$include,]
            ## should split and rejoin before continuing
            effects1 <- split(effects, effects$name)
            if (inherits(data, "sienaGroup"))
                depvarnames <- names(data[[1]]$depvars)
            else
                depvarnames <- names(data$depvars)
            effects1order <- match(names(effects1), depvarnames)
            effects <- do.call(rbind, effects1[effects1order])
            row.names(effects) <- 1:nrow(effects)
            z$theta <- effects$initialValue
            z$fixed <- effects$fix
            z$test <- effects$test
            z$pp <- length(z$test)
            z$posj <- rep(FALSE,z$pp)
            z$posj[effects$basicRate] <- TRUE
            z$BasicRateFunction <- z$posj
            ## fix up effect names
            assort <- grep("assortativity", effects$effectName)
            if (length(assort) > 0 &  any(effects$parm[assort] != 2))
            {
                effects$functionName[assort] <-
                    ifelse(effects$parm[assort] == 2,
                           effects$functionName[assort],
                           sub("^(1/2)", "", effects$functionName[assort],
                               fixed=TRUE))
                effects$effectName[assort] <-
                    ifelse(effects$parm[assort] == 2,
                           effects$effectName[assort],
                           sub("^(1/2)", "", effects$effectName[assort],
                               fixed=TRUE))
            }
            if (inherits(data, 'sienaGroup'))
            {
                nGroup <- length(data)
            }
            else
            {
                nGroup <- 1
                data <- sienaGroupCreate(list(data), singleOK=TRUE)
            }
            if (x$cconditional)
            {
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
                z$condvar <- (1:nrow(effects))[effects$name==z$condname][1:
                                               observations]
                z$theta<- z$theta[-z$condvar]
                z$fixed<- z$fixed[-z$condvar]
                z$test<- z$test[-z$condvar]
                z$pp<- z$pp-length(z$condvar)
                z$scale<- z$scale[-z$condvar]
                z$BasicRateFunction <- z$posj[-z$condvar]
                z$posj <- z$posj[-z$condvar]
                z$theta[z$posj] <-
                    z$theta[z$posj] / effects$initialValue[z$condvar]
                z$ntim<- matrix(NA, nrow=x$n3, ncol=observations)
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
            if (x$cconditional)
            {
                attr(f, "change") <-
                    sapply(f, function(xx)attr(xx$depvars[[z$condname]],
                                               'distance'))
                attr(f,"condEffects") <- effects[z$condvar,]
                effects <- effects[-z$condvar, ]
            }
            ## see if we can use the original dfra
            if (!is.null(prevAns) && inherits(prevAns, "sienaFit"))
            {
                if (all(names(prevAns$dfra) == effects$shortNames)
                    && !is.null(prevAns$sf))
                {
                    z$haveDfra <- TRUE
                    z$dfra <- prevAns$dfra
                    z$sf <- prevAns$sf
                }
            }
            z$effects <- effects
        }
        pData <- .Call('setupData', PACKAGE="RSiena",
                       lapply(f, function(x)(as.integer(x$observations))),
                       lapply(f, function(x)(x$nodeSets)))
        ans<- .Call('OneMode', PACKAGE="RSiena",
                    pData, lapply(f, function(x)x$nets))
        ans <- .Call('Behavior', PACKAGE="RSiena",
                     pData, lapply(f, function(x)x$behavs))
        ans<-.Call('ConstantCovariates', PACKAGE="RSiena",
                   pData, lapply(f, function(x)x$cCovars))
        ans<-.Call('ChangingCovariates', PACKAGE="RSiena",
                   pData,lapply(f, function(x)x$vCovars))
        ans<-.Call('DyadicCovariates', PACKAGE="RSiena",
                   pData,lapply(f, function(x)x$dycCovars))
        ans<-.Call('ChangingDyadicCovariates', PACKAGE="RSiena",
                   pData, lapply(f, function(x)x$dyvCovars))
        ans<-.Call('ExogEvent', PACKAGE="RSiena",
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
            myeffects <- f$myeffects
        }
        ans<- .Call('effects', PACKAGE="RSiena",
                    pData, myeffects)
        pModel <- ans[[1]][[1]]
        for (i in 1:length(ans[[2]])) ## ans[[2]] is a list of lists of
            ## pointers to effects. Each list corresponds to one
            ## dependent variable
        {
            effectPtr <- ans[[2]][[i]]
            myeffects[[i]]$effectPtr <- effectPtr
        }
        if (!initC)
        {
            ans <- .Call('getTargets', PACKAGE="RSiena",
                         pData, pModel, myeffects)
            z$targets <- rowSums(ans)
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
        if (x$cconditional)
        {
            CONDVAR <- z$condname
            CONDTARGET <- attr(f, "change")
        }
        else
        {
            CONDVAR <- NULL
            CONDTARGET <- NULL
        }
        ans <- .Call("setupModelOptions", PACKAGE="RSiena",
                     pData, pModel, MAXDEGREE, CONDVAR, CONDTARGET,
                     profileData)
        f$myeffects <- myeffects
        z$effects <- effects
        if (!initC)
        {
            DataReport(z, x, f)
        }
        f$observations <- attr(f, "observations") + 1
        f$randomseed2 <- z$randomseed2
        FRANstore(f) ## store f in FRANstore
        return(z)
    }
    if (TERM)
    {
        if (x$cconditional)
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
        z$f <- f
        PrintReport(z, x)
        if (sum(z$test))
        {
            z$fra <- colMeans(z$sf, na.rm=TRUE)
            z <- ScoreTest(z, x)
            TestOutput(z, x)
        }
        dimnames(z$dfra)[[1]] <- z$effects$shortName
        return(z)
    }
    ## iteration entry point
    f <- FRANstore()

    if (is.null(f$randomseed2))
    {
        randomseed2 <- NULL
    }
    else
    {
        randomseed2 <- as.integer(f$randomseed2)
    }
    ans <- .Call('model', PACKAGE="RSiena",
                 z$Deriv, f$pData, f$seeds,
                 fromFiniteDiff, f$pModel, f$myeffects, z$theta,
                 randomseed2)
    ## browser()
    if (!fromFiniteDiff)
    {
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
    f$randomseed2 <- ans[[5]]
    FRANstore(f)
    list(sc = sc, fra = fra, ntim0 = ntim, feasible = TRUE, OK = TRUE)
}
clearData <- function(pData)
{
    ans <- .Call('deleteData', PACKAGE="RSiena",
                 pData)
}
clearModel <- function(pModel)
{
    ans <- .Call('deleteModel', PACKAGE="RSiena",
                 pModel)
}
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
    mat2 <- mat2[mat2[, 1] != mat2[, 2], ]
    ## mat3 structurals
    struct <- mat1[,3] %in% c(10, 11)
    mat1[struct,3] <- mat1[struct,3] - 10
    mat3 <- mat1[struct, , drop=FALSE]
    mat1 <- mat1[!mat1[,3] == 0, ] ##remove any zeros just created
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
        require(Matrix)
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
                netmiss[[i]][netmiss[[i]][, 1] != netmiss[[i]][, 2], ]
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
convertToStructuralZeros <- function()
{
}
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

unpackData <- function(data)
{
    f <- NULL
    observations<- data$observations
    types <- sapply(data$depvars, function(x) attr(x, 'type'))
    f$nDepvars <- length(data$depvars)
    oneModes <- data$depvars[types == 'oneMode']
    Behaviors <- data$depvars[types == 'behavior']
    f$nets <- lapply(oneModes, function(x, n, comp) unpackOneMode(x, n, comp),
                     n = observations, comp=data$compositionChange)
    names(f$nets) <- names(oneModes)
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
    if (is.null(f$nets))
        f$nActors <- nrow(f$behavs[[1]])
    else
        f$nActors <- attr(f$nets[[1]][[1]]$mat1,'nActors')
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
