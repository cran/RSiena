#/******************************************************************************
# * SIENA: Simulation Investigation for Empirical Network Analysis
# *
# * Web: http://www.stats.ox.ac.uk/~snidjers/siena
# *
# * File: sienaDataCreate.r
# *
# * Description: This module contains the code to create
# * Siena data object and group data objects.
# *****************************************************************************/
sienaDataCreate<- function(..., nodeSets=NULL)
{
    validNodeSet <- function(nodeSetName, n)
    {
        if (!nodeSetName == 'Actors')
        {
            sub <- match(nodeSetName, nodeSetNames)
            if (is.na(sub))
                stop('node set not found')
            n == length(nodeSets[[sub]])
        }
        else
            TRUE
    }
    narg <- nargs()
    dots <- as.list(substitute(list(...)))[-1] ##first entry is the word 'list'
    if (length(dots) == 0)
        stop('need some networks')
    nm <- names(dots)
    if (is.null(nm))
        fixup <- seq(along=dots)
    else
        fixup <- nm == ''
    dep <- sapply(dots[fixup], function(x) deparse(x)[1])
    if (is.null(nm))
        nm <-dep
    else if (length(dep) > 0)
        nm[fixup] <- dep
    dots <- list(...)
    names(dots)<- nm
    if (any(duplicated(nm)))
        stop('names must be unique')
    ## process the inputs: check dimensions,
    ## sort out missings and structural zeros and symmetric etc
    ## check sizes match the corresponding nodeSets
    observations <- 0
    depvars <- vector('list',narg)
    cCovars <- vector('list',narg)
    vCovars <- vector('list',narg)
    dycCovars <- vector('list',narg)
    dyvCovars <- vector('list',narg)
    compositionChange <- vector('list',narg)
    v1 <- 0; v2 <- 0; v3 <- 0; v4 <- 0; v5 <- 0; v6 <- 0
    for (i in seq(along = dots))
        switch(class(dots[[i]]),
               sienaNet = {
                   if (attr(dots[[i]],'sparse'))
                   {
                       require(Matrix)
                       netdims <- c(dim(dots[[i]][[1]]), length(dots[[i]]))
                   }
                   else
                       netdims <- dim(dots[[i]])
                   if (observations == 0)
                       observations <- netdims[3]
                   else if (observations != netdims[3])
                       stop('differing number of observations')
                   v1 <- v1+1
                   depvars[[v1]] <- dots[[i]]
                   names(depvars)[v1] <- nm[i]
               },
               coCovar = {
                   v2<- v2+1
                   cCovars[[v2]] <- dots[[i]]
                   names(cCovars)[v2] <- nm[i]
              },
               varCovar = {
                   v3<- v3+1
                   vCovars[[v3]]<- dots[[i]]
                   names(vCovars)[v3]<- nm[i]
               },
               coDyadCovar={
                   v4<- v4+1
                   dycCovars[[v4]]<- dots[[i]]
                   names(dycCovars)[v4]<- nm[i]
               },
               varDyadCovar={
                   v5<- v5+1
                   dyvCovars[[v5]]<- dots[[i]]
                   names(dyvCovars)[v5]<- nm[i]
               },
               compositionChange={
                   v6 <- v6 + 1
                   compositionChange[[v6]] <- dots[[i]]
                   names(compositionChange)[v6] <- nm[i]
               },
               stop(paste('invalid object in sienaDataCreate',
                          class(dots[[i]])),call.=FALSE)
               )
    if (v1 == 0)
        stop('need a network')
    depvars <- depvars[1:v1]
    if (is.null(nodeSets))
        nodeSets <- list(sienaNodeSet(attr(depvars[[1]],'netdims')[1]))
    nodeSetNames <- sapply(nodeSets,function(x) attr(x,'nodeSetName'))
    if (v2 == 0)
        cCovars <- list()
    else
        cCovars <- cCovars[1:v2]
    if (v3 == 0)
        vCovars <- list()
    else
        vCovars <- vCovars[1:v3]
    if (v4 == 0)
        dycCovars <- list()
    else
        dycCovars <- dycCovars[1:v4]
    if (v5 == 0)
        dyvCovars <- list()
    else
        dyvCovars <- dyvCovars[1:v5]
    if (v6 == 0)
        compositionChange <- list()
    else
        compositionChange <- compositionChange[1:v6]
    ##now can check dimensions and find ranges
    for (i in seq(along = cCovars))
    {
        if (!validNodeSet(attr(cCovars[[i]], 'nodeSet'), length(cCovars[[i]])))
            stop('constant covariate incorrect size: ', names(cCovars)[i])
        varmean <- mean(cCovars[[i]] ,na.rm=TRUE)
        range2 <- range(cCovars[[i]], na.rm=TRUE)
        attr(cCovars[[i]], 'moreThan2') <- length(table(cCovars[[i]])) > 2
        vartotal <- sum(cCovars[[i]], na.rm=TRUE)
        nonMissingCount <- sum(!is.na(cCovars[[i]]))
        cCovars[[i]] <- cCovars[[i]] - varmean
        attr(cCovars[[i]], 'mean') <- varmean
        rr <- rangeAndSimilarity(cCovars[[i]], range2)
        if (rr$range[2] == rr$range[1] && !any(is.na(cCovars[[i]])))
                attr(cCovars[[i]], 'poszvar') <- FALSE
        else
            attr(cCovars[[i]], 'poszvar') <- TRUE
        attr(cCovars[[i]], 'range') <- rr$range[2] - rr$range[1]
        storage.mode(attr(cCovars[[i]], 'range')) <- 'double'
        attr(cCovars[[i]], 'range2') <- range2
       # attr(cCovars[[i]], 'simTotal') <- rr$simTotal
        attr(cCovars[[i]], 'simMean') <- rr$simMean
       # attr(cCovars[[i]], 'simCnt') <- rr$simCnt
        attr(cCovars[[i]], "name") <- names(cCovars)[i]
        attr(cCovars[[i]], "vartotal") <- vartotal
        attr(cCovars[[i]], "nonMissingCount") <- nonMissingCount
    }
    for (i in seq(along=vCovars)) ## note that behaviour variables are not here!
    {
        if (!validNodeSet(attr(vCovars[[i]], 'nodeSet'), nrow(vCovars[[i]])))
            stop('changing covariate incorrect size: ', names(vCovars)[i])
        if (ncol(vCovars[[i]]) < (observations - 1))
            stop('changing covariate not enough columns')
        if (ncol(vCovars[[i]]) != (observations - 1))
        {
            tmpatt <- attributes(vCovars[[i]])
            vCovars[[i]] <- vCovars[[i]][, 1:(observations - 1), drop=FALSE]
            attnames <- names(tmpatt)
            for (att in seq(along=attnames))
            {
                if (!attnames[att] %in% c('dim', 'dimnames'))
                {
                    attr(vCovars[[i]], attnames[att]) <- tmpatt[[att]]
                }
            }
        }
        tmpmat <- vCovars[[i]]
        varmean <- mean(tmpmat, na.rm=TRUE)
        vartotal <- sum(vCovars[[i]], na.rm=TRUE)
        nonMissingCount <- sum(!is.na(vCovars[[i]]))
        attr(vCovars[[i]], "rangep") <- apply(tmpmat, 2, range,
                                              na.rm=TRUE)
        attr(vCovars[[i]], "meanp") <- colMeans(tmpmat, na.rm=TRUE)
        cr <- range(tmpmat, na.rm=TRUE)
        attr(vCovars[[i]], 'range') <- cr[2] - cr[1]
        storage.mode(attr(vCovars[[i]], 'range')) <- 'double'
        attr(vCovars[[i]], 'mean') <- varmean
        vCovars[[i]] <- vCovars[[i]] - varmean
        rr <- rangeAndSimilarity(tmpmat, cr)
        if (rr$range[2] == rr$range[1] && !any(is.na(tmpmat)))
            attr(vCovars[[i]], 'poszvar') <- FALSE
        else
            attr(vCovars[[i]], 'poszvar') <- TRUE
        attr(vCovars[[i]], 'simMean') <- rr$simMean
        attr(vCovars[[i]], 'moreThan2') <- length(unique(vCovars[[i]])) > 2
        attr(vCovars[[i]], 'name') <- names(vCovars)[i]
        attr(vCovars[[i]], "vartotal") <- vartotal
        attr(vCovars[[i]], "nonMissingCount") <- nonMissingCount
    }
    for (i in seq(along=dycCovars))
    {
        nattr <- attr(dycCovars[[i]], 'nodeSet')
        if (!validNodeSet(nattr[1], nrow(dycCovars[[i]])))
            stop('dyadic covariate incorrect nbr rows', names(dycCovars)[i])
        if (!validNodeSet(nattr[2], ncol(dycCovars[[i]])))
             stop('dyadic covariate incorrect nbr columns', names(dycCovars)[i])
        diag(dycCovars[[i]]) <- NA
        varmean <- mean(dycCovars[[i]], na.rm=TRUE)
        attr(dycCovars[[i]],'mean') <- varmean
        rr<-  range(dycCovars[[i]], na.rm=TRUE)
        attr(dycCovars[[i]],'range') <- rr[2] - rr[1]
        storage.mode(attr(dycCovars[[i]], 'range')) <- 'double'
        attr(dycCovars[[i]],'range2') <- rr
        attr(dycCovars[[i]], 'name') <- names(dycCovars)[i]
        diag(dycCovars[[i]]) <- 0
    }
    for (i in seq(along=dyvCovars))
    {
        nattr <- attr(dyvCovars[[i]],'nodeSet')
        if (!validNodeSet(nattr[1], dim(dyvCovars[[i]])[1]))
            stop('dyadic changing covariate incorrect nbr rows',
                 names(dyvCovars)[i])
        if (!validNodeSet(nattr[2], dim(dyvCovars[[i]])[2]))
            stop('dyadic changing covariate incorrect nbr columns',
                 names(dyvCovars)[i])
        if (dim(dyvCovars[[i]])[3] < (observations - 1))
            stop('Dyadic changing covariate not enough observations')
         if (dim(dyvCovars[[i]])[3] != (observations - 1))
        {
            tmpatt <- attributes(dyvCovars[[i]])
            dyvCovars[[i]] <- dyvCovars[[i]][, 1:(observations - 1)]
            attnames <- names(tmpatt)
            for (att in seq(along=attnames))
            {
                if (attnames[att] != "dim")
                {
                    attr(dyvCovars[[i]], attnames[att]) <- tmpatt[[att]]
                }
            }
          }
        for (obs in 1:dim(dyvCovars[[i]])[3])
            diag(dyvCovars[[i]][,,obs]) <- NA
        varmean <- mean(dyvCovars[[i]], na.rm=TRUE)
        attr(dyvCovars[[i]],'mean') <- mean(dyvCovars[[i]], na.rm=TRUE)
        rr <-  range(dyvCovars[[i]], na.rm=TRUE)
        attr(dyvCovars[[i]],'range') <- rr[2] - rr[1]
        storage.mode(attr(dyvCovars[[i]], 'range')) <- 'double'
        attr(dyvCovars[[i]], 'name') <- names(dyvCovars)[i]
        for (obs in 1:dim(dyvCovars[[i]])[3])
            diag(dyvCovars[[i]][,,obs]) <- 0
    }
    compnodesets <- sapply(compositionChange, function(x) attr(x, 'nodeSet'))
    if (any(duplicated(compnodesets)))
        stop('Only one composition change allowed for each nodeSet')
    for (i in seq(along = compositionChange))
    {
        thisNodeSet <- attr(compositionChange[[i]], 'nodeSet')
        nodeSetSize <- length(compositionChange[[i]])
        if (!validNodeSet(thisNodeSet, nodeSetSize))
            stop('composition change incorrect size: ',
                 names(compositionChange)[i])
        if (any(sapply(compositionChange[[i]], function(x)
                       any(x < 1.0 | x > observations))))
            stop("invalid times of composition change")
        if (!all(sapply(compositionChange[[i]], length) %% 2 == 0))
            stop(" Each composition change entry must have an ",
                 "even number of digits")
        ## generate events and active flags
        activeStart <- matrix(FALSE, nrow=nodeSetSize, ncol=observations)
        action <- matrix(0, nrow=nodeSetSize, ncol=observations)
        events <- vector("list", nodeSetSize * 2 * observations)
        evSubs <- 1
        for (j in 1:nodeSetSize)
        {
            xsubs <- 1
            x <- compositionChange[[i]][[j]]
            repeat
            {
                ##process one interval
                start <- x[xsubs]
                end <- x[xsubs+1]
                startIndex <- ceiling(x[xsubs])
                endIndex <- trunc(x[xsubs + 1])
              #  if (startIndex < observations && startIndex <= activeEndIndex)
              #  {
                    activeStart[j, startIndex:endIndex] <- TRUE
              #  }
                if (x[xsubs] > 1.0)
                {
                    period <- trunc(x[xsubs])
                    evTime <- x[xsubs] - period
                    events[[evSubs]] <- data.frame(event="join",
                                                   period=period,
                                                   actor = j, time=evTime)
                    evSubs <- evSubs + 1
                }
                if (x[xsubs+1] < observations)
                {
                    period <- trunc(x[xsubs+1])
                    evTime <- x[xsubs+1] - period
                    events[[evSubs]] <- data.frame(event="leave",
                                                   period=period,
                                                   actor = j, time=evTime)
                    evSubs <- evSubs + 1
                }
                xsubs <- xsubs + 2
                if (xsubs > length(x))
                {
                    break
                }
            }
         #   cat(j, 'active',activeStart[j,],  x, '\n')
            if (any(!activeStart[j,]))
            {
                notActive <- which(!activeStart[j,])
                for (jj in notActive)
                {
                    precActive <- jj > 1 && sum(activeStart[j, 1:(jj-1)]) > 0
                    len <- length(activeStart[i,])
                    succActive <- (jj < len) &&
                    (sum(activeStart[j, (jj + 1):len]) > 0)

                    if (!precActive)
                    {
                        action[j, jj] <- 1
                    }
                    else if (!succActive)
                    {
                        action[j, jj] <- 2
                    }
                    else
                    {
                        action[j, jj] <- 3
                    }
                }
            }
        }
        if (evSubs > 1)
        {
            events <- do.call(rbind, events[1:(evSubs-1)])
        }
        else
        {
            events <- data.frame(event="join", period=1, actor = 1, time=0)[-1,]
        }
        events$event <- factor(events$event, levels=c('join','leave'))
        storage.mode(events$period) <- "integer"
        storage.mode(events$actor) <- "integer"
        attr(compositionChange[[i]], "events") <- events
        attr(compositionChange[[i]], "activeStart") <- activeStart
        attr(compositionChange[[i]], "action") <- action
    }
    for (i in 1:v1)
    {
        nattr <- attr(depvars[[i]], 'nodeSet')
        netdims <- attr(depvars[[i]], 'netdims')
        type <- attr(depvars[[i]], 'type')
        sparse <- attr(depvars[[i]], 'sparse')
        myarray <- depvars[[i]]
        if (!validNodeSet(nattr[1],netdims[1]))
            stop('1st net dimension wrong')
        if (type =='bipartite')
            if (!validNodeSet(nattr[2], netdims[2]))
                stop('2nd net dimension wrong')
        attr(depvars[[i]], 'uponly') <- rep(FALSE, observations - 1)
        attr(depvars[[i]], 'downonly') <- rep(FALSE, observations - 1)
        attr(depvars[[i]], 'distance') <- rep(FALSE, observations - 1)
        attr(depvars[[i]], 'vals') <- vector("list", observations)
        attr(depvars[[i]], 'nval') <- rep(NA, observations)
        attr(depvars[[i]], 'noMissing') <- rep(NA, observations)
        if (type == 'behavior')
        {
            attr(depvars[[i]], 'noMissing') <- FALSE
            for (j in 1:(observations - 1))
            {
                myvector1 <- myarray[, , j]
                myvector2 <- myarray[, , j+1]
                mydiff <- abs(myvector1 - myvector2)
                attr(depvars[[i]], "distance")[j] <- sum(mydiff, na.rm=TRUE)
                attr(depvars[[i]], "vals")[[j]] <- table(myvector1,
                                                         useNA="always")
                attr(depvars[[i]], "vals")[[j+1]] <- table(myvector2,
                                                           useNA="always")

                attr(depvars[[i]], "nval")[j] <-  sum(!is.na(myvector1))
                attr(depvars[[i]], "nval")[j + 1] <-  sum(!is.na(myvector2))
                attr(depvars[[i]], 'noMissing')[j] <- sum(is.na(myvector1))
                attr(depvars[[i]], 'noMissing')[j+1] <- sum(is.na(myvector2))
            }
            rr <- range(depvars[[i]], na.rm=TRUE)
            if (rr[2] == rr[1] && !any(is.na(depvars[[i]])))
                attr(depvars[[i]], 'poszvar') <- FALSE
            else
                attr(depvars[[i]], 'poszvar') <- TRUE
            crange <- rr[2] - rr[1]
            attr(depvars[[i]],'range') <- crange
            attr(depvars[[i]], "range2") <- rr
            attr(depvars[[i]],'moreThan2') <- length(unique(depvars[[i]])) > 2
            modes <- apply(depvars[[i]][, 1, ], 2, function(z)
                       {
                           taba <- table(round(z))
                           as.numeric(names(which.max(taba)))
                       }
                           )
            attr(depvars[[i]],'modes') <- modes
            attr(depvars[[i]], 'missing') <- any(is.na(depvars[[i]]))
            tmpmat <- depvars[[i]][, 1, ]
            rr <- rangeAndSimilarity(tmpmat[, - ncol(tmpmat)], rr)
        #    attr(depvars[[i]], 'simTotal') <- rr$simTotal
         #   attr(depvars[[i]], 'simCnt') <- rr$simCnt
            attr(depvars[[i]], 'simMean') <- rr$simMean
            attr(depvars[[i]], 'structural') <- FALSE
        }
        else
        {
            for (j in 1:(observations - 1))
            {
                if (sparse)
                {
                    mymat1 <- myarray[[j]]
                    mymat2 <- myarray[[j + 1]]
                    ##remove structural values
                    x1 <- mymat1@x
                    x2 <- mymat2@x
                    x1[x1 %in% c(10, 11)] <- NA
                    x2[x2 %in% c(10, 11)] <- NA
                    mymat1@x <- x1
                    mymat2@x <- x2
                    mydiff <- mymat2 - mymat1
                    attr(depvars[[i]], 'distance')[j] <- sum(mydiff != 0,
                                                             na.rm = TRUE)
                    if (all(mydiff@x >= 0, na.rm=TRUE))
                        attr(depvars[[i]], 'uponly')[j] <- TRUE
                    if (all(mydiff@x <= 0, na.rm=TRUE))
                        attr(depvars[[i]], 'downonly')[j] <- TRUE
                }
                else
                {
                    mymat1 <- myarray[,,j, drop=FALSE]
                    mymat2 <- myarray[,,j + 1,drop=FALSE]
                    ##remove structural values
                    mymat1[mymat1 %in% c(10,11)] <- NA
                    mymat2[mymat2 %in% c(10,11)] <- NA
                    mydiff <- mymat2 - mymat1
                    attr(depvars[[i]], 'distance')[j] <- sum(mydiff != 0,
                                                             na.rm = TRUE)
                    if (all(mydiff >= 0, na.rm=TRUE))
                        attr(depvars[[i]], 'uponly')[j] <- TRUE
                    if (all(mydiff <= 0, na.rm=TRUE))
                        attr(depvars[[i]], 'downonly')[j] <- TRUE
                }
            }
            if (type == 'oneMode')
            {
                attr(depvars[[i]], 'balmean') <- calcBalmean(depvars[[i]])
                attr(depvars[[i]], 'symmetric') <- TRUE
                attr(depvars[[i]], 'missing') <- FALSE
                attr(depvars[[i]], 'structural') <- FALSE
                for (j in 1:observations)
                {
                    if (sparse)
                    {
                        mymat <- myarray[[j]]
                    }
                    else
                    {
                        mymat <- myarray[, , j]
                    }
                    if (!isSymmetric(mymat))
                    {
                        attr(depvars[[i]], 'symmetric') <- FALSE
                    }
                    if (any(is.na(mymat)))
                    {
                        attr(depvars[[i]], 'missing') <- TRUE
                    }
                    if (any(!is.na(mymat) & (mymat == 10 | mymat == 11)))
                    {
                        attr(depvars[[i]], "structural") <- TRUE
                    }
                    if (sparse)
                    {
                        nonZeros <- table(mymat@x, useNA="always")
                        Zeros <- nrow(mymat) * nrow(mymat - 1) - sum(nonZeros)
                        attr(depvars[[i]], "vals")[[j]] <-
                            c(table(rep(0, Zeros)), nonZeros)
                    }
                    else
                    {
                        attr(depvars[[i]], "vals")[[j]] <- table(mymat,
                                                                 useNA="always")
                    }
                    attr(depvars[[i]], "nval")[j] <-
                        sum(!is.na(mymat[row(mymat) != col(mymat)]))
                }
                ### need to exclude the structurals here
                if (sparse)
                {
                    vals <- lapply(depvars[[i]], function(x)
                                   c(x@x[!(is.na(x@x) |
                                           x@x %in% c(10, 11))] , 0))
                    attr(depvars[[i]], "range") <-
                        do.call(range, vals)
                }
                else
                {
                    tmp <- depvars[[i]]
                    attr(depvars[[i]], "range") <-
                        range(tmp[!(is.na(tmp) | tmp %in% c(10, 11))])
                }
            }
            else #type=='bipartite' not sure what we need here
            {
            }
        }
        attr(depvars[[i]], 'name') <- names(depvars)[i]
    }
    ## create the object
    z<- NULL
    z$nodeSets <- nodeSets
    z$observations <-observations
    z$depvars <- depvars
    z$cCovars <- cCovars
    z$vCovars <- vCovars
    z$dycCovars <- dycCovars
    z$dyvCovars <- dyvCovars
    z$compositionChange <- compositionChange
    class(z) <- 'siena'
    z
}

rangeAndSimilarity<- function(vals, rvals=NULL)
{
    vals <- as.matrix(vals)
    dims <- dim(vals)
    observations <- dims[2]
    n <- dims[1]
    if (is.null(rvals))
    {
        rvals <- range(vals,na.rm=TRUE)
    }
    rvals1 <- rvals[2] - rvals[1]
    ## use of outer creates a potentially huge matrix unnecessarily
    ## tmp3 <- apply(vals,2,function(x)
    ##          {
    ##              y <- 1 - abs(outer(x, x, "-")) / rvals1
    ##              diag(y) <- NA
    ##              y
    ##          })
    ## tmp3 <- array(tmp3, dim=c(n, n, observations))
    tmp <- apply(vals, 2, function(v){
    sapply(1: length(v), function(x, y, r){
        z <- y
        z[x] <- NA
        ##browser()
        tmp1 <- 1 - abs(y[x] - z) / r
        list(sum(tmp1, na.rm=TRUE), sum(!is.na(tmp1)))
         },
           y=v, r=rvals1)})

    tmp <- unlist(tmp)
    raw <- tmp[seq(1, length(tmp), by=2)]
    cnts <- tmp[seq(2, length(tmp), by=2)]
    simTotal <- sum(raw)
    simCnt <- sum(cnts)
    simMean <- simTotal/simCnt
    list(simTotal=simTotal, simMean=simMean, range=rvals, simCnt=simCnt)
}

groupRangeAndSimilarityAndMean <- function(group)
{
    atts <- attributes(group)
    behavs <- atts$types == "behavior"
    netnames <- atts$netnames
    bRange <- rep(NA, length(netnames))
    behRange <- matrix(NA, ncol=length(netnames), nrow=2)
    bSim <- rep(NA, length(netnames))
    bPoszvar <- rep(NA, length(netnames))
    bMoreThan2 <- rep(NA, length(netnames))
    bAnyMissing <- rep(FALSE, length(netnames))
    for (net in which(atts$types == "behavior"))
    {
        simTotal <- 0
        simCnt <- 0
        anyMissing <- FALSE
        bPoszvar[net] <- TRUE
        thisrange <- matrix(NA, ncol=length(group), nrow=2)
        for (i in 1:length(group))
        {
            j <- match(netnames[net], names(group[[i]]$depvars))
            if (is.na(j))
                stop("network names not consistent")
            depvar <- group[[i]]$depvars[[j]][, 1, ]
            ## this should be a matrix
            thisrange[, i] <- range(depvar, na.rm=TRUE)
            if (any(is.na(depvar)))
            {
                anyMissing <- TRUE
            }
        }
        behRange[, net] <- range(thisrange, na.rm=TRUE)
        bRange[net] <- behRange[, net][2] - behRange[, net][1]
        if (behRange[, net][2] == behRange[, net][1] && !anyMissing)
        {
            bPoszvar[net] <- FALSE
        }
        values <- NULL
        for (i in 1:length(group))
        {
            j <- match(netnames[net], names(group[[i]]$depvars))
            depvar <- group[[i]]$depvars[[j]][, 1, ]
            ## this should be a matrix

            tmp <- rangeAndSimilarity(depvar[, -ncol(depvar)],
                                      behRange[, net])
            simTotal <- simTotal + tmp$simTotal
            simCnt <- simCnt + tmp$simCnt
            values <- c(values, unique(depvar))
        }
        simMean <- simTotal/simCnt
        bSim[net] <- simMean
        bMoreThan2 <- length(unique(values)) > 2
        if (anyMissing)
            bAnyMissing[net] <- TRUE
    }
    ## constant ones will not exist unless there is only one data object
    cCovarRange <- rep(NA, length(atts$cCovars))
    cCovarSim <- rep(NA, length(atts$cCovars))
    cCovarPoszvar <- rep(TRUE, length(atts$cCovars))
    cCovarMoreThan2 <- rep(FALSE, length(atts$cCovars))
    cCovarMean <- rep(NA, length(atts$cCovars))
    cCovarRange2 <- matrix(NA, nrow=length(atts$cCovars), ncol=2)
    for (covar in seq(along=atts$cCovars))
    {
        if (length(group) > 1)
            stop("group create constant covariate error")
        simTotal <- 0
        simCnt <- 0
        anyMissing <- FALSE
        ## first find the range
        thisrange <- matrix(NA, ncol=length(group),nrow=2)
        for (i in 1:length(group))
        {

            j <- match(atts$cCovars[covar], names(group[[i]]$cCovars))
            if (is.na(j))
            {
                stop("inconsistent covariate names")
            }
            thisrange[, i] <- range(group[[i]]$cCovars[[j]],
                                   na.rm=TRUE)
            if (any(is.na(group[[i]]$cCovars[[j]])))
            {
                anyMissing <- TRUE
            }
        }
        rr <- range(thisrange, na.rm=TRUE)
        cCovarRange[covar] <- rr[2] - rr[1]
        if (rr[2] == rr[1] && !anyMissing)
                cCovarPoszvar[covar] <- FALSE
        ##then calculate similarity
        for (i in 1:length(group))
        {

            j <- match(atts$cCovars[covar], names(group[[i]]$cCovars))
            tmp <- rangeAndSimilarity(group[[i]]$cCovars[[covar]],
                                      rr)
            simTotal <- simTotal + tmp$simTotal
            simCnt <- simCnt + tmp$simCnt
        }
        simMean <- simTotal/simCnt
        cCovarSim[covar] <- simMean
        cCovarMoreThan2[covar] <- attr(group[[i]]$cCovars[[covar]], "moreThan2")
        cCovarMean[covar] <- attr(group[[i]]$cCovars[[covar]], "mean")
        cCovarRange2[covar, ] <- attr(group[[i]]$cCovars[[covar]], "range2")
   }

    vCovarRange <- rep(NA, length(atts$vCovar))
    vCovarSim <- rep(NA, length(atts$vCovar))
    vCovarPoszvar <- rep(TRUE, length(atts$vCovar))
    vCovarMoreThan2 <- rep(FALSE, length(atts$vCovar))
    vCovarMean <- rep(NA, length(atts$vCovar))
    for (covar in seq(along=atts$vCovars))
    {
        vartotal <- 0
        nonMissingCount <- 0
        ## need to re-centre these values
        for (i in 1:length(group))
        {
            j <- match(atts$vCovars[covar], names(group[[i]]$vCovars))
            if (is.na(j))
            {
                stop("inconsistent covariate names")
            }
            vartotal <- vartotal + attr(group[[i]]$vCovars[[j]], "vartotal")
            nonMissingCount <- nonMissingCount +
                attr(group[[i]]$vCovars[[j]], "nonMissingCount")
            group[[i]]$vCovars[[j]] <- group[[i]]$vCovars[[j]] +
                attr(group[[i]]$vCovars[[j]], "vartotal") /
                    attr(group[[i]]$vCovars[[j]], "nonMissingCount")
        }
        varmean <- vartotal / nonMissingCount
        for (i in 1:length(group))
        {
            j <- match(atts$vCovars[covar], names(group[[i]]$vCovars))
            if (is.na(j))
            {
                stop("inconsistent covariate names")
            }

            group[[i]]$vCovars[[j]] <- group[[i]]$vCovars[[j]] -
                varmean
        }
        simTotal <- 0
        simCnt <- 0
        anyMissing <- FALSE
        ## first find the range
        thisrange <- matrix(NA, ncol=length(group), nrow=2)
        values <- NULL
        for (i in 1:length(group))
        {

            j <- match(atts$vCovars[covar], names(group[[i]]$vCovars))
            if (is.na(j))
            {
                stop("inconsistent covariate names")
            }
            thisrange[, i] <- range(group[[i]]$vCovars[[j]],
                                    na.rm=TRUE)
            if (any(is.na(group[[i]]$vCovars[[j]])))
            {
                anyMissing <- TRUE
            }
            values <- c(values, unique(group[[i]]$vCovars[[j]]))
        }
        rr <- range(thisrange, na.rm=TRUE)
        if (rr[2] == rr[1] && !anyMissing)
        {
            vCovarPoszvar[covar] <-  FALSE
        }
        vCovarRange[covar] <- rr[2] - rr[1]
        ##then calculate similarity. Note ignore final observation
        for (i in 1:length(group))
        {

            j <- match(atts$vCovars[covar], names(group[[i]]$vCovars))
            tmpmat <- group[[i]]$vCovars[[covar]]
            tmp <- rangeAndSimilarity(tmpmat, rr)
            simTotal <- simTotal + tmp$simTotal
            simCnt <- simCnt + tmp$simCnt
        }
        simMean <- simTotal/simCnt
        vCovarSim[covar] <- simMean
        ## storage.mode(attr(vCovars[[i]], 'range')) <- 'double'
        vCovarMean[covar] <- varmean
        vCovarMoreThan2[covar] <- length(unique(values)) > 2
    }
    dycCovarMean <- rep(NA, length(atts$dycCovar))
    dycCovarRange <- rep(NA, length(atts$dycCovar))
    dycCovarRange2 <- matrix(NA, length(atts$dycCovar), 2)
    for (covar in seq(along=atts$dycCovar))
    {
        if (length(group) > 1)
            stop("error in dyadic constant covariate, group create")
        j <- match(atts$dycCovars[covar], names(group[[1]]$dycCovars))
        if (is.na(j))
        {
            stop("inconsistent covariate names")
        }
        dycCovarMean[covar] <- attr(group[[1]]$dycCovars[[j]], "mean")
        dycCovarRange[covar] <- attr(group[[1]]$dycCovars[[j]], "range")
        dycCovarRange2[covar, ] <- attr(group[[1]]$dycCovars[[j]], "range2")
    }
    dyvCovarMean <- rep(NA, length(atts$dyvCovar))
    dyvCovarRange <- rep(NA, length(atts$dyvCovar))
    for (covar in seq(along=atts$dyvCovar))
    {
        vartotal <- 0
        nonMissingCount <- 0
        thisrange <- matrix(NA, ncol=length(group), nrow=2)
        for (i in 1:length(group))
        {
            j <- match(atts$dyvCovars[covar], names(group[[i]]$dyvCovars))
            if (is.na(j))
            {
                stop("inconsistent covariate names")
            }
            vartotal <- vartotal + sum(group[[i]]$dyvCovars[[j]], na.rm=TRUE)
            nonMissingCount <- nonMissingCount +
                sum(!is.na(group[[i]]$dyvCovars[[j]]))
            thisrange[, i] <- range(group[[i]]$dyvCovars[[j]],
                                    na.rm=TRUE)
        }
        dyvCovarMean[covar] <- vartotal / nonMissingCount
        rr <- range(thisrange, na.rm=TRUE)
        dyvCovarRange[covar] <- rr[2] - rr[1]
   }
    names(bSim) <- netnames
    names(bRange) <- netnames
    colnames(behRange) <- netnames
    names(vCovarRange) <- names(vCovarSim) <- atts$vCovars
    names(cCovarRange) <- names(cCovarSim) <- atts$cCovars
    names(dycCovarMean) <- atts$dycCovars
    names(dyvCovarMean) <- atts$dyvCovars
    attr(group, "bRange") <- bRange
    attr(group, "behRange") <- behRange
    attr(group, "bSim") <- bSim
    attr(group, "bPoszvar") <- bPoszvar
    attr(group, "bMorethan2") <- bMoreThan2
    attr(group, "bAnyMissing") <- bAnyMissing
    attr(group, "cCovarPoszvar") <- cCovarPoszvar
    attr(group, "cCovarMoreThan2") <- cCovarMoreThan2
    attr(group, "cCovarRange") <- cCovarRange
    attr(group, "cCovarRange2") <- cCovarRange2
    attr(group, "cCovarSim") <- cCovarSim
    attr(group, "cCovarMean") <- cCovarMean
    attr(group, "vCovarRange") <- vCovarRange
    attr(group, "vCovarSim") <- vCovarSim
    attr(group, "vCovarMoreThan2") <- vCovarMoreThan2
    attr(group, "vCovarPoszvar") <- vCovarPoszvar
    attr(group, "vCovarMean") <- vCovarMean
    attr(group, "dycCovarMean") <- dycCovarMean
    attr(group, "dycCovarRange") <- dycCovarRange
    attr(group, "dycCovarRange2") <- dycCovarRange2
    attr(group, "dyvCovarRange") <- dyvCovarRange
    attr(group, "dyvCovarMean") <- dyvCovarMean
    group
}
sienaGroupCreate <- function(objlist, singleOK=FALSE)
{
    copyAttributes <- function(x, y)
    {
        atts <- attributes(y)
        attr(x, "rangep") <- matrix(atts$range2, ncol=ncol(x), nrow=2)
        attr(x, "meanp") <- rep(atts$mean, ncol(x))
        attr(x, "range") <- atts$range
        attr(x, 'mean') <- atts$mean
        attr(x, 'vartotal') <- atts$vartotal
        attr(x, 'nonMissingCount') <- atts$nonMissingCount
        rr <- rangeAndSimilarity(x, atts$range2)
        attr(x, 'poszvar') <- atts$poszvar
        attr(x, 'similarity') <- rr$simValues
        attr(x, 'simMean') <- rr$simMean
        attr(x, 'moreThan2') <- atts$moreThan2
        attr(x, 'name') <- atts$name
        ## storage.mode(attr(vCovars[[i]], 'range')) <- 'double'
       x
    }
    if (!is.list(objlist))
    {
        stop('Need a list of objects')
    }
    if (any (sapply(objlist, function(x) !inherits(x, 'siena'))))
    {
        stop('Not a list of valid siena objects')
    }
    if (length(objlist) == 1 && !singleOK)
    {
        stop('Need more than one siena object')
    }
    ## get hold of the attributes from the networks and lists of periods etc.
    group <- objlist
    ## get the names of the nets and other objects

    netnames <-  names(objlist[[1]]$depvars)
    cCovars <- names(objlist[[1]]$cCovars)
    if (is.null(cCovars))
    {
        cCovars <- character(0)
    }
    vCovars <- names(objlist[[1]]$vCovars)
    if (is.null(vCovars))
    {
        vCovars <- character(0)
    }
    dycCovars <-  names(objlist[[1]]$dycCovars)
    if (is.null(dycCovars))
    {
        dycCovars <- character(0)
    }
    dyvCovars <- names(objlist[[1]]$dyvCovars)
    if (is.null(dyvCovars))
    {
        dyvCovars <- character(0)
    }
    nNetnames <- length(netnames)
    anyMissing <- rep(FALSE, nNetnames)
    symmetric <- rep(TRUE, nNetnames)
    structural <- rep(FALSE, nNetnames)
    allUpOnly <- rep(TRUE, nNetnames)
    allDownOnly <- rep(TRUE, nNetnames)
    anyUpOnly <- rep(FALSE, nNetnames)
    anyDownOnly <- rep(FALSE, nNetnames)
    types <- rep(NA, nNetnames)
    nodeSets <- vector('list', nNetnames)
    ccnodeSets <- rep(NA, length(cCovars))
    cvnodeSets <- rep(NA, length(vCovars))
    dycnodeSets <- vector("list", length(dycCovars))
    dyvnodeSets <- vector("list", length(dyvCovars))
    observations <- 0
    periodNos <- rep(NA, 2)
    for (i in 1:length(objlist))
    {
        for (j in 1:length(objlist[[i]]$depvars))
        {
            varname <- names(objlist[[i]]$depvars)[j]
            netnamesub <- match(varname, netnames)
            if (is.na(netnamesub))
            {
                stop('object names inconsistent')
            }
            attribs <- attributes(objlist[[i]]$depvars[[j]])
            if (is.na(types[netnamesub]))
            {
                types[netnamesub] <- attribs[['type']]
            }
            else if (types[netnamesub] != attribs[['type']])
            {
                stop('Inconsistent network types')
            }
            if (is.null(nodeSets[[netnamesub]]))
            {
                nodeSets[[netnamesub]] <- attribs[['nodeSet']]
            }
            else if (any(nodeSets[[netnamesub]] != attribs[['nodeSet']]))
            {
                stop('Inconsistent node Sets')
            }
            if (attribs[['type']] == 'oneMode')
            {
                if (!attribs[['symmetric']])
                    symmetric[netnamesub] <- FALSE
            }
            if (any(attribs[['uponly']]))
            {
                anyUpOnly[netnamesub] <- TRUE
            }
            if (any(attribs[['downonly']]))
            {
                anyDownOnly[netnamesub] <- TRUE
            }
            if (!all(attribs[['uponly']]))
            {
                allUpOnly[netnamesub] <- FALSE
            }
            if (!all(attribs[['downonly']]))
            {
                allDownOnly[netnamesub] <- FALSE
            }
            if (attribs[["missing"]])
            {
                anyMissing[netnamesub] <- TRUE
            }
            if (attribs[["structural"]])
            {
                structural[netnamesub] <- TRUE
            }
       }
        for (j in seq(along=objlist[[i]]$cCovars))
        {
            varname <- names(objlist[[i]]$cCovars)[j]
            covarsub <- match(varname, cCovars)
            if (is.na(covarsub))
            {
                stop('covariate names inconsistent')
            }
            attribs <- attributes(objlist[[i]]$cCovars[[j]])
            if (is.na(ccnodeSets[covarsub]))
            {
                ccnodeSets[covarsub] <- attribs[['nodeSet']]
            }
            else if (ccnodeSets[covarsub] != attribs[['nodeSet']])
            {
                stop('Inconsistent node Sets')
            }
        }
        for (j in seq(along=objlist[[i]]$vCovars))
        {
            varname <- names(objlist[[i]]$vCovars)[j]
            covarsub <- match(varname, vCovars)
            if (is.na(covarsub))
            {
                stop('covariate names inconsistent')
            }
            attribs <- attributes(objlist[[i]]$vCovars[[j]])
            if (is.na(cvnodeSets[covarsub]))
            {
                cvnodeSets[covarsub] <- attribs[['nodeSet']]
            }
            else if (cvnodeSets[covarsub] != attribs[['nodeSet']])
            {
                stop('Inconsistent node Sets')
            }
        }

         for (j in seq(along=objlist[[i]]$dycCovars))
        {
            varname <- names(objlist[[i]]$dycCovars)[j]
            covarsub <- match(varname, dycCovars)
            if (is.na(covarsub))
            {
                stop('covariate names inconsistent')
            }
            attribs <- attributes(objlist[[i]]$dycCovars[[j]])
            if (is.null(dycnodeSets[covarsub]))
            {
                dycnodeSets[[covarsub]] <- attribs[['nodeSet']]
            }
            else if (any(dycnodeSets[[covarsub]] != attribs[['nodeSet']]))
            {
                stop('Inconsistent node Sets')
            }
        }
        for (j in seq(along=objlist[[i]]$dyvCovars))
        {
            varname <- names(objlist[[i]]$dyvCovars)[j]
            covarsub <- match(varname, dyvCovars)
            if (is.na(covarsub))
            {
                stop('covariate names inconsistent')
            }
            attribs <- attributes(objlist[[i]]$dyvCovars[[j]])
            if (is.null(dyvnodeSets[[covarsub]]))
            {
                dyvnodeSets[[covarsub]] <- attribs[['nodeSet']]
            }
            else if (any(dyvnodeSets[[covarsub]] != attribs[['nodeSet']]))
            {
                stop('Inconsistent node Sets')
            }
        }
        newobs <- objlist[[i]]$observations
        periodNos[observations + (1 : (newobs - 1))] <-
                  observations + i - 1 + (1 : (newobs - 1))
        observations <- observations + objlist[[i]]$observations -1
    }
    ## if more than one object, now create the group proper
    if (length(objlist) > 1)
    {
        for (i in 1:length(objlist))
        {
            ## constant covars turn into changing ones, probably because the
            ## actors change. Change on the individual data object.
            const <- objlist[[i]]$cCovars
            vars <- objlist[[i]]$vCovars
            oldnames <- names(vars)
            nVCovar <- length(vars)
            for (j in seq(along=const))
            {
                newcovar <-
                    varCovar(matrix(const[[j]],
                                    ncol = objlist[[i]]$observations - 1,
                                    nrow=length(const[[j]])),
                             nodeSet=attr(const[[j]], "nodeSet"))
                newcovar <- copyAttributes(newcovar, const[[j]])
                nVCovar <- nVCovar + 1
                vars[[nVCovar]] <- newcovar
            }
            names(vars) <- c(oldnames, names(const))
            objlist[[i]]$vCovars <- vars
            objlist[[i]]$cCovars <- list()

            ## now the dyadic ones similarly.

            const <- objlist[[i]]$dycCovars
            vars <- objlist[[i]]$dyvCovars
            oldnames <- names(vars)
            nVCovar <- length(vars)
            for (j in seq(along=const))
            {
                dim3 <- objlist[[i]]$observations - 1
                newcovar <-
                    varDyadCovar(array(const[[j]], dim=c(dim(const[[j]]),
                                                   dim3)),
                                 nodeSet=attr(const[[j]], "nodeSet"))
                 attr(newcovar, "vartotal") <- attr(const[[j]], "vartotal")
                 attr(newcovar, "nonMissingCount") <-
                     attr(const[[j]], "nonMissingCount")
               attr(newcovar, "mean") <- attr(const[[j]], "mean")
                attr(newcovar, "range") <- attr(const[[j]], "range")
                attr(newcovar, "rangep") <- rep(attr(const[[j]], "range"),
                                                dim3)
                attr(newcovar, "range2") <- matrix(attr(const[[j]], "range2"),
                                                   ncol=dim3, nrow=2)
                attr(newcovar, 'name') <- attr(const, "name")
                nVCovar <- nVCovar + 1
                vars[[nVCovar]] <- newcovar
            }
            names(vars) <- c(oldnames, names(const))
            objlist[[i]]$dyvCovars <- vars
            objlist[[i]]$dycCovars <- list()
     }
        ## update the working lists as well
        cCovars <- names(objlist[[1]]$cCovars)
        dycCovars <- names(objlist[[1]]$dycCovars)
        vCovars <- names(objlist[[1]]$vCovars)
        dyvCovars <- names(objlist[[1]]$dyvCovars)
    }
    symmetric[types=='behavior'] <- NA
    names(symmetric) <- netnames
    group <-  objlist
    attr(group, 'netnames') <- netnames
    attr(group, 'symmetric') <- symmetric
    attr(group, 'structural') <- structural
    attr(group, 'allUpOnly') <- allUpOnly
    attr(group, 'allDownOnly') <- allDownOnly
    attr(group, 'anyUpOnly') <- anyUpOnly
    attr(group, 'anyDownOnly') <- anyDownOnly
    attr(group, 'types') <- types
    attr(group, 'observations') <- observations
    attr(group, 'periodNos') <- periodNos
    attr(group, 'netnodeSets') <- nodeSets
    attr(group, 'cCovars') <- cCovars
    attr(group, 'vCovars') <- vCovars
    attr(group, 'dycCovars') <- dycCovars
    attr(group, 'dyvCovars') <- dyvCovars
    attr(group, 'ccnodeSets') <- ccnodeSets
    attr(group, 'cvnodeSets') <- cvnodeSets
    attr(group, 'dycnodeSets') <- dycnodeSets
    attr(group, 'dyvnodeSets') <- dyvnodeSets
    attr(group, 'compositionChange') <-
        any( sapply(objlist, function(x) length(x$compositionChange)>0))
    ## in fact assume all the same - validate earlier!
    if (attr(group, 'compositionChange'))
    {
        exooptions <- sapply(objlist[[1]]$compositionChange, function(x)
                             attr(x, "ccOption"))
        names(exooptions) <- sapply(objlist[[1]]$compositionChange, function(x)
                                    attr(x, "nodeSet"))
        attr(group, 'exooptions') <- exooptions
    }
    else
    {
        attr(group, 'exooptions') <- character(0)
    }
    if (is.null(names(group)))
    {
        names(group) <- paste('Data', 1:length(group), sep="")
    }
    class(group)<- c("sienaGroup", "siena")
    balmeans <- calcBalmeanGroup (group)
    names(balmeans) <- netnames
    attr(group, "balmean") <- balmeans
    group <- groupRangeAndSimilarityAndMean(group)
    bAnyMissing <- attr(group, "bAnyMissing")
    attr(group, "anyMissing") <- anyMissing | bAnyMissing
    attr(group, "bAnyMissing") <- NULL
    tmp <- getGroupNetRanges(group)
    colnames(tmp) <- netnames
    attr(group, "netRanges") <- tmp
    ## copy the global attributes down to individual level where appropriate
    copyGroupAttributes(group, "depvars", "balmean", "balmean")
    copyGroupAttributes(group, "depvars", "symmetric", "symmetric")
    copyGroupAttributes(group, "depvars", "bSim", "simMean")
    copyGroupAttributes(group, "depvars", "bposzvar", "poszvar")
    copyGroupAttributes(group, "depvars", "brange", "range")
    copyGroupAttributes(group, "depvars", "behRange", "behRange")
    copyGroupAttributes(group, "depvars", "MoreThan2", "moreThan2")
    copyGroupAttributes(group, "depvars", "missing", "missing")
    copyGroupAttributes(group, "depvars", "structural", "structural")

    copyGroupAttributes(group, "vCovars", "vCovarSim", "simMean")
    copyGroupAttributes(group, "vCovars", "vCovarRange", "range", TRUE)
    copyGroupAttributes(group, "vCovars", "vCovarMean", "mean", TRUE)
    copyGroupAttributes(group, "vCovars", "vCovarPoszvar", "poszvar")
    copyGroupAttributes(group, "vCovars", "vCovarMoreThan2", "moreThan2")
    copyGroupAttributes(group, "vCovars", "vCovarRange", "range", TRUE)
    copyGroupAttributes(group, "dycCovars", "dycCovarMean", "mean")
    copyGroupAttributes(group, "dycCovars", "dycCovarRange2", "range2", TRUE)
    copyGroupAttributes(group, "dyvCovars", "dyvCovarMean", "mean")
    copyGroupAttributes(group, "dyvCovars", "dyvCovarRange", "range", TRUE)
    group
}
copyGroupAttributes <- function(group, vartype, groupattrname, attrname,
                                storage.mode=FALSE)
{
    attributeValue <- attr(group, groupattrname)
    for (i in names(group))
    {
        for (j in names(group[[i]][[vartype]]))
        {
            attr(group[[i]][[vartype]][[j]], attrname) <-
                as.vector(attributeValue[j])
            if (storage.mode)
                 storage.mode(attr(group[[i]][[vartype]][[j]], attrname)) <-
                     "double"
        }
    }
    group
}

calcBalmeanGroup <- function(data)
{
    atts <- attributes(data)
    netnames <- atts$netnames
    types <- atts$types
   ## cat(types,'\n')
    balmeans <- rep(NA, length(netnames))
    for (net in seq(along=netnames))
    {
        if (types[net] == "oneMode")
        {
            tempra <- 0
            temprb <- 0
            for (i in 1: length(data))
            {
                j <- match(netnames[net], names(data[[i]]$depvars))
                if (is.na(j))
                    stop("network names not consistent")
                depvar <- data[[i]]$depvars[[j]]
                dims <- attr(depvar, "netdims")
                for (k in 1 : (dims[3] - 1))
                {
                    ## remove structural values here?
                    if (attr(depvar, "sparse"))
                    {
                        tmp <- depvar[[k]]
                        diag(tmp)  <-  NA ## just in case
                        x1 <- tmp@x
                        struct <- !is.na(x1) | x1 %in% c(10, 11)
                        x1[struct] <- x1[struct] - 10
                        tmp@x <- x1
                        tempra <- tempra +
                            sum(2 * colSums(tmp, na.rm=TRUE) *
                                colSums(1 - tmp, na.rm=TRUE),
                                na.rm=TRUE)
                        tmp1 <- colSums(is.na(tmp))
                        tmp2 <- colSums(!is.na(tmp))
                        temprb <- temprb + sum(tmp2 * (tmp2 - 1))
                    }
                    else
                    {
                        tmp <- depvar[, , k]
                        diag(tmp) <- NA ## just in case
                        struct <- !is.na(tmp) & tmp %in% c(10, 11)
                        tmp[struct] <- tmp[struct] - 10
                        tempra <- tempra +
                            sum(2 * colSums(tmp, na.rm=TRUE) *
                                colSums(1 - tmp, na.rm=TRUE),
                                na.rm=TRUE)
                        tmp1 <- colSums(is.na(tmp))
                        tmp2 <- colSums(!is.na(tmp))
                        temprb <- temprb + sum(tmp2 * (tmp2 - 1))
                    }
                }
            }
            balmeans[net] <- tempra / temprb
        }
    }
    balmeans
}
calcBalmean <- function(depvar)
{
    tempra <- 0
    temprb <- 0
    dims <- attr(depvar, "netdims")
    for (k in 1 : (dims[3] - 1))
    {
        if (attr(depvar, "sparse"))
        {
            tmp <- depvar[[k]]
            diag(tmp)  <-  NA ## just in case
            x1 <- tmp@x
            struct <- !is.na(x1) | x1 %in% c(10, 11)
            x1[struct] <- x1[struct] - 10
            tmp@x <- x1
            tempra <- tempra +
                sum(2 * colSums(tmp, na.rm=TRUE) *
                    colSums(1 - tmp, na.rm=TRUE),
                    na.rm=TRUE)
            tmp1 <- colSums(is.na(tmp))
            tmp2 <- colSums(!is.na(tmp))
            temprb <- temprb + sum(tmp2 * (tmp2 - 1))
     }
        else
        {
            ##extract this observation
            tmp <- depvar[, , k]
            ##remove diagonal
            diag(tmp) <- NA ## just in case
            ##subtract 10 from structurals
            struct <- !is.na(tmp) & tmp %in% c(10, 11)
            tmp[struct] <- tmp[struct] - 10
            ## colSums here counts how many 1's or 0's in tmp
            tempra <- tempra +
                sum(2 * colSums(tmp, na.rm=TRUE) *
                    colSums(1 - tmp, na.rm=TRUE),
                    na.rm=TRUE)
            tmp1 <- colSums(is.na(tmp)) ## counts missings in column
            tmp2 <- colSums(!is.na(tmp)) ## counts non-missings in column
            temprb <- temprb + sum(tmp2 * (tmp2 - 1))
        }
    }
    balmean <- tempra / temprb
    #  cat(tempra, temprb,balmean,'\n')
    balmean
}
getGroupNetRanges <- function(data)
{
    atts <- attributes(data)
    netnames <- atts$netnames
    types <- atts$types
   ## cat(types,'\n')
    ranges <- matrix(NA, ncol=length(netnames), nrow=2)
    varmin <- NA
    varmax <- NA
    for (net in seq(along=netnames))
    {
        if (types[net] == "oneMode")
        {
            varmin <- NA
            varmax <- NA
            for (i in 1: length(data))
            {
                j <- match(netnames[net], names(data[[i]]$depvars))
                if (is.na(j))
                    stop("network names not consistent")
                depvar <- data[[i]]$depvars[[j]]
                ### need to exclude structural values
                if (attr(depvar, "sparse"))
                {
              ##   varmin <- min(varmin, sapply(depvar, function(x)
              ##                               {
              ##                                   tmp <- x@x
              ##                                   min(tmp[!(is.na(tmp) |
              ##                                         tmp %in% c(10, 11))])
              ##                               }), na.rm=TRUE)
                    varmin <- 0 ## in sparse matrices 0's are not there
                    varmax <- max(varmax, sapply(depvar, function(x)
                                             {
                                                 tmp <- x@x
                                                 max(tmp[!(is.na(tmp) |
                                                           tmp %in% c(10, 11))])
                                             }), na.rm=TRUE)
                }
                else
                {
                    varmin <- min(varmin, depvar[!(is.na(depvar) |
                                                   depvar %in% c(10,11))],
                                  na.rm=TRUE)
                    varmax <- max(varmax, depvar[!(is.na(depvar) |
                                                   depvar %in% c(10,11))],
                                  na.rm=TRUE)
                }
            }
            ranges[, net] <- c(varmin, varmax)
        }
    }
    ranges
}
