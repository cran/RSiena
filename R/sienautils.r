#/******************************************************************************
# * SIENA: Simulation Investigation for Empirical Network Analysis
# *
# * Web: http://www.stats.ox.ac.uk/~snidjers/siena
# *
# * File: sienautils.r
# *
# * Description: This module contains utilities for creating siena special
# * objects from vectors etc
# *****************************************************************************/

sienaCompositionChange <- function(changelist, nodeSet='Actors', option=1)
{
    if (!is.list(changelist))
        stop('changelist must be a list')
    warn <- getOption('warn')
    options(warn = -1)
    changelist <- lapply(changelist, as.numeric)
    options(warn = warn)
    if (any(sapply(changelist,function(x)any(is.na(x)))))
        stop('Non numeric data')
    out <- changelist
    class(out) <- 'compositionChange'
    if (is.vector(nodeSet) && length(nodeSet) > 1)
        stop ('only one node set relevant')
    if (!is.character(nodeSet))
        stop ('nodeset should be a character string')
    if (!is.numeric(option))
        stop("option should be numeric")
    if (option < 1 || option > 4)
        stop("option should be between 1 and 4")
    attr(out, 'nodeSet') <- nodeSet
    attr(out, "ccOption") <- option
    out
}
sienaCompositionChangeFromFile <- function(filename, nodeSet='Actors',
                                           fileobj=NULL, option=1)
{
    if (is.null(fileobj))
        tmp <- readLines(filename)
    else
        tmp <- fileobj
    changelist <- lapply(tmp, function(x)
                     {
                     #  cat("x1",x)
                         x <- sub("^ +", "", x)
                      #  cat("x2",x, "x")
                         x <- unlist(strsplit(x, " +"))
                     #   print(x)
                         x
                     })
    warn <- getOption('warn')
    options(warn = -1)
    changelist <- lapply(changelist, as.numeric)
    options(warn = warn)
    if (any(sapply(changelist,function(x)any(is.na(x)))))
        stop('Non numeric data')
    out <- changelist
    class(out) <- 'compositionChange'
    if (is.vector(nodeSet) && length(nodeSet) > 1)
        stop ('only one node set relevant')
    if (!is.character(nodeSet))
        stop ('nodeset should be a character string')
    if (!is.numeric(option))
        stop("option should be numeric")
    if (option < 1 || option > 4)
        stop("option should be between 1 and 4")
    attr(out, 'nodeSet') <- nodeSet
    attr(out, "ccOption") <- option
    out
}
sienaNodeSet <- function(n, nodeSetName='Actors', names=NULL)
{
    if (!is.numeric(n))
        stop('n must be numeric')
    if (!is.character(nodeSetName))
        stop('node set name must be a character string')
    if (!is.null(names) && length(names) != n)
        stop('node names incorrect length')
    out <- 1:n
    attributes(out) <- list(class='sienaNodeSet', nodeSetName=nodeSetName)
    if (!is.null(names))
        names(out) <- names
    out
}

coCovar <- function(val, nodeSet='Actors')
{
    ##vector, numeric or factor
    if (!is.vector(val))
        stop('val must be a vector')
    if (!(is.numeric(val) || is.factor(val)))
        stop('val must be numeric or a factor')
    out <- val
    class(out) <- 'coCovar'
    attr(out, 'nodeSet') <- nodeSet
    out
}

varCovar<- function(val, nodeSet='Actors')
{
    ##matrix, numeric or factor, nrow = nactors and cols = observations-1
    if (!is.matrix(val))
        stop('val must be a matrix')
    if (!(is.numeric(val) || is.factor(val)))
        stop('val must be numeric or a factor')
    out <- val
    class(out) <- 'varCovar'
    attr(out, 'nodeSet') <- nodeSet
    out
}

coDyadCovar<- function(val, nodeSets=c('Actors','Actors'))
{
    ##matrix, numeric or factor, dims= those of net - must validate later
    if (!is.matrix(val))
        stop('val must be a matrix')
    if (!(is.numeric(val) || is.factor(val)))
        stop('val must be numeric or a factor')
    if (length(nodeSets) > 2)
        stop('nodeSets may only have one or two elements')
     if (!is.character(nodeSets))
        stop('nodeSets must be a vector of character strings')
    out<- val
    class(out) <- 'coDyadCovar'
    attr(out, 'nodeSet') <- nodeSets
    attr(out, 'sparse') <- FALSE ### for now!
    out
}
varDyadCovar<- function(val, nodeSets=c('Actors','Actors'))
{
    ##array, numeric or factor, dims= those of net by observations-1 -
    ##must validate later
    if (!is.array(val) || !(length(dim(val)) == 3))
        stop('val must be a 3d array')
    if (!(is.numeric(val) || is.factor(val)))
        stop('val must be numeric or a factor')
    if (length(nodeSets) > 2)
        stop('nodeSets may only have one or two elements')
    if (!is.character(nodeSets))
        stop('nodeSets must be a vector of character strings')
    if (length(nodeSets) == 1)
        nodeSets <- c(nodeSets, nodeSets)
    out <- val
    class(out) <- 'varDyadCovar'
    attr(out, 'nodeSet') <- nodeSets
    attr(out, 'sparse') <- FALSE ### for now!
    out
}

sienaNet<- function(netarray, type=c('oneMode','bipartite','behavior'),
                    nodeSet='Actors', sparse=is.list(netarray))
{
    if (!sparse)
    {
        if (!is.array(netarray))
            stop('netarray must be an array or a list of sparse matrices')
        netdims <- dim(netarray)
        if (length(netdims) == 2) ## assume behavior network
        {
            dim(netarray) <- c(netdims[1], 1, netdims[2])
            netdims <- dim(netarray)
        }
        if (!is.numeric(netarray))
            stop('entries must be numeric')
    }
    else
    {
        require(Matrix)
        if (!is.list(netarray))
            stop('netarray must be an array or a list of sparse matrices')
        if (!all(sapply(netarray, function(x) is(x,'dgTMatrix'))))
            stop('not a list of sparse triples matrices')
        netdims <- sapply(netarray, dim) ## dimensions of network in columns
        if (any(netdims != netdims[, 1]))
            stop('all matrices must have the same dimension')
        netdims <- netdims[, 1]
        netdims[3] <- length(netarray)
    }
    observations <- netdims[3]
    if (observations < 2)
        stop('netarray must have at least two observations')
    if (netdims[2] == 1)
    {
        if (missing(type))
            type <- 'behavior'
        else
        {
            type <- match.arg(type)
            if (type != 'behavior')
                stop('incorrect type')
        }
    }
    else if (netdims[1] != netdims[2])
    {
        if (missing(type))
            type <- 'bipartite'
        else
        {
            type <- match.arg(type)
            if (type != 'bipartite')
                stop('incorrect type')
        }
    }
    else
    {
        type <- match.arg(type)
        if (type == 'behavior')
            stop('incorrect type')
    }
    if (type == 'bipartite')
    {
        if (!is.vector(nodeSet) || length(nodeSet) != 2)
            stop ('need 2 node sets for a bipartite network')
        if (!is.character(nodeSet[[1]]) || ! is.character(nodeSet[[2]]))
            stop ('nodesets should be character strings')
    }
    else
     {
        if (is.vector(nodeSet) && length(nodeSet) > 1)
            stop ('only one node set for this network')
        if (!is.character(nodeSet))
            stop ('nodeset should be a character string')
    }
    if (type != 'behavior')
        if (sparse)
        {
             if (!all(sapply(netarray, function(x)
                            {
                                tmp <- x@x
                                all(is.na(tmp) | tmp == 1 )
                            }
                                )))
                 stop('entries in networks must be 0 or 1')
         }
        else
        {
            if (!all(netarray %in% c(0,1) | is.na(netarray)))
                stop('entries in networks must be 0 or 1')
        }

    obj <- netarray
    class(obj) <- ('sienaNet')
    attr(obj, 'type') <- type
    attr(obj, 'sparse') <- sparse
    attr(obj, 'nodeSet') <- nodeSet
    attr(obj, 'netdims') <- netdims
    obj
}

validateSienaNet <- function(net)
{
    if (!inherits(net,'sienaNet'))
        stop ('Not a siena Net')
    if (!attr(net, 'type') %in% c('oneMode', 'bipartite', 'behavior'))
        stop ('invalid type in net')
}
