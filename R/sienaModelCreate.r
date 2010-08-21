#/******************************************************************************
# * SIENA: Simulation Investigation for Empirical Network Analysis
# *
# * Web: http://www.stats.ox.ac.uk/~snidjers/siena
# *
# * File: sienaModelCreate.r
# *
# * Description: This module contains the function for creating model objects.
# *
# *****************************************************************************/

##@sienaModelCreate DataCreate
sienaModelCreate <-
    function(fn,
             projname="Siena", MaxDegree=0, useStdInits=FALSE,
             n3=1000, nsub=4, maxlike=FALSE, diag=!maxlike,
             condvarno=0, condname='',
             firstg=0.2, cond=NA, findiff=FALSE,  seed=NULL,
             pridg=0.1, prcdg=0.1, prper=0.3, pripr=0.25, prdpr=0.25,
             prirms=0.0, prdrms=0.0, maximumPermutationLength=40,
             minimumPermutationLength=2, initialPermutationLength=20)
{
    model <- NULL
    model$projname <- projname
    model$useStdInits <- useStdInits
    model$checktime <- TRUE
    model$n3 <- n3
    model$firstg <- firstg
    model$maxrat <- 1.0
    model$maxmaxrat <- 10.0
    model$maxlike <-  maxlike
    model$FRANname <- deparse(substitute(fn))
    if (maxlike)
    {
        if (missing(fn))
        {
            model$FRANname <- "maxlikec"
        }
        if (is.na(cond))
        {
            cond <- FALSE
        }
        if (cond)
        {
            stop("Conditional estimation is not possible with",
                  "maximum likelihood estimation")
        }
        if (findiff)
        {
            stop("Finite differences estimation of derivatives",
                 "is not possible with maximum likelihood estimation")
        }
    }
    else
    {
        if (missing(fn))
        {
            model$FRANname <- "simstats0c"
        }
    }
    model$cconditional <- cond
    if (!is.na(cond) && cond && condvarno == 0 && condname == "")
    {
        model$condvarno <-  1
        model$condname <- ""
    }
    else
    {
        model$condvarno <-  condvarno
        model$condname <- condname
    }
    model$FinDiff.method <-  findiff
    model$nsub <- nsub
    model$diag <- diag
    model$ModelType <- 1
    model$MaxDegree <- MaxDegree
    model$randomSeed <- seed
    model$pridg <- pridg
    model$prcdg <- prcdg
    model$prper <- prper
    model$pripr <- pripr
    model$prdpr <- prdpr
    model$prirms <- prirms
    model$prdrms <- prdrms
    model$maximumPermutationLength <- maximumPermutationLength
    model$minimumPermutationLength <- minimumPermutationLength
    model$initialPermutationLength <- initialPermutationLength
    class(model) <- "sienaModel"
    model
}

model.create <- sienaModelCreate
