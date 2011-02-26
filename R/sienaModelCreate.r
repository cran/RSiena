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

ModelTypeStrings <- c("Standard actor-oriented model",
                      "Forcing model",
                      "Initiative model",
                      "Pairwise forcing model",
                      "Pairwise mutual model",
                      "Pairwise joint model")

##@sienaModelCreate DataCreate
sienaModelCreate <-
    function(fn,
             projname="Siena", MaxDegree=0, useStdInits=FALSE,
             n3=1000, nsub=4, maxlike=FALSE, diag=!maxlike,
             condvarno=0, condname='',
             firstg=0.2, cond=NA, findiff=FALSE,  seed=NULL,
             pridg=0.05, prcdg=0.05, prper=0.2, pripr=0.3, prdpr=0.3,
             prirms=0.05, prdrms=0.05, maximumPermutationLength=40,
             minimumPermutationLength=2, initialPermutationLength=20,
             modelType=1, mult=4)
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
    model$modelType <- modelType
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
    model$mult <- mult
    class(model) <- "sienaModel"
    model
}

model.create <- sienaModelCreate
