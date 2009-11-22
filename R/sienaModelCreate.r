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
sienaModelCreate<- function(fn=simstats0c, usesimstats0c=TRUE,
                        projname="Siena", MaxDegree=0, useStdInits=FALSE,
                        n3=1000, nsub=4, maxlike=FALSE, diag=TRUE,
                        condvarno=0, condname='',
                        firstg=0.2, cond=NA, findiff=FALSE,  seed=NULL)
{
    model <- NULL
    model$projname <- projname
    model$useStdInits <- useStdInits
    model$checktime <- TRUE
    model$n3 <- n3
    model$firstg <- firstg
    model$maxrat <- 1.0
    model$maxmaxrat <- 10.0
    model$FRAN <- fn
    model$maxlike <-  maxlike
    model$cconditional <- cond
    model$condvarno <-  condvarno
    model$condname <- condname
    model$FinDiff.method <-  findiff
    model$nsub <- nsub
    model$diag <- diag
    model$ModelType <- 1
    model$MaxDegree <- MaxDegree
    model$randomSeed <- seed
    if (deparse(substitute(fn)) == "simstats0c")
        model$simstats0c <- TRUE
    else
        model$simstats0c <- usesimstats0c
    class(model) <- "sienaModel"
    model
}

model.create <- sienaModelCreate
