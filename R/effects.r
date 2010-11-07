#/******************************************************************************
# * SIENA: Simulation Investigation for Empirical Network Analysis
# *
# * Web: http://www.stats.ox.ac.uk/~snidjers/siena
# *
# * File: effects.r
# *
# * Description: This module contains the code for the creation of the
# * effects object to go with a Siena data object or group object.
# *****************************************************************************/
##@getEffects DataCreate
getEffects<- function(x, nintn = 10, behNintn=4, getDocumentation=FALSE)
{
    ##@duplicateDataFrameRow internal getEffects Put period numbers in
    duplicateDataFrameRow <- function(x, n)
    {
        tmp <- NULL
        for (i in 1:n)
        {
            xx <- x
            xx[, c("effectName", "functionName", "period")] <-
                sub("nnnnnn", periodNos[i], xx[, c("effectName", "functionName",
                                        "period")])
            tmp <-  rbind(tmp, xx)
        }
        tmp
    }
    ##@substituteNames internal getEffects replace xxxxxx, yyyyyy, zzzzzz
    substituteNames <- function(nameVectors, xName=NULL, yName=NULL, zName=NULL)
    {
        effects <- nameVectors[, c("effectName", "functionName",
                                   "interaction1", "interaction2")]
        if (!is.null(xName))
        {
            effects <- sapply(effects, function(x)
                                  gsub("xxxxxx", xName, x))
        }
        if (!is.null(yName))
        {
            effects <- sapply(effects, function(x)
                                  gsub("yyyyyy", yName, x))
        }
         if (!is.null(zName))
        {
            effects <- sapply(effects, function(x)
                                  gsub("zzzzzz", zName, x))
        }
        nameVectors[, c("effectName", "functionName",
                        "interaction1", "interaction2")] <- effects
        nameVectors
    }
    ##@createEffects internal getEffects Extract required rows and change text
    createEffects <- function(effectGroup, xName=NULL, yName=NULL)
    {
        effects <- allEffects[allEffects$effectGroup == effectGroup, ]
        if (nrow(effects) == 0)
        {
            stop("empty effect group")
        }
        if (any(is.na(effects$effectName)))
        {
            stop("missing effect name")
        }
        effects <- substituteNames(effects, xName, yName)
        effects
    }
    ##@networkRateEffects internal getEffects create a set of rate effects
    networkRateEffects <- function(depvar, varname, symmetric, bipartite)
    {
        if (symmetric)
        {
            rateEffects <- createEffects("symmetricRate", varname)
        }
        else if (bipartite)
        {
            rateEffects <- createEffects("bipartiteRate", varname)
        }
        else
        {
            rateEffects <- createEffects("nonSymmetricRate", varname)
        }
        if (observations == 1)
        {
            rateEffects <- rateEffects[-2, ] ## remove the extra period
        }
        else
        {
            ## get correct number of rows
            rateEffects <- rbind(duplicateDataFrameRow(rateEffects[2, ],
                                                       observations),
                                 rateEffects[-c(1, 2), ])
        }
        rateEffects
    }
    ##@oneModeNet internal getEffects
    oneModeNet <- function(depvar, varname)
    {
        symmetric <- attr(depvar, "symmetric")
        nodeSet <- attr(depvar, 'nodeSet')

        rateEffects <- networkRateEffects(depvar, varname, symmetric=symmetric,
                                          bipartite=FALSE)

        if (symmetric)
        {
            objEffects <- createEffects("symmetricObjective", varname)
        }
        else
        {
            objEffects <- createEffects("nonSymmetricObjective", varname)
        }
        for (j in seq(along = xx$dycCovars))
        {
            if (attr(xx$dycCovars[[j]], 'nodeSet')[1] == nodeSet)
            {
                objEffects <- rbind(objEffects,
                                    createEffects("dyadObjective",
                                                  names(xx$dycCovars)[j]))
            }
        }
        for (j in seq(along = xx$dyvCovars))
        {
            if (attr(xx$dyvCovars[[j]], 'nodeSet')[1] == nodeSet)
            {
                objEffects <- rbind(objEffects,
                                    createEffects("dyadObjective",
                                                  names(xx$dyvCovars)[j]))
            }
        }
        for (j in seq(along = xx$cCovars))
        {
            if (attr(xx$cCovars[[j]], 'nodeSet') == nodeSet)
            {
                tmp <- covarOneModeEff(names(xx$cCovars)[j],
                                     attr(xx$cCovars[[j]], 'poszvar'),
                                     attr(xx$cCovars[[j]], 'moreThan2'),
                                   symmetric)
                objEffects <-  rbind(objEffects, tmp$objEff)
                rateEffects <- rbind(rateEffects, tmp$rateEff)
            }
        }
        for (j in seq(along=xx$depvars))
        {
            if (types[j] == 'behavior' &&
                attr(xx$depvars[[j]], 'nodeSet') == nodeSet)
            {
                tmp <- covarOneModeEff(names(xx$depvars)[j],
                                     poszvar=TRUE,
                                     attr(xx$depvars[[j]], 'moreThan2'),
                                   symmetric)
                objEffects <- rbind(objEffects, tmp$objEff)
                rateEffects <- rbind(rateEffects, tmp$rateEff)
            }
        }
        for (j in seq(along=xx$vCovars))
        {
            if (attr(xx$vCovars[[j]], 'nodeSet') == nodeSet)
            {
                tmp <- covarOneModeEff(names(xx$vCovars)[j],
                                     attr(xx$vCovars[[j]], 'poszvar'),
                                     attr(xx$vCovars[[j]], 'moreThan2'),
                                   symmetric)
                objEffects <- rbind(objEffects,tmp$objEff)
                rateEffects<- rbind(rateEffects,tmp$rateEff)
            }
        }

### not sure we need this: if so then check relevant combinations of nodesets
        if (length(xx$cCovars) + length(xx$vCovars) +
            length(xx$dycCovars) + length(xx$dyvCovars) +
            length(types=='behavior') > 0)
        {
            interaction <- createEffects("unspecifiedNetInteraction")
            objEffects <-  rbind(objEffects, interaction[rep(1, nintn), ])
        }

        for (j in seq(along=xx$depvars))
        {
            otherName <- names(xx$depvars)[j]
            if (types[j] == 'oneMode' &&
                attr(xx$depvars[[j]], 'nodeSet') == nodeSet &&
                varname != otherName)
            {
                if (attr(xx$depvars[[j]], "symmetric"))
                {
                    objEffects <-
                        rbind(objEffects,
                              createEffects("nonSymmetricSymmetricObjective",
                                            otherName))
                }
                else
                {
                    objEffects <-
                        rbind(objEffects,
                              createEffects("nonSymmetricNonSymmetricObjective",
                                            otherName))
                }
            }
            if (types[j] == 'bipartite' &&
               (attr(xx$depvars[[j]], 'nodeSet')[1] == nodeSet))
            {
                objEffects <-
                    rbind(objEffects,
                          createEffects("nonSymmetricBipartiteObjective",
                                        otherName))
            }
            if (types[j] != "behavior" && varname != otherName)
            {
                for (k in seq(along=xx$cCovars))
                {
                    if (attr(xx$cCovars[[k]], 'nodeSet') == nodeSet)
                    {
                        objEffects <-
                            rbind(objEffects,
                                  createEffects("covarNetNetObjective",
                                                otherName, names(xx$cCovars)[k]))
                    }
                }
                 for (k in seq(along=xx$vCovars))
                {
                    if (attr(xx$vCovars[[k]], 'nodeSet') == nodeSet)
                    {
                        objEffects <-
                            rbind(objEffects,
                                  createEffects("covarNetNetObjective",
                                                otherName, names(xx$vCovars)[k]))
                    }
                }
                  for (k in seq(along=xx$depvars))
                {
                    if (types[j] == 'behavior' &&
                        attr(xx$depvars[[j]], 'nodeSet') == nodeSet)
                    {
                        objEffects <-
                            rbind(objEffects,
                                  createEffects("covarNetNetObjective",
                                                otherName,
                                                names(xx$depvars)[k]))
                    }
                }

            }
        }
        if ((nOneModes + nBipartites) > 1) ## add the network name
        {
            objEffects$functionName <- paste(varname, ': ',
                                             objEffects$functionName, sep = '')
            objEffects$effectName <- paste(varname, ': ',
                                           objEffects$effectName, sep = '')
        }
        ## now create the real effects, extra rows for endowment effects etc
        objEffects <- createObjEffectList(objEffects, varname)
        rateEffects <- createRateEffectList(rateEffects, varname)

        ## replace the text for endowment effects
        tmp <- objEffects$functionName[objEffects$type =='endow']
        tmp <- paste('Lost ties:', tmp)
        objEffects$functionName[objEffects$type == 'endow'] <- tmp

        ## get starting values
        starts <- getNetworkStartingVals(depvar)
        ##set defaults
        rateEffects[1:noPeriods, "include"] <- TRUE
        rateEffects[1:noPeriods, "initialValue"] <-  starts$startRate
        rateEffects$basicRate[1:observations] <- TRUE

        objEffects$untrimmedValue <- rep(0, nrow(objEffects))
        if (attr(depvar,'symmetric'))
        {
            objEffects[objEffects$shortName == "density" &
                       objEffects$type == "eval",
                       c('include', "initialValue", "untrimmedValue")] <-
                           list(TRUE, starts$degree, starts$untrimmed)
            objEffects[objEffects$shortName=='transTriads' &
                       objEffects$type=='eval','include'] <- TRUE
        }
        else
        {
            if (!(attr(depvar,'allUpOnly') || attr(depvar, 'allDownOnly')))
            {
                objEffects[objEffects$shortName == "density" &
                           objEffects$type == 'eval',
                           c('include', "initialValue", "untrimmedValue")] <-
                               list(TRUE, starts$degree, starts$untrimmed)
            }
            else
            {
                objEffects <-
                    objEffects[!objEffects$shortName == "density", ]
            }
            objEffects[objEffects$shortName == 'recip'&
                       objEffects$type == 'eval', 'include'] <- TRUE
        }
        rateEffects$basicRate[1:observations] <- TRUE
        list(effects=rbind(rateEffects = rateEffects, objEffects = objEffects),
             starts=starts)
    }

    ##@behaviornet internal getEffects
    behaviorNet <- function(depvar, varname)
    {
        nodeSet <- attr(depvar,'nodeSet')

        rateEffects <- createEffects("behaviorRate", varname)
        if (observations == 1)
        {
            rateEffects <- rateEffects[-2, ] ## remove the extra period
        }
        else
        {
            ## get correct number of rows
            rateEffects <- rbind(duplicateDataFrameRow(rateEffects[2, ],
                                                       observations),
                                 rateEffects[-c(1, 2), ])
        }

        objEffects <- createEffects("behaviorObjective", varname)

        for (j in seq(along=xx$depvars))
        {
            if (types[j] == 'oneMode' &&
                attr(xx$depvars[[j]], 'nodeSet') == nodeSet)
            {
                objEffects <- rbind(objEffects,
                                    createEffects("behaviorOneModeObjective",
                                               varname, names(xx$depvars)[j]))
                rateEffects <- rbind(rateEffects,
                                        createEffects("behaviorOneModeRate",
                                               varname, names(xx$depvars)[j]))
            }
            if (types[j] == 'bipartite' &&
                (attr(xx$depvars[[j]], 'nodeSet')[1] == nodeSet))
            {
                 objEffects <- rbind(objEffects,
                                    createEffects("behaviorBipartiteObjective",
                                               varname, names(xx$depvars)[j]))
                rateEffects <- rbind(rateEffects,
                                        createEffects("behaviorBipartiteRate",
                                               varname, names(xx$depvars)[j]))
            }
        }

        for (j in seq(along = xx$cCovars))
        {
            if (attr(xx$cCovars[[j]], 'nodeSet') == nodeSet)
            {
                tmp <- covBehEff(varname, names(xx$cCovars)[j], nodeSet,
                                 type='')
                objEffects<- rbind(objEffects, tmp$objEff)
                rateEffects<- rbind(rateEffects, tmp$rateEff)
           }
        }
        for (j in seq(along=xx$depvars))
        {
            if (types[j] == 'behavior' &&
                attr(xx$depvars[[j]], 'nodeSet') == nodeSet)
            {
                tmp <- covBehEff(varname, names(xx$depvars)[j], nodeSet, j==i,
                                 type='Beh')
                objEffects<- rbind(objEffects, tmp$objEff)
                rateEffects<- rbind(rateEffects, tmp$rateEff)
           }
        }
        for (j in seq(along=xx$vCovars))
        {
            if (attr(xx$vCovars[[j]], 'nodeSet') == nodeSet)
            {
                tmp <- covBehEff(varname, names(xx$vCovars)[j], nodeSet,
                                 type='Var')
                objEffects<- rbind(objEffects, tmp$objEff)
                rateEffects<- rbind(rateEffects, tmp$rateEff)
            }
        }
        for (j in seq(along=xx$depvars))
        {
            if (types[j] == 'oneMode' &&
                attr(xx$depvars[[j]], 'nodeSet') == nodeSet)
            {
                 objEffects <- rbind(objEffects,
                                    createEffects("behaviorOneModeObjective2",
                                               varname, names(xx$depvars)[j]))
            }
            if (types[j] == 'bipartite' &&
                attr(xx$depvars[[j]], 'nodeSet')[1] == nodeSet)
            {
                 objEffects <- rbind(objEffects,
                                    createEffects("behaviorBipartiteObjective2",
                                               varname, names(xx$depvars)[j]))
            }
        }
        interaction <- createEffects("unspecifiedBehaviorInteraction",
                                     varname)
        objEffects <- rbind(objEffects, interaction[rep(1, behNintn),])

        ## now create the real effects, extra rows for endowment effects etc
        objEffects <- createObjEffectList(objEffects, varname)
        rateEffects <- createRateEffectList(rateEffects, varname)

        ## get starting values
        starts <- getBehaviorStartingVals(depvar)
        ## set defaults
        if (!(attr(depvar,'allUpOnly') || attr(depvar, 'allDownOnly')))
        {
            objEffects[grepl("linear shape", objEffects$effectName) &
                       objEffects$type == 'eval',
                       c('include', 'initialValue','untrimmedValue')]  <-
                           list(TRUE, starts$tendency, starts$untrimmed)
        }
        else
        {
            objEffects <- objEffects[objEffects$shortName != "linear", ]
        }
        if (attr(depvar, "range") > 2)
        {
            objEffects[grepl("quadratic shape", objEffects$effectName) &
                       objEffects$type == 'eval','include']  <- TRUE
            ## no starting value for quadratic effect
        }


        rateEffects[1:observations, 'include'] <- TRUE
        rateEffects[1:noPeriods, 'initialValue'] <-  starts$startRate
        rateEffects$basicRate[1:observations] <- TRUE

        ## alter the text for endowment names
        objEffects$effectName[objEffects$type == 'endow'] <-
            sub('behavior', 'dec. beh.',
                objEffects$effectName[objEffects$type == 'endow'])

        list(effects = rbind(rateEffects = rateEffects,
             objEffects = objEffects), starts=starts)
    }

    ##@bipartiteNet internal getEffects
    bipartiteNet <- function(depvar, varname)
    {
        nodeSets <- attr(depvar, 'nodeSet')

        rateEffects <- networkRateEffects(depvar, varname, symmetric=FALSE,
                                          bipartite=TRUE)

        objEffects <- createEffects("bipartiteObjective", varname)

        for (j in seq(along = xx$dycCovars))
        {
            if (all(nodeSets == attr(xx$dycCovars[[j]], 'nodeSet')))
            {
                objEffects <- rbind(objEffects,
                                    createEffects("dyadBipartiteObjective",
                                                  names(xx$dycCovars)[j] ))
            }
        }
        for (j in seq(along = xx$dyvCovars))
        {
            if (all(nodeSets == attr(xx$dyvCovars[[j]], 'nodeSet')))
            {
                objEffects <- rbind(objEffects,
                                    createEffects("dyadBipartiteObjective",
                                                  names(xx$dyvCovars)[j]))
            }
        }
        for (j in seq(along = xx$cCovars))
        {
            covNodeset <- match(attr(xx$cCovars[[j]], "nodeSet"),
                                nodeSets)
            if (!is.na(covNodeset))
            {
                tmp <- covarBipartiteEff(names(xx$cCovars)[j],
                                      attr(xx$cCovars[[j]],
                                           'poszvar'),
                                      attr(xx$cCovars[[j]],
                                           'moreThan2'),
                                      covNodeset)
                objEffects <- rbind(objEffects, tmp$objEff)
                rateEffects <- rbind(rateEffects, tmp$rateEff)
            }
        }
        for (j in seq(along=xx$depvars))
        {
            if (types[j] == "behavior")
            {
                covNodeset <- match(attr(xx$depvars[[j]], "nodeSet"),
                                    nodeSets)
                if (!is.na(covNodeset))
                {
                    tmp <- covarBipartiteEff(names(xx$depvars)[j],
                                           poszvar=TRUE,
                                           attr(xx$depvars[[j]],
                                                'moreThan2'),
                                           covNodeset)
                    objEffects <- rbind(objEffects,tmp$objEff)
                    rateEffects <- rbind(rateEffects,tmp$rateEff)
                }
            }
        }
        for (j in seq(along=xx$vCovars))
        {
            covNodeset <- match(attr(xx$vCovars[[j]], "nodeSet"),
                                nodeSets)
            if (!is.na(covNodeset))
            {
                tmp <- covarBipartiteEff(names(xx$vCovars)[j],
                                        attr(xx$vCovars[[j]],
                                             'poszvar'),
                                        attr(xx$vCovars[[j]],
                                             'moreThan2'),
                                        covNodeset)
                objEffects <- rbind(objEffects, tmp$objEff)
                rateEffects <- rbind(rateEffects, tmp$rateEff)
            }
        }
### not sure we need this: if so then check relevant combinations of nodesets
        if (length(xx$cCovars) + length(xx$vCovars) +
            length(xx$dycCovars) + length(xx$dyvCovars) +
            length(types=='behavior') > 0)
        {
            interaction <- createEffects("unspecifiedNetInteraction")
            objEffects <-  rbind(objEffects, interaction[rep(1, nintn), ])
        }

         for (j in seq(along=xx$depvars))
        {
            otherName <- names(xx$depvars)[j]
            if (types[j] == 'oneMode' &&
                attr(xx$depvars[[j]], 'nodeSet') ==  nodeSets[1] )
            {
                if (attr(xx$depvars[[j]], "symmetric"))
                {
                    objEffects <-
                        rbind(objEffects,
                              createEffects("bipartiteSymmetricObjective",
                                            names(xx$depvars)[[j]]))
                }
                else
                {
                    objEffects <-
                        rbind(objEffects,
                              createEffects("bipartiteNonSymmetricObjective",
                                            names(xx$depvars)[[j]]))
                }
            }
            if (types[j] == 'bipartite' &&
                (attr(xx$depvars[[j]], 'nodeSet')[1] == nodeSets[1]) &&
                varname != otherName)
            {
                    objEffects <-
                        rbind(objEffects,
                              createEffects("bipartiteBipartiteObjective",
                                            names(xx$depvars)[[j]]))
            }
       }
        if ((nOneModes + nBipartites) > 1) ## add the network name
        {
            objEffects$functionName <- paste(varname, ': ',
                                             objEffects$functionName, sep = '')
            objEffects$effectName <- paste(varname, ': ',
                                           objEffects$effectName, sep = '')
        }
        ## now create the real effects, extra rows for endowment effects etc
        objEffects <- createObjEffectList(objEffects, varname)
        rateEffects <- createRateEffectList(rateEffects, varname)

        objEffects$functionName[objEffects$type == 'endow'] <-
            paste('Lost ties:',
                  objEffects$functionName[objEffects$type =='endow'])

        ## get starting values
        starts <- getBipartiteStartingVals(depvar)
        ##set defaults
        rateEffects[1:observations, 'include'] <- TRUE
        rateEffects[1:noPeriods, 'initialValue'] <-  starts$startRate
        rateEffects$basicRate[1:observations] <- TRUE

        if (!(attr(depvar,'allUpOnly') || attr(depvar, 'allDownOnly')))
        {
            objEffects[objEffects$shortName =='density' &
                       objEffects$type == 'eval',
                       c('include', 'initialValue', 'untrimmedValue')] <-
                           list(TRUE, starts$degree, starts$untrimmed)
        }
        else
        {
            objEffects <-
                objEffects[!objEffects$shortName == "density", ]
        }
        if (attr(xx$depvars[[i]],'uponly') || attr(xx$depvars[[i]],
                                                   'downonly'))
        {
            objEffects <-
                objEffects[!objEffects$shortName == "density", ]
        }

        rateEffects$basicRate[1:observations] <- TRUE

        list(effects=rbind(rateEffects = rateEffects, objEffects = objEffects),
             starts=starts)
    }

    ##@covarOneModeEff internal getEffects
    covarOneModeEff<- function(covarname, poszvar, moreThan2, symmetric)
    {
        if (symmetric)
        {
            covObjEffects <- createEffects("covarSymmetricObjective", covarname)
            covRateEffects <- createEffects("covarSymmetricRate", covarname)
        }
        else
        {
             covObjEffects <- createEffects("covarNonSymmetricObjective",
                                            covarname)
             covRateEffects <- createEffects("covarNonSymmetricRate", covarname)
         }

        if (!poszvar)
        {
            if (symmetric)
            {
                covObjEffects <- covObjEffects[1:2, ]
            }
            else
            {
                covObjEffects <- covObjEffects[3, ]
            }
        }
        if (!moreThan2)
        {
            covObjEffects <- covObjEffects[-2, ]
        }

        list(objEff=covObjEffects, rateEff=covRateEffects)
    }
    ##@covarBipartiteEff internal getEffects
    covarBipartiteEff<- function(covarname, poszvar, moreThan2, nodesetNbr)
    {
        covRateEffects  <-  NULL
        if (nodesetNbr == 1)
        {
            covObjEffects <-
                createEffects("covarBipartiteObjective", covarname)[3, ]
            covRateEffects <- createEffects("covarBipartiteRate", covarname)
        }
        else if (poszvar)
        {
            covObjEffects <- createEffects("covarBipartiteObjective", covarname)
            covObjEffects <- covObjEffects[1:2, ]
            if (!moreThan2)
            {
                covObjEffects <- covObjEffects[-2, ]
            }
        }
        else
        {
            covObjEffects <- NULL
        }


        list(objEff=covObjEffects, rateEff=covRateEffects)
    }
   ##@covBehEff internal getEffects
    covBehEff<- function(varname, covarname, nodeSet, same=FALSE,
        ## same indicates that varname and covarname are the same:
        ## just one rate effect required
                           type=c('', 'Var', 'Beh'))
    {
        covarBehObjEffects <- createEffects("covarBehaviorObjective", varname,
                                            covarname)
        covObjEffects <-  NULL
        if (!same)
        {
            covObjEffects<- covarBehObjEffects[1, ]
        }

        covRateEffects <- createEffects("covarBehaviorRate", varname, covarname)

        if (!same)
        {
            for (j in seq(along=xx$depvars))
            {
                if (types[j] == 'oneMode' &&
                    attr(xx$depvars[[j]], 'nodeSet') == nodeSet)
                {
                    covObjEffects <-
                        rbind(covObjEffects,
                              substituteNames(covarBehObjEffects[2, ],
                                              zName=names(xx$depvars)[j]))
                }
                if ((types[j] =="oneMode" &&
                    attr(xx$depvars[[j]], 'nodeSet') == nodeSet)
                || (types[j] == "bipartite" &&
                    attr(xx$depvars[[j]], 'nodeSet')[2] == nodeSet))
                {
                    covObjEffects <-
                        rbind(covObjEffects,
                              substituteNames(covarBehObjEffects[3, ],
                                              zName=names(xx$depvars)[j]))
                }
            }
        }
     #   if (!is.null(covObjEffects))
     #   {
     #       usestr <- paste("effFrom", type, sep="")
     #       covObjEffects$shortName <-
     #           sub("effFrom", usestr, covObjEffects$shortName)
     #   }

        list(objEff=covObjEffects, rateEff=covRateEffects)
    }
    ##@covarNetNetEff internal getEffects
    covarNetNetEff<- function(othernetname,
                              covarname, poszvar, moreThan2)
    {
        if (poszvar && (!moreThan2))
        {
            createEffects("covarNetNetObjective", othernetname,
                          covarname)
        }
    }
  ##@createObjEffectList internal getEffects
    createObjEffectList<- function(effects, varname)
    {
        nn <- nrow(effects)
        effectsname <- rep(varname, nn)
        effects <- data.frame(name=effectsname, effects, stringsAsFactors=FALSE)
        neweffects <- effects[rep(1:nn,
                                  times=(1 + as.numeric(effects$endowment))), ]
        neweffects$type <-  unlist(lapply(effects$endowment, function(x) if (x)
                               c('eval', 'endow') else 'eval'))
        nn <- nrow(neweffects)

        effectFn <- vector('list', nn)
        statisticFn <- vector('list', nn)
        neweffects$effectFn <- effectFn
        neweffects$statisticFn <- statisticFn
        neweffects
    }
    ##@createRateEffectList internal getEffects
    createRateEffectList<- function(effects, varname)
    {
        nn <- nrow(effects)
        effectsname <- rep(varname, nn)
        effects <- data.frame(name=effectsname, effects, stringsAsFactors=FALSE)
        effectFn <- vector('list', nn)
        statisticFn <- vector('list', nn)
        effects$effectFn <- effectFn
        effects$statisticFn <- statisticFn
        effects
    }
    ## start of function getEffects
    if (getDocumentation)
    {
        tt <- getInternals()
        return(tt)
    }
    if (!inherits(x, 'sienaGroup') && !inherits(x, 'siena'))
    {
        stop('Not a valid siena data object or group')
    }
    if (inherits(x, 'sienaGroup'))
    {
        groupx <- TRUE
    }
    else
    {
        groupx <- FALSE
    }
    ## validate the object?
    ## find the total number of periods to be processed = local var observations
    ## then process the first or only data object. Fill in starting values
    ## for other periods from the other objects, if any.
    if (groupx)
    {
        groupNames <- names(x)
        observations <- attr(x, 'observations')
        xx <- x[[1]]
        periodNos <- attr(x, "periodNos")
   }
    else
    {
        groupNames <- 'Group1'
        observations <- x$observations - 1
        xx <- x
        periodNos <- 1:observations
    }
    n <- length(xx$depvars)
    types <- sapply(xx$depvars, function(x)attr(x, 'type'))
    sparses <- sapply(xx$depvars, function(x)attr(x, 'sparse'))
   ## if (any(sparses))
   ##     require(Matrix)
    nOneModes <- sum(types == 'oneMode')
    nBehaviors <- sum(types == 'behavior')
    nBipartites <- sum(types =='bipartite')
    effects <- vector('list',n)
    nodeSetNames <- sapply(xx$nodeSets, function(x)attr(x, 'nodeSetName'))
    names(effects) <- names(xx$depvars)
    for (i in 1:n)
    {
        varname<- names(xx$depvars)[i]
        noPeriods <- xx$observations - 1 ## xx is a single object
        depvar <- xx$depvars[[i]]
        if (groupx)
        {
            netnamesub <- match(varname, attr(x, 'netnames'))
            attr(depvar, 'allUpOnly') <- attr(x, 'allUpOnly')[netnamesub]
            attr(depvar, 'anyUpOnly') <- attr(x, 'anyUpOnly')[netnamesub]
            attr(depvar, 'anyDownOnly') <- attr(x, 'anyDownOnly')[netnamesub]
            attr(depvar, 'allDownOnly') <- attr(x, 'allDownOnly')[netnamesub]
            if (types[i] == 'oneMode')
                attr(depvar, 'symmetric') <- attr(x, 'symmetric')[netnamesub]
        }
        else
        {
            attr(depvar, 'allUpOnly') <- all(attr(depvar, 'uponly'))
            attr(depvar, 'anyUpOnly') <- any(attr(depvar, 'uponly'))
            attr(depvar, 'anyDownOnly') <- any(attr(depvar, 'downonly'))
            attr(depvar, 'allDownOnly') <- all(attr(depvar, 'downonly'))
        }
        switch(types[i],
               behavior =
           {
               tmp <- behaviorNet(depvar, varname)
               effects[[i]] <- tmp$effects
               effects[[i]]$netType <- 'behavior'
               effects[[i]]$effectName <- as.character(effects[[i]]$effectName)
               attr(effects[[i]], 'starts') <- tmp$starts
           },
               oneMode =
           {
               tmp <- oneModeNet(depvar, varname)
               effects[[i]] <- tmp$effects
               effects[[i]]$netType <- 'oneMode'
               effects[[i]]$effectName <- as.character(effects[[i]]$effectName)
               attr(effects[[i]], 'starts') <- tmp$starts
           },
               bipartite =
           {
               tmp <- bipartiteNet(depvar, varname)
               effects[[i]] <- tmp$effects
               effects[[i]]$netType <- 'bipartite'
               effects[[i]]$effectName <- as.character(effects[[i]]$effectName)
               attr(effects[[i]], 'starts') <- tmp$starts

           },
               stop('error type'))
        effects[[i]]$groupName <- groupNames[1]
        effects[[i]]$group <- 1
    }
    ## add starting values for the other objects
    if (groupx && length(x) > 1)
    {
        period <-  xx$observations ##periods used so far

        for (group in 2:length(x))
        {
            xx <- x[[group]]
            n <- length(xx$depvars)
            types <- sapply(xx$depvars, function(x)attr(x, 'type'))
            noPeriods <- xx$observations - 1
            for (i in 1:n)
            {
                varname<- names(xx$depvars)[i]
                depvar <- xx$depvars[[i]]
                netnamesub <- match(varname, attr(x, 'netnames'))
                if (types[i] == 'oneMode')
                    attr(depvar, 'symmetric') <-
                        attr(x, 'symmetric')[netnamesub]
                switch(types[i],
                       behavior =
                   {
                       starts <-  getBehaviorStartingVals(depvar)
                       ## find the appropriate set of effects
                       eff <- match(varname, names(effects))
                       if (is.na(eff))
                           stop("depvars don't match")
                       effectname <- paste('rate ', varname,' (period ',
                                           period + 1:noPeriods,
                                           ')',sep='')
                       use <- effects[[eff]]$effectName %in%
                                      effectname
                       effects[[eff]][use, c('include','initialValue',
                                             'groupName', 'group', 'period')] <-
                                                 list(TRUE, starts$startRate,
                                                      groupNames[group], group,
                                                      1:noPeriods)
                       ## now sort out the tendency and update the
                       ## attribute on the effects list:
                       newdif <- c(starts$dif,
                                   attr(effects[[eff]], "starts")$dif)
                       meandif <- mean(newdif, na.rm=TRUE)
                       vardif <- var(as.vector(newdif), na.rm=TRUE)
                       if (meandif < 0.9 * vardif)
                       {
                           tendency <- 0.5 * log((meandif + vardif)/
                                                 (vardif - meandif))
                       }
                       else
                       {
                           tendency <- meandif / (vardif + 1)
                       }
                       untrimmed <- tendency
                       tendency <- ifelse(tendency < -3.0, -3.0,
                                          ifelse(tendency > 3/0, 3.0, tendency))
                       use <- (effects[[eff]]$shortName == "linear" &
                               effects[[eff]]$type == "eval")
                       effects[[eff]][use, c("include", "initialValue",
                                             "untrimmedValue")] <-
                                                 list(TRUE, tendency,
                                                      untrimmed)
                       attr(effects[[eff]], 'starts')$dif <- newdif
                   },
                       oneMode =
                   {
                       starts <- getNetworkStartingVals(depvar)
                       ## find the appropriate set of effects
                       eff <- match(varname, names(effects))
                       if (is.na(eff))
                       {
                           stop("depvars don't match")
                       }
                       effectname <- paste('constant ', varname,
                                           ' rate (period ',
                                           period + 1:noPeriods,')', sep='')
                       use <- effects[[eff]]$effectName %in% effectname
                       effects[[eff]][use, c('include', 'initialValue',
                                             "groupName", "group", "period")] <-
                                                 list(TRUE, starts$startRate,
                                                      groupNames[group], group,
                                                      1:noPeriods)
                       ## now sort out the degree and
                       ## update the attribute on the effects list
                       oldstarts <- attr(effects[[eff]], "starts")
                       alpha <- c(oldstarts$alpha, starts$alpha)
                       prec <-  c(oldstarts$prec, starts$prec)
                       degree <- sum(alpha * prec) / sum(prec)
                       untrimmed <- degree
                       degree <- ifelse (degree < -3, -3,
                                         ifelse(degree > 3, 3, degree))
                       attr(effects[[eff]], "starts")$alpha <- alpha
                       attr(effects[[eff]], "starts")$prec <-  prec
                       if (attr(depvar, 'symmetric'))
                       {
                           effects[[eff]][effects[[eff]]$shortName ==
                                          'density' &
                                          effects[[eff]]$type == 'eval',
                                          c('initialValue','untrimmedValue')] <-
                                              list(degree, untrimmed)
                       }
                       else
                       {
                           if (!(attr(x,'anyUpOnly') || attr(x, 'anyDownOnly')))
                           {
                               effects[[eff]][effects[[eff]]$shortName ==
                                              'density' &
                                              effects[[eff]]$type == 'eval',
                                              c('initialValue',
                                                "untrimmedValue")] <-
                                                    list(degree, untrimmed)
                           }
                       }
                       effects

                   },
                       bipartite =
                   {
                       starts <- getBipartiteStartingVals(depvar)
                       ## find the appropriate set of effects
                       eff <- match(varname, names(effects))
                       if (is.na(eff))
                       {
                           stop("depvars don't match")
                       }
                       effectname <- paste('constant ', varname,
                                           ' rate (period ',
                                           period + 1:noPeriods,')', sep='')
                       use <- effects[[eff]]$effectName %in% effectname
                       effects[[eff]][use, c('include', 'initialValue',
                                             'groupName', 'group',
                                             'period')] <-
                                                 list(TRUE, starts$startRate,
                                                      groupNames[group],
                                                      group, 1:noPeriods)
                    ## now sort out the degree and
                       ## update the attribute on the effects list
                       oldstarts <- attr(effects[[eff]], "starts")
                       alpha <- c(oldstarts$alpha, starts$alpha)
                       prec <-  c(oldstarts$prec, starts$prec)
                       degree <- sum(alpha * prec) / sum(prec)
                       untrimmed <- degree
                       degree <- ifelse (degree < -3, -3,
                                         ifelse(degree > 3, 3, degree))
                       attr(effects[[eff]], "starts")$alpha <- alpha
                       attr(effects[[eff]], "starts")$prec <-  prec
                       if (!(attr(x,'anyUpOnly') || attr(x, 'anyDownOnly')))
                       {
                           effects[[eff]][effects[[eff]]$shortName ==
                                          'density' &
                                          effects[[eff]]$type == 'eval',
                                          c('initialValue',
                                            "untrimmedValue")] <-
                                                list(degree, untrimmed)
                       }
                       effects

                   },
                       stop('error type'))
            }
            period <-  period + xx$observations ##periods used so far
        }
    }
    effects <- do.call(rbind, effects)
    effects$endowment <- NULL
    effects$effectGroup <- NULL
    attr(effects, "starts") <- NULL
    effects <- cbind(effects, effectNumber=1:nrow(effects))
    cl <- class(effects)
    if (groupx)
        class(effects) <- c('sienaGroupEffects','sienaEffects', cl)
    else
        class(effects) <- c('sienaEffects', cl)
    myrownames <- paste(sapply(strsplit(row.names(effects), ".", fixed=TRUE),
                               function(x)paste(x[1:2], collapse='.')),
                        effects$type, sep='.')
    myrownames <- paste(myrownames,
                         as.vector(unlist(sapply(table(myrownames),
                                                 function(x)1:x))), sep=".")
    myrownames <- sub("Effects", "", myrownames)

    effects
}
##@getBehaviorStartingVals DataCreate
getBehaviorStartingVals <- function(depvar)
{
    drange <- attr(depvar, 'range')
    depvar <- depvar[,1,]
    nactors <- nrow(depvar)
    ## rate
    if (round(drange) == 2)
    {
        rr <- range(depvar, na.rm=TRUE)
        dif <- t(diff(t(depvar))) ##calculate column differences
        tmp <- sapply(1:ncol(dif), function(x, y, z){
            mintab <- c("FALSE"=0, "TRUE" = 0)
            mintab[1] <- 1 + sum(dif[z[, x] == rr[1], x] == 0, na.rm=TRUE)
            mintab[2] <- 1 + sum(dif[z[, x] == rr[1], x] > 0, na.rm=TRUE)
            maxtab <- c("FALSE"=0, "TRUE" = 0)
            maxtab[1] <- 1 + sum(dif[z[, x] == rr[2], x] == 0, na.rm=TRUE)
            maxtab[2] <- 1 + sum(dif[z[, x] == rr[2], x] < 0, na.rm=TRUE)
            val <- mintab[2] / sum(mintab) + maxtab[2]  / sum(maxtab)
            if (val > 0.9) val <- 0.5
            c(-log(1 - val), mintab, maxtab)
        }, z = depvar, y = dif)
        startRate <- tmp[1, ]
        ##tendency
        tmp <- rowSums(tmp[-1, , drop=FALSE])
        tendency <- log((tmp[2] * (tmp[3] + tmp[4])) /
                        (tmp[4] * (tmp[1] + tmp[2])))
        untrimmed <- tendency
        tendency <- ifelse(tendency < -2, -2, ifelse (tendency > 2, 2,
                                                      tendency))
    }
    else
    {
        dif <- t(round(diff(t(depvar)))) ## get column differences
        absdif <- colSums(abs(dif), na.rm=TRUE) ##sum abs over columns
       # sqrdif <- colSums(dif*dif, na.rm=TRUE) ##sum sqr over columns
       # dif <- colSums(dif,na.rm=TRUE) ## sum over columns
       # startRate <- nactors / (nactors-1) * (sqrdif/nactors -
       #                                       dif * dif/nactors/nactors)
        startRate <- apply(dif, 2, var, na.rm=TRUE)
        startRate <- pmax(startRate, 0.1 + absdif/nactors)
        ## tendency
      ##  dif<- sum(dif/nactors,na.rm=TRUE)
      ##  dif2<- nactors/(nactors-1)*(sum(sqrdif,na.rm=TRUE)/nactors-dif^2)
        dif1 <- mean(dif, na.rm=TRUE)
        dif2 <- var(as.vector(dif), na.rm=TRUE)
        if (dif1 < 0.9 * dif2)
        {
            tendency <- 0.5 * log((dif1+dif2)/(dif2-dif1))
        }
        else
        {
            tendency <- dif1/(dif2+1)
        }
    }
    untrimmed <- tendency
    tendency <- ifelse(tendency < -3, -3,
                       ifelse (tendency > 3, 3, tendency))
    list(startRate=startRate, tendency=tendency, untrimmed = untrimmed, dif=dif)
}
##@getNetworkStartingVals DataCreate
getNetworkStartingVals <- function(depvar)
{
    noPeriods <- attr(depvar, "netdims")[3] - 1
    ##rate
    ##get distance and number of valid links
    if (!attr(depvar,'sparse'))
    {
        nactors <- nrow(depvar)
        use <- !is.na(depvar) & (depvar == 10 | depvar == 11)
        depvar[use] <- depvar[use] - 10  ## remove structural values
        tmp <- sapply(1:noPeriods, function(x, z){
            diag(z[ , , x]) <- NA
            diag(z[, , x + 1]) <- NA
            matdiff <- sum(z[, , x + 1] != z[, , x], na.rm=TRUE)
            matchange <- table(z[, , x + 1], z[, , x])
            matcnt <- nactors * nactors -
                sum(is.na(z[, , x + 1]) | is.na(z[, , x]))
            tmp <- c(matcnt=matcnt, matdiff=matdiff, matchange=matchange)
            names(tmp) <- c("matcnt", "matdiff", "matchangeFrom0To0",
                            "matchangeFrom0To1",
                            "matchangeFrom1To0", "matchangeFrom1To1")
            tmp
        }, z=depvar)
    }
    else
    {
        nactors <- nrow(depvar[[1]])
        matdiff<- rep(NA, noPeriods)
        matcnt<- rep(NA, noPeriods)
        matchange<- matrix(NA, nrow=4, ncol=noPeriods)
        for (i in 1: noPeriods)
        {
            mymat1 <- depvar[[i]]
            mymat2 <- depvar[[i+1]]
            use <- mymat1@x %in% c(10, 11)
            mymat1@x[use] <- mymat1@x[use] - 10
            use <- mymat2@x %in% c(10, 11)
            mymat2@x[use] <- mymat2@x[use] - 10
            mymat1 <- drop0(mymat1)
            mymat2 <- drop0(mymat2)
            diag(mymat1) <- NA
            diag(mymat2) <- NA
            mydif <- mymat2 - mymat1
            matdiff[i] <- sum(abs(mydif), na.rm=TRUE)
            tmp <- table(mydif@x)
            tmp00 <- nactors * nactors - length(mydif@x)
            tmp <- c(tmp00, tmp[c(3, 1, 2)])
            matchange[, i] <- tmp
            matcnt[i] <- sum(tmp)
        }
        matchange <- data.frame(matchange)
        tmp <- as.matrix(rbind(matcnt=matcnt, matdiff=matdiff,
                              matchange=matchange))
        row.names(tmp) <- c("matcnt", "matdiff", "matchangeFrom0To0",
                            "matchangeFrom0To1",
                            "matchangeFrom1To0", "matchangeFrom1To1")
    }
    distance <- attr(depvar, "distance" )
    if (attr(depvar,'symmetric'))
        startRate <- nactors * (0.2 + distance)/(tmp['matcnt',] %/% 2 +1)
    else
        startRate <- nactors * (0.2 + 2 * distance)/(tmp['matcnt',] + 1)
    startRate <- pmax(0.1, startRate)
    startRate <- pmin(100, startRate)
    ##degree
    matchange<- as.matrix(tmp[grep("matchange", rownames(tmp)),,drop=FALSE])
    if (attr(depvar,'symmetric'))
    {
        matchange <- matchange %/% 2
    }
    p01 <- ifelse (matchange["matchangeFrom0To0", ] +
                   matchange["matchangeFrom0To1", ] >=1,
                   matchange["matchangeFrom0To1", ] /
                   (matchange["matchangeFrom0To0", ] +
                    matchange["matchangeFrom0To1", ]), 0.5)
    p10 <- ifelse (matchange["matchangeFrom1To0", ]
                   + matchange["matchangeFrom1To1", ] >=1,
                   matchange["matchangeFrom1To0", ] /
                   (matchange["matchangeFrom1To0", ] +
                    matchange["matchangeFrom1To1", ]), 0.5)
    p01 <- pmax(0.02, p01)
    p10 <- pmax(0.02, p10)
    p01 <- pmin(0.98, p01)
    p10 <- pmin(0.98, p10)
    alpha <- 0.5 * log(p01 / p10)
    p00 <- ifelse (matchange["matchangeFrom0To0", ] +
                   matchange["matchangeFrom0To1", ] >=1,
                   matchange["matchangeFrom0To0", ] /
                   (matchange["matchangeFrom0To0", ] +
                    matchange["matchangeFrom0To1", ]), 0.0)
    p11 <- ifelse (matchange["matchangeFrom1To0", ]
                   + matchange["matchangeFrom1To1", ] >=1,
                   matchange["matchangeFrom1To1", ] /
                   (matchange["matchangeFrom1To0", ] +
                    matchange["matchangeFrom1To1", ]), 0.0)
    p00 <- pmax(0.02, p00)
    p11 <- pmax(0.02, p11)
    p00 <- pmin(0.98, p00)
    p11 <- pmin(0.98, p11)
    prec <- ifelse(matchange["matchangeFrom0To1", ] *
                   matchange["matchangeFrom1To0", ] >= 1,
                   4 / ((p00 / matchange["matchangeFrom0To1", ]) +
                       (p11 / matchange["matchangeFrom1To0", ])), 1e-6)
    alphaf1 <- sum(alpha * prec / sum(prec))
    untrimmed <- alphaf1
    alphaf1 <- ifelse(alphaf1 < -3, -3, ifelse(alphaf1 > 3, 3, alphaf1))
    list(startRate=startRate, degree=alphaf1, alpha=alpha, prec=prec, tmp=tmp,
        untrimmed = untrimmed)
}
##@getBipartiteStartingVals DataCreate
getBipartiteStartingVals <- function(depvar)
{
    noPeriods <- attr(depvar, "netdims")[3] - 1
    ##rate
    ##get distance and number of valid links
    if (!attr(depvar,'sparse'))
    {
        nsenders<- nrow(depvar)
        nreceivers <- ncol(depvar)
            use <- !is.na(depvar) & (depvar == 10 | depvar == 11)
            depvar[use] <- depvar[use] - 10  ## remove structural values
        tmp <- sapply(1:noPeriods, function(x, z){
            matdiff <- sum(z[, , x + 1] != z[, , x], na.rm=TRUE)
            matchange <- table(z[, , x + 1], z[, , x])
            matcnt <- nsenders * nreceivers -
                sum(is.na(z[, , x + 1]) | is.na(z[, , x]))
            tmp <- c(matcnt=matcnt, matdiff=matdiff, matchange=matchange)
            names(tmp) <- c("matcnt", "matdiff", "matchangeFrom0To0",
                            "matchangeFrom0To1",
                            "matchangeFrom1To0", "matchangeFrom1To1")
            tmp
        }, z=depvar)
    }
    else
    {
        nsenders <- nrow(depvar[[1]])
        nreceivers <- ncol(depvar[[2]])
        matdiff<- rep(NA, noPeriods)
        matcnt<- rep(NA, noPeriods)
        matchange<- matrix(NA, nrow=4, ncol=noPeriods)
        for (i in 1: noPeriods)
        {
            mymat1 <- depvar[[i]]
            mymat2 <- depvar[[i+1]]
            use <- mymat1@x %in% c(10, 11)
            mymat1@x[use] <- mymat1@x[use] - 10
            use <- mymat2@x %in% c(10, 11)
            mymat2@x[use] <- mymat2@x[use] - 10
            mymat1 <- drop0(mymat1)
            mymat2 <- drop0(mymat2)
            mydif <- mymat2 - mymat1
            matdiff[i] <- sum(abs(mydif), na.rm=TRUE)
            tmp <- table(mydif@x)
            tmp00 <- nsenders * nreceivers - length(mydif@x)
            tmp <- c(tmp00, tmp[c(3, 1, 2)])
            matchange[,i] <- tmp
            matcnt[i] <- sum(tmp)
        }
        matchange <- data.frame(matchange)
        tmp <-as.matrix(rbind(matcnt=matcnt, matdiff=matdiff,
                              matchange=matchange))
        row.names(tmp) <- c("matcnt", "matdiff", "matchangeFrom0To0",
                            "matchangeFrom0To1",
                            "matchangeFrom1To0", "matchangeFrom1To1")
    }
    distance <- attr(depvar, "distance" )
    startRate <- nsenders * (0.2 + 2 * distance)/(tmp['matcnt',] + 1)
    startRate <- pmax(0.1, startRate)
    startRate <- pmin(100, startRate)
    ##degree
    matchange<- as.matrix(tmp[grep("matchange", rownames(tmp)),,drop=FALSE])
    p01 <- ifelse (matchange["matchangeFrom0To0", ] +
                   matchange["matchangeFrom0To1", ] >=1,
                   matchange["matchangeFrom0To1", ] /
                   (matchange["matchangeFrom0To0", ] +
                    matchange["matchangeFrom0To1", ]), 0.5)
    p10 <- ifelse (matchange["matchangeFrom1To0", ]
                   + matchange["matchangeFrom1To1", ] >=1,
                   matchange["matchangeFrom1To0", ] /
                   (matchange["matchangeFrom1To0", ] +
                    matchange["matchangeFrom1To1", ]), 0.5)
    p01 <- pmax(0.02,p01)
    p10 <- pmax(0.02,p10)
    p01 <- pmin(0.98,p01)
    p10 <- pmin(0.98,p10)
    alpha <- 0.5 * log(p01/p10)
    p00 <- ifelse (matchange["matchangeFrom0To0", ] +
                   matchange["matchangeFrom0To1", ] >=1,
                   matchange["matchangeFrom0To0", ] /
                   (matchange["matchangeFrom0To0", ] +
                    matchange["matchangeFrom0To1", ]), 0.0)
    p11 <- ifelse (matchange["matchangeFrom1To0", ]
                   + matchange["matchangeFrom1To1", ] >=1,
                   matchange["matchangeFrom1To1", ] /
                   (matchange["matchangeFrom1To0", ] +
                    matchange["matchangeFrom1To1", ]), 0.0)
    p00 <- pmax(0.02,p00)
    p11 <- pmax(0.02,p11)
    p00 <- pmin(0.98,p00)
    p11 <- pmin(0.98,p11)
    prec <- ifelse(matchange["matchangeFrom0To1", ] *
                   matchange["matchangeFrom1To0", ] >= 1,
                   4 / ((p00 / matchange["matchangeFrom0To1", ]) +
                       (p11 / matchange["matchangeFrom1To0", ])), 1e-6)
    alphaf1 <- sum(alpha*prec/sum(prec))
    ## }
    untrimmed <- alphaf1
    alphaf1 <- ifelse(alphaf1 < -3, -3, ifelse(alphaf1 > 3, 3, alphaf1))
    list(startRate=startRate, degree=alphaf1, alpha=alpha, prec=prec, tmp=tmp,
        untrimmed = untrimmed)
}

