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
getEffects<- function(x, nintn = 10)
{
    oneModeNet <- function(depvar, varname)
    {
        nodeSet <- attr(depvar, 'nodeSet')
        if (attr(depvar, 'symmetric'))
        {
            if (observations > 1)
            {
                rateEffects <- paste('constant', varname,' rate (period ',
                                        periodNos, ')', sep = '')
                rateFunctions <- paste("Amount of network change in period",
                                       periodNos)
                rateShortNames <- rep('Rate', observations)
                ratePeriods <- 1:observations
                rateTypes <- rep(NA, observations)
            }
            else
            {
                rateEffects <- paste('basic rate parameter', varname)
                rateShortNames <- 'Rate'
                rateFunctions <- "Amount of network change"
                ratePeriods <- 1
                rateTypes <- NA
           }
            rateEffects <- c(rateEffects,
                             paste(symmetricRateEffects[-(1:2), 1], varname))
            rateFunctions <- c(rateFunctions, symmetricRateEffects[-(1:2), 2])
            rateShortNames <- symmetricRateEffects[, 3]
            ratePeriods <- c(ratePeriods, rep(NA, nrow(symmetricRateEffects)-2))
            rateTypes <- c(rateTypes, rep('structural',
                                          nrow(symmetricRateEffects)-2))
            objEffects <- symmetricObjEffects[, 1]
            objFunctions <- symmetricObjEffects[, 2]
            objEndowment <- symmetricObjEffects[, 3]
            objShortnames <- symmetricObjEffects[, 4]
            objParms <- symmetricObjEffects[, 5]
            objEffects <- createObjEffectList(objEffects, objFunctions,
                                              objEndowment, objShortNames,
                                              objParms, varname)
            rateEffects <- createRateEffectList(rateEffects, rateFunctions,
                                              rateShortNames, ratePeriods,
                                                rateTypes,
                                                varname)
            for (j in seq(along = xx$dycCovars))
            {
                if (attr(x$dycCovars[[j]], 'nodeSet')[1] == nodeSet)
                {
                    tmp <- dyadNetObjEff(names(xx$dycCovars)[j],
                                         symmetric=TRUE)
                    objEffects <- rbind(objEffects,   tmp$objEffects)
                }
            }
            for (j in seq(along = xx$dyvCovars))
            {
                if (attr(x$dvvCovars[[j]], 'nodeSet')[1] == nodeSet)
                {
                    tmp <- dyadNetObjEff(names(xx$dyvCovars)[j],
                                         symmetric = TRUE)
                    objEffects <- rbind(objEffects,   tmp$objEffects)
                }
            }
            for (j in seq(along = xx$cCovars))
            {
                if (attr(x$cCovars[[j]], 'nodeSet') == nodeSet)
                {
                    tmp <- covSymmNetEff(names(xx$cCovars)[j],
                                         attr(xx$cCovars[[j]], 'poszvar'),
                                         attr(xx$cCovars[[j]], 'moreThan2'))
                    objEffects <-  rbind(objEffects, tmp$objEff)
                    rateEffects <- rbind(rateEffects, tmp$rateEff)
               }
            }
            for (j in seq(along=x$depvars))
            {
                if (types[j] == 'behavior' &&
                    attr(x$depvars[[j]], 'nodeSet') == nodeSet)
                {
                    tmp <- covSymmNetEff(names(xx$depvars)[j],
                                        poszvar=TRUE,
                                        attr(xx$depvars[[j]], 'moreThan2'))
                    objEffects <- rbind(objEffects, tmp$objEff)
                    rateEffects <- rbind(rateEffects, tmp$rateEff)
               }
            }
            for (j in seq(along=x$vCovars))
            {
                if (attr(x$cCovars[[j]], 'nodeSet') == nodeSet)
                {
                    tmp <- covSymmNetEff(names(xx$vCovars)[j],
                                        attr(xx$vCovars[[j]], 'poszvar'),
                                        attr(xx$vCovars[[j]], 'moreThan2'))
                    objEffects <- rbind(objEffects,tmp$objEff)
                    rateEffects<- rbind(rateEffects,tmp$rateEff)
                }
            }
        }
        else ##not symmetric
        {
            if (observations > 1)
            {
                rateEffects <- paste('constant ', varname,' rate (period ',
                                     periodNos, ')', sep = '')
                rateFunctions <- paste("Amount of network change in period",
                                       periodNos)
                rateShortNames <- rep('Rate', observations)
                ratePeriods <- 1:observations
                rateTypes <- rep(NA, observations)
            }
            else
            {
                rateEffects <- paste('basic rate parameter', varname)
                rateFunctions <- "Amount of network change"
                rateShortNames <- 'Rate'
                ratePeriods <- 1
                rateTypes <- NA
            }
            rateEffects <- c(rateEffects, nonSymmetricRateEffects[-(1:2), 1])
            ratePeriods <- c(ratePeriods,
                             rep(NA, nrow(nonSymmetricRateEffects) - 2))
            rateTypes <- c(rateTypes, rep('structural',
                                          nrow(nonSymmetricRateEffects) - 2))
            objEffects <- nonSymmetricObjEffects[, 1]
            rateFunctions <- c(rateFunctions, nonSymmetricRateEffects[-(1:2),2])
            rateShortNames <- c(rateShortNames,
                                nonSymmetricRateEffects[-c(1:2), 3])
            objFunctions <- nonSymmetricObjEffects[, 2]
            objEndowment <- nonSymmetricObjEffects[, 3]
            objShortNames <- nonSymmetricObjEffects[, 4]
            objParms <- nonSymmetricObjEffects[, 5]

            objEffects <- createObjEffectList(objEffects, objFunctions,
                                              objEndowment, objShortNames,
                                              objParms, varname)
            rateEffects <- createRateEffectList(rateEffects, rateFunctions,
                                              rateShortNames, ratePeriods,
                                              rateTypes, varname)
            for (j in seq(along = xx$dycCovars))
            {
                if (attr(xx$dycCovars[[j]], 'nodeSet')[1] == nodeSet)
                {
                    tmp <- dyadNetObjEff(names(xx$dycCovars)[j],
                                         symmetric = FALSE)
                    objEffects <- rbind(objEffects, tmp$objEff)
                }
            }
            for (j in seq(along = xx$dyvCovars))
            {
                if (attr(xx$dyvCovars[[j]], 'nodeSet')[1] == nodeSet)
                {
                    tmp <- dyadNetObjEff(names(xx$dyvCovars)[j],
                                         symmetric = FALSE)
                    objEffects <- rbind(objEffects, tmp$objEff)
                }
            }
            for (j in seq(along = xx$cCovars))
            {
                if (attr(xx$cCovars[[j]], 'nodeSet') == nodeSet)
                {
                    tmp<- covNonSymmNetEff(names(xx$cCovars)[j],
                                           attr(xx$cCovars[[j]],
                                                'poszvar'),
                                           attr(xx$cCovars[[j]],
                                                'moreThan2'))
                    objEffects <- rbind(objEffects, tmp$objEff)
                    rateEffects <- rbind(rateEffects, tmp$rateEff)
                }
            }
            for (j in seq(along=xx$depvars))
            {
                if (types[j] == 'behavior' &&
                    attr(xx$depvars[[j]], 'nodeSet') == nodeSet)
                {
                    tmp <- covNonSymmNetEff(names(xx$depvars)[j],
                                            poszvar=TRUE,
                                            attr(xx$depvars[[j]],
                                                 'moreThan2'))
                    objEffects <- rbind(objEffects,tmp$objEff)
                    rateEffects <- rbind(rateEffects,tmp$rateEff)
                }
            }
            for (j in seq(along=xx$vCovars))
            {
                if (attr(xx$vCovars[[j]], 'nodeSet') == nodeSet)
                {
                    tmp <- covNonSymmNetEff(names(xx$vCovars)[j],
                                            attr(xx$vCovars[[j]],
                                                 'poszvar'),
                                            attr(xx$vCovars[[j]],
                                                 'moreThan2'))
                    objEffects <- rbind(objEffects, tmp$objEff)
                    rateEffects <- rbind(rateEffects, tmp$rateEff)
                }
            }
        }
### not sure we need this: if so then check relevant combinations of nodesets
       if (length(xx$cCovars) + length(xx$vCovars) +
            length(xx$dycCovars) + length(xx$dyvCovars) +
            length(types=='behavior') > 0)
        {
            objEff <- rep('unspecified interaction effect', nintn)
            objEnd <-  rep(TRUE, nintn)
            objFun <-  rep('unspecified interaction statistic', nintn)
            objSho <- rep('unspInt', nintn)
            objParms <- rep(0, nintn)
            objEffects <- rbind(objEffects, createObjEffectList(objEff, objFun,
                                              objEnd, objSho, objParms, varname))
        }
        if (nOneModes > 1)
        {
            rateEffects$functionName <- paste(varname, ': ',
                                              rateEffects$functionName,
                                              sep = '')
            objEffects$functionName <- paste(varname, ': ',
                                             objEffects$functionName, sep = '')
        }
        tmp <- objEffects$functionName[objEffects$type =='endow']
        tmp <- paste('Lost ties:', tmp)
        objEffects$functionName[objEffects$type == 'endow'] <- tmp
        starts <- getNetworkStartingVals(depvar)
        ##set defaults
        if (observations == 1)
            effectname <- paste('basic rate parameter', varname)
        else
            effectname <- paste('constant ', varname,' rate (period ',
                                1:noPeriods,')',sep='')
        rateEffects[rateEffects$effectName %in%
                    effectname, 'include'] <- TRUE
        rateEffects[rateEffects$effectName %in% effectname,
                    'initialValue'] <-  starts$startRate
        rateEffects$basicRate[1:observations] <- TRUE
        objEffects$untrimmedValue <- rep(0, nrow(objEffects))
        if (attr(depvar,'symmetric'))
        {
            objEffects[objEffects$effectName == 'degree (density)' &
                       objEffects$type == 'eval', 'include'] <- TRUE
            objEffects[objEffects$effectName =='degree (density)' &
                       objEffects$type == 'eval', 'initialValue'] <-
                           starts$degree
            objEffects[objEffects$effectName =='degree (density)' &
                       objEffects$type == 'eval', 'untrimmedValue'] <-
                           starts$untrimmed
            objEffects[objEffects$effectName=='transitive triads' &
                       objEffects$type=='eval','include'] <- TRUE
        }
        else
        {
            if (!(attr(depvar,'anyUpOnly') || attr(depvar, 'anyDownOnly')))
            {
                objEffects[objEffects$effectName =='outdegree (density)'&
                           objEffects$type == 'eval', 'include'] <- TRUE
                objEffects[objEffects$effectName ==
                           'outdegree (density)' &
                           objEffects$type == 'eval', 'initialValue'] <-
                               starts$degree
                                objEffects[objEffects$effectName ==
                           'outdegree (density)' &
                           objEffects$type == 'eval', 'untrimmedValue'] <-
                               starts$untrimmed

            }
            objEffects[objEffects$effectName == 'reciprocity'&
                       objEffects$type == 'eval','include'] <- TRUE
            ##if (attr(x$depvars[[i]],'uponly') ||attr(x$depvars[[i]],
            ##'downonly'))
            ##effects[['outdegree (density)']]$eval$fix <- TRUE
            ## maybe when you run it in siena07!
        }
        rateEffects$basicRate[1:observations] <- TRUE
        rateEffects$untrimmedValue <- rep(0, nrow(rateEffects))
        list(effects=rbind(rateEffects = rateEffects, objEffects = objEffects),
             starts=starts)
    }

    behaviorNet <- function(depvar, varname)
    {
        nodeSet <- attr(depvar,'nodeSet')
        objEffects <- paste('behavior', varname,
                            behaviorObjEffects[1:2, 1])
        objFunctions <- paste('beh.', varname,
                              behaviorObjEffects[1:2, 2])
        objEndowment <- behaviorObjEffects[1:2, 3]
        objShortNames <- behaviorObjEffects[1:2, 4]
        objParms <- rep(0, length(objEffects))
        if (observations==1)
        {
            rateEffects <- paste('rate', varname)
            rateFunctions <- "Amount of behavioral change"
            rateShortNames <- 'Rate'
            ratePeriods <- 1
            rateTypes <- NA
       }
        else
        {
            rateEffects <- paste('rate ',varname,' (period ',
                                 periodNos, ')', sep='')
            rateFunctions <- paste("Amount of behavioral change in period",
                                   periodNos, 'on', varname)
            rateShortNames <- rep('Rate', observations)
            ratePeriods <- 1:observations
            rateTypes <- rep(NA, observations)
       }
        objEffects <- createObjEffectList(objEffects, objFunctions,
                                          objEndowment, objShortNames,
                                          objParms, varname)
        rateEffects <- createRateEffectList(rateEffects, rateFunctions,
                                            rateShortNames, ratePeriods,
                                            rateTypes, varname)
        for (j in seq(along=xx$depvars))
        {
            if (types[j] == 'oneMode' &&
                attr(xx$depvars[[j]], 'nodeSet') == nodeSet)
            {
                tmp <- netBehEff(varname, names(xx$depvars)[j])
                objEffects<- rbind(objEffects, tmp$objEff)
                rateEffects<- rbind(rateEffects, tmp$rateEff)
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
                netObjEffect <- paste('behavior ', varname,
                                      ': infl. one-sided ? x ', varname,
                                      ' alter', sep='')
                netObjFunction <- paste('beh. ', varname,
                                      ': infl. interaction? x ', varname,
                                      ' alter', sep='')
                netShortName <- 'behInfl1sid'
                objEff <-  createObjEffectList(netObjEffect, netObjFunction,
                                              TRUE, netShortName, 0,
                                               varname,
                                               varname2=names(xx$depvars)[j])
                objEffects<- rbind(objEffects, objEff)
            }
        }
        objEff <-  rep(paste('behavior ', varname,
                                  ': unspecified interaction', sep=''), 4)
        objFun <-  rep(paste('behavior ', varname,
                                  ': unspecified interaction', sep=''), 4)
        objEnd <-  rep(TRUE, 4)
        objShortNames <- rep('behUnspInt', 4)
        objParms <- rep(0, 4)
        objEffects <- rbind(objEffects,
                            createObjEffectList(objEff, objFun,
                                          objEnd, objShortNames, objParms, varname))
        objEffects$untrimmedValue <- rep(0, nrow(objEffects))
        starts <- getBehaviorStartingVals(depvar)
        if (!(attr(depvar,'anyUpOnly') || attr(depvar, 'anyDownOnly')))
        {
            effectname <- paste('behavior', varname, 'linear shape')
            objEffects[objEffects$effectName == effectname &
                       objEffects$type == 'eval','include']  <- TRUE
            objEffects[objEffects$effectName == effectname &
                       objEffects$type=='eval','initialValue']  <-
                           starts$tendency
             objEffects[objEffects$effectName == effectname &
                       objEffects$type=='eval','untrimmedValue']  <-
                           starts$untrimmed
           effectname <- paste('behavior', varname, 'quadratic shape')
            objEffects[objEffects$effectName == effectname &
                       objEffects$type == 'eval','include']  <- TRUE
            ## no starting value yet for quadratic effect
        }
        if (observations == 1)
            effectname <- paste('rate', varname)
        else
            effectname <- paste('rate ', varname,' (period ',
                                1:noPeriods, ')', sep='')
        rateEffects[rateEffects$effectName %in%
                    effectname, 'include'] <- TRUE
        rateEffects[rateEffects$effectName %in% effectname,
                    'initialValue'] <-  starts$startRate
        rateEffects$basicRate[1:observations] <- TRUE
        rateEffects$untrimmedValue <- rep(0, nrow(rateEffects))
        tmp <- objEffects$effectName[objEffects$type == 'endow']
        tmp <- sub('behavior', 'dec. beh.', tmp)
        objEffects$effectName[objEffects$type == 'endow'] <- tmp
        list(effects = rbind(rateEffects = rateEffects,
             objEffects = objEffects), starts=starts)
    }

    dyadNetObjEff<- function(covarname, symmetric)
    {
        if (symmetric)
        {
            covObjEffects <- covarname
            covObjFunctions <- paste('Sum of ties x', covarname)
        }
        else
        {
            covObjEffects <- c(paste(covarname,c('','x reciprocity')))
            covObjFunctions <- paste(c('Sum of ties x',
                                       'Sum reciprocated ties x'), covarname)
        }
        covObjEffects <- c(covObjEffects, paste(dyadObjEffects[, 1],
                                                covarname))
        covObjFunctions <- c(covObjFunctions, paste(dyadObjEffects[, 1],
                                                    covarname))
        covObjShortNames <- c('X', 'Xrecip', dyadObjEffects[,2])
        objEff <-  createObjEffectList(covObjEffects, covObjFunctions,
                                       rep(TRUE, length(covObjEffects)),
                                       covObjShortNames,
                                       rep(0, length(covObjEffects)),
                                       varname, varname2=covarname)
        list(objEff=objEff)
    }
    covSymmNetEff<- function(covarname, poszvar, moreThan2)
    {
        covEffects <- paste(covarname,covarSymmetricObjEffects[, 1])
        covEffects <- c(covEffects, paste('same', covarname))
        covEffects <- c(covEffects, paste(covarname, 'ego x', covarname,
                                          'alter'))
        covEffects <- c(covEffects, paste(covarname, 'of indirect ties'))
        covShortNames <- c(covarSymmetricObjEffects[, 2],
                           'sameX', 'egoXAlt')
        if (!poszvar)
        {
            covEffects <- covEffects[1:2]
            covShortNames <- covShortNames[1:2]
        }
        if (!moreThan2)
        {
            covEffects <- covEffects[-2]
            covShortNames <- covShortNames[-2]
        }

        rateEffects<- paste('effect', covarname, 'on rate')
        rateShortNames <- 'RateX'
        ratePeriods <- NA
        rateTypes <- 'covariate'
        objEff <-  createObjEffectList(covEffects, covEffects,
                                       rep(TRUE, length(covEffects)),
                                       covShortNames, 0,
                                       varname, varname2=covarname)
        rateEff <-  createRateEffectList(rateEffects, rateEffects,
                                         rateShortNames, ratePeriods,
                                         rateTypes, varname, varname2=covarname)
        list(objEff=objEff, rateEff=rateEff)
    }
    covNonSymmNetEff<- function(covarname, poszvar, moreThan2)
    {
        covEffects<- paste(covarname, covarNonSymmetricObjEffects[, 1])
        covFunctions<- paste(covarNonSymmetricObjEffects[, 2], covarname)
        covEffects<- c(covEffects, paste('same',covarname, c('','x reciprocity')))
        covEffects<- c(covEffects, paste(covarname,'ego x',covarname,
                                         c('alter', 'alter x recipr.')))
        covEffects <- c(covEffects, paste('higher',covarname))
        covEffects<- c(covEffects, paste(covarname,'of indirect ties'))
        covFunctions <- c(covFunctions,
                          paste(c('Same values on',
                                  'Same values x reciprocity on'), covarname))
        covFunctions <- c(covFunctions, paste ('Sum', covarname, 'ego x',
                                               covarname,
                                               c('alter', 'rec.alter')))
         covFunctions <- c(covFunctions,
                           paste(c('ego > alter for', 'Sum'),
                                 covarname, c('', 'of indirect ties')))
        covShortNames <- c(covarNonSymmetricObjEffects[, 3],
                           "sameX", "sameXRecip", "egoXaltX",
                           "egoXaltXRecip", "higher", "IndTies")
        endow <- rep(TRUE, length(covEffects))
        endow[length(endow)] <- FALSE
        if (!poszvar)
        {
            covEffects <- covEffects[3]
            covFunctions <- covFunctions[3]
            covShortNames <- covShortNames[3]
            endow <- endow[3]
        }
        else if (!moreThan2)
        {
            covEffects <- covEffects[-2]
            covFunctions <- covFunctions[-2]
            covShortNames <- covShortNames[-2]
            endow <- endow[-2]
        }
        rateEffects <- paste('effect', covarname,'on rate')
        rateFunctions <- paste('Amount of change x', covarname)
        rateShortNames <- 'RateX'
        ratePeriods <- NA
        rateTypes <- 'covariate'
        objEff <-  createObjEffectList(covEffects, covFunctions,
                                       endow, covShortNames,
                                       rep(0, length(covEffects)),
                                       varname, varname2=covarname)
        rateEff <-  createRateEffectList(rateEffects, rateFunctions,
                                         rateShortNames, ratePeriods,
                                         rateTypes, varname, varname2=covarname)
        list(objEff=objEff, rateEff=rateEff)
    }
    covBehEff<- function(varname, covarname, nodeSet, same=FALSE,
        ## same indicates that varname and covarname are the same:
        ## just one rate effect required
                           type=c('', 'Var', 'Beh'))
    {
        if (!same)
        {
            covObjEffects<- paste('behavior ', varname,': ',
                                  covarBehObjEffects[1, 1],' ',
                                  covarname, sep='')
            covObjFunctions<- paste('beh. ', varname, ' ',
                                    covarBehObjEffects[1, 2],' ',
                                covarname, sep='')
            covShortNames <- paste(covarBehObjEffects[1, 3], type, sep='')
        }
        covRateEffects<- paste('effect', covarname,'on rate', varname)
        covRateFunctions <- paste('Amount of change on', varname, 'x',
                                  covarname)
        rateShortNames <- 'RateX'
         ratePeriods <- NA
        rateTypes <- 'covariate'
        covname3 <- ""
        if (!same)
        {
            for (j in seq(along=xx$depvars))
            {
                if (types[j] == 'oneMode' &&
                    attr(xx$depvars[[j]], 'nodeSet') == nodeSet)
                {
                    covObjEffects <- c(covObjEffects,
                                       paste('behavior ', varname, ': ',
                                             covarBehObjEffects[2, 1], ' ',
                                             covarname, sep=''))
                    covObjFunctions <- c(covObjFunctions,
                                         paste('beh. ', varname, ' ',
                                               covarBehObjEffects[2, 2], ' ',
                                               covarname, sep=''))
                    covShortNames <- c(covShortNames,
                                       paste(covarBehObjEffects[2, 3],
                                             type, sep=''))
                    covname3 <- c(covname3, names(xx$depvars)[j])
                }
            }
        }
        if (!same)
        {
            objEff <-  createObjEffectList(covObjEffects, covObjFunctions,
                                           rep(TRUE, length(covObjEffects)),
                                           covShortNames,
                                           rep(0, length(covObjEffects)),
                                           varname, varname2=covarname,
                                           varname3=covname3)
        }
        else
            objEff <- NULL
        rateEff <-  createRateEffectList(covRateEffects, covRateFunctions,
                                         rateShortNames, ratePeriods,
                                         rateTypes, varname, varname2=covarname)
        list(objEff=objEff, rateEff=rateEff)
    }
    netBehEff<- function(varname, netname)
    {
        netObjEffects <- paste('behavior', varname,
                               behaviorObjEffects[-c(1, 2), 1])
        netObjFunctions <- paste('beh.', varname,
                                 behaviorObjEffects[-c(1, 2), 2])
        netShortNames <- behaviorObjEffects[-c(1, 2), 4]
        netRateEffects <- paste(behaviorRateEffects[-1, 1], varname)
        netRateFunctions <- paste('Amount of change on', varname,
                                  behaviorRateEffects[-1, 2])
        netRateShortNames <- behaviorRateEffects[-1, 3]
        netRatePeriods <- NA
        netRateTypes <- NA
        objEff <-  createObjEffectList(netObjEffects, netObjFunctions,
                                       rep(TRUE, length(netObjEffects)),
                                       netShortNames,
                                       rep(0, length(netObjEffects)),
                                       varname, varname2=netname)
         rateEff <-  createRateEffectList(netRateEffects, netRateFunctions,
                                             netRateShortNames, netRatePeriods,
                                             netRateTypes, varname,
                                          varname2=netname)
       list(objEff=objEff, rateEff=rateEff)
    }
    createObjEffectList<- function(effectnames, functionnames, endowment,
                                   shortnames, parms, varname, varname2="",
                                   varname3=NULL)
    {
        nn <- length(effectnames)
        effectnames <- rep(effectnames, times=(1 + as.numeric(endowment)))
        functionnames <- rep(functionnames, times=(1 + as.numeric(endowment)))
        shortnames <- rep(shortnames, times=(1 + as.numeric(endowment)))
        parms <- rep(parms, times=(1 + as.numeric(endowment)))
        if (!is.null(varname3))
            varname3 <- rep(varname3, times=(1 + as.numeric(endowment)))
        type <-  unlist(lapply(endowment, function(x) if (x)
                               c('eval', 'endow') else 'eval'))
        nn <- length(effectnames)
        if (is.null(varname3))
            varname3 <- rep("", nn)
        tmp <-data.frame(name=rep(varname, nn),
                         effectName=effectnames,
                         functionName=functionnames,
                         shortName=shortnames,
                         interaction1=rep(varname2, nn),
                         interaction2=varname3,
                         type=type,
                         basicRate=rep(FALSE, nn),
                         include=rep(FALSE, nn),
                         randomEffects=rep(FALSE, nn),
                         fix=rep(FALSE, nn),
                         test=rep(FALSE, nn),
                         initialValue=rep(0, nn),
                         parm=parms,
                         functionType=rep('objective', nn),
                         period = rep(NA, nn),
                         rateType = rep(NA, nn),
                         stringsAsFactors=FALSE)
        effectFn <- vector('list', nn)
        statisticFn <- vector('list', nn)
        tmp$effectFn <- effectFn
        tmp$statisticFn <- statisticFn
        tmp
    }
    createRateEffectList<- function(effectnames, functionnames, shortnames,
                                    ratePeriods, rateTypes,
                                    varname, varname2="")
    {
       # cat(ratePeriods, '\n')
        nn <- length(effectnames)
        tmp <- data.frame(name=rep(varname, nn),
                          effectName=effectnames,
                          functionName=functionnames,
                          shortName=shortnames,
                          interaction1=rep(varname2, nn),
                          interaction2=rep("", nn),
                          type=rep('rate', nn),
                          basicRate=rep(FALSE, nn),
                          include=rep(FALSE, nn),
                          randomEffects=rep(FALSE, nn),
                          fix=rep(FALSE, nn),
                          test=rep(FALSE, nn),
                          initialValue=rep(0, nn),
                          parm=rep(0, nn),
                          functionType=rep('rate', nn),
                          period = ratePeriods,
                          rateType = rateTypes,
                         stringsAsFactors=FALSE)
        effectFn <- vector('list', nn)
        statisticFn <- vector('list', nn)
        tmp$effectFn <- effectFn
        tmp$statisticFn <- statisticFn
        tmp
    }
#### start of function createEffects
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
### validate the object?
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
    if (any(sparses))
        require(Matrix)
    nOneModes <- sum(types == 'oneMode')
    nBehaviors <- sum(types == 'behavior')
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
               bipartite = {},
               stop('error type'))
        effects[[i]]$groupName <- groupNames[1]
        effects[[i]]$group <- 1
    }
    ## add starting values for the other objects
    if (groupx && length(x) > 1)
    {
        period <- xx$observations   ### periods used so far
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
                       effects[[eff]][use, 'include'] <- TRUE
                       effects[[eff]][use, 'initialValue'] <-
                                          starts$startRate
                       effects[[eff]][use, 'groupName'] <- groupNames[group]
                       effects[[eff]][use, 'group' ]<- group
                       effects[[eff]][use, 'period'] <-  period - 1 +
                          1:noPeriods
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
                       effects[[eff]][use, 'include'] <- TRUE
                       effects[[eff]][use, 'initialValue'] <- tendency
                       effects[[eff]][use, "untrimmedValue"] <- untrimmed
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
                       effects[[eff]][use, 'include'] <- TRUE
                       effects[[eff]][use, 'initialValue'] <-
                                               starts$startRate
                       effects[[eff]][use, 'groupName'] <- groupNames[group]
                       effects[[eff]][use, 'group'] <- group
                       effects[[eff]][use, 'period'] <- 1:noPeriods
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
                           effects[[eff]][effects[[eff]]$effectName ==
                                          'degree (density)' &
                                          effects[[eff]]$type == 'eval',
                                          'initialValue'] <- degree
                            effects[[eff]][effects[[eff]]$effectName ==
                                          'degree (density)' &
                                          effects[[eff]]$type == 'eval',
                                          'untrimmedValue'] <- untrimmed
                      }
                       else
                       {
                           if (!(attr(x,'anyUpOnly') || attr(x, 'anyDownOnly')))
                           {
                               effects[[eff]][effects[[eff]]$effectName ==
                                              'outdegree (density)' &
                                              effects[[eff]]$type == 'eval',
                                              'initialValue'] <- degree
                                effects[[eff]][effects[[eff]]$effectName ==
                                              'outdegree (density)' &
                                              effects[[eff]]$type == 'eval',
                                              'untrimmedValue'] <- untrimmed
                          }
                       }
                       effects

                   },
                       bipartite = {},
                       stop('error type'))
            }
        }
    }
    effects <- do.call(rbind, effects)
    attr(effects, "starts") <- NULL
    cl <- class(effects)
    if (groupx)
        class(effects) <- c('groupEffects','effects', cl)
    else
        class(effects) <- c('effects', cl)
    effects
}

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
            mintab <- table(dif[z[, x] == rr[1], x] > 0)
            maxtab <- table(dif[z[, x] == rr[2], x] < 0)
            val <- (mintab[2] + 1) / (sum(mintab) + 2) +
                (maxtab[2] + 1) / (sum(maxtab) + 2)
            if (val > 0.9) val <- 0.5
            c(-log(1 - val), mintab, maxtab)
        }, z = depvar, y = dif)
        startRate <- tmp[1, ]
        ##tendency
        tmp <- rowSums(tmp[-1, ]) + 2
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
getNetworkStartingVals <- function(depvar, structValid=TRUE)
{
    noPeriods <- attr(depvar, "netdims")[3] - 1
    ##rate
    ##get distance and number of valid links
    if (!attr(depvar,'sparse'))
    {
        nactors <- nrow(depvar)
        if (structValid)
        {
            use <- !is.na(depvar) & (depvar == 10 | depvar == 11)
            depvar[use] <- depvar[use] - 10  ## remove structural values
        }
        else
        {
            depvar[depvar==10 | depvar==11] <- NA ## remove structural values
        }
        tmp <- sapply(1:noPeriods, function(x, z){
            diag(z[ , , x]) <- NA
            diag(z[, , x + 1]) <- NA
            matdiff <- sum(z[, , x + 1] != z[, , x], na.rm=TRUE)
            matchange <- table(z[, , x + 1], z[, , x])
            matcnt <- nactors * nactors -
                sum(is.na(z[, , x + 1]) | is.na(z[, , x]))
            c(matcnt=matcnt, matdiff=matdiff, matchange=matchange)
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
            if (structValid)
            {
                use <- mymat1@x %in% c(10, 11)
                mymat1@x[use] <- mymat1@x[use] - 10
                use <- mymat2@x %in% c(10, 11)
                mymat2@x[use] <- mymat2@x[use] - 10
            }
            else
            {
                mymat1@x[mymat1@x==10] <- NA
                mymat1@x[mymat1@x==11] <- NA
                mymat2@x[mymat2@x==10] <- NA
                mymat2@x[mymat2@x==11] <- NA
            }
            diag(mymat1) <- NA
            diag(mymat2) <- NA
            mydif <- mymat2-mymat1
            matdiff[i] <-sum(abs(mydif), na.rm=TRUE)
            tmp <- table(mydif@x)
            tmp00 <- nactors * nactors - length(mydif@x)
            tmp <- c(tmp00, tmp[c(3, 1, 2)])
            matchange[,i] <- tmp
            matcnt[i] <- sum(tmp)
        }
        matchange <- data.frame(matchange)
        tmp <-as.matrix(rbind(matcnt=matcnt, matdiff=matdiff,
                              matchange=matchange))
    }
    distance <- attr(depvar, "distance" )
    if (attr(depvar,'symmetric'))
        startRate<- nactors * (0.2 + distance)/(tmp['matcnt',]+1)
    else
        startRate<- nactors * (0.2 + 2 * distance)/(tmp['matcnt',]+1)
    startRate <- pmax(0.1, startRate)
    startRate <- pmin(100, startRate)
    ##degree
    matchange<- as.matrix(tmp[grep("matchange", rownames(tmp)),,drop=FALSE])
    if (attr(depvar,'symmetric'))
    {
        matchange <- matchange %/% 2
        matcnt <- matcnt %/% 2
    }
    p01 <- ifelse (matchange[1,] + matchange[2,] >=1,
                   matchange[2,]/(matchange[1,]+matchange[2,]),0.5)
    p10 <- ifelse (matchange[3,] + matchange[4,] >=1,
                   matchange[3,]/(matchange[3,]+matchange[4,]),0.5)
    p01 <- pmax(0.02,p01)
    p10 <- pmax(0.02,p10)
    p01 <- pmin(0.98,p01)
    p10 <- pmin(0.98,p10)
    alpha <- 0.5 * log(p01/p10)
    ##  if (observations == 2) ##more observations may come later!
    ##       alphaf1 <- alpha
    ##  else
    ## {
    p00 <- ifelse (matchange[1,] + matchange[2,] >=1,
                   matchange[1,]/(matchange[1,]+matchange[2,]),0.0)
    p11 <- ifelse (matchange[3,] + matchange[4,] >=1,
                   matchange[4,]/(matchange[3,]+matchange[4,]),0.0)
    p00 <- pmax(0.02,p00)
    p11 <- pmax(0.02,p11)
    p00 <- pmin(0.98,p00)
    p11 <- pmin(0.98,p11)
    prec <- ifelse(matchange[2,] * matchange[3,] >=1,
                   4 /((p00/matchange[2,]) +
                       (p11/matchange[3,])),1e-6)
    alphaf1 <- sum(alpha*prec/sum(prec))
    ## }
    untrimmed <- alphaf1
    alphaf1 <- ifelse(alphaf1 < -3, -3, ifelse(alphaf1 > 3, 3, alphaf1))
    list(startRate=startRate, degree=alphaf1, alpha=alpha, prec=prec, tmp=tmp,
        untrimmed = untrimmed)
}
