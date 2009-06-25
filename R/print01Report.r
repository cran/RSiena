print01Report <- function(data, myeff, modelname="Siena", session=NULL)
{
    reportDataObject1 <- function(x)
    {
        Report(c(x$observations, "observations,\n"), outf)
        ##
        if (length(x$nodeSets) > 1)
        {
            Report("Node Sets:\n", outf)
            lapply(x$nodeSets, function(z)
               {
                   Report(c(" ", format(attr(z, "nodeSetName"), width=15),
                            ":",
                            format(length(z), width=3), "nodes\n"), outf)
               })
            Report("\n", outf)
        }
        else
        {
            Report(c(length(x$nodeSets[[1]]), "actors\n"), outf)
        }
    }
    reportDataObject <- function(x, periodFromStart=0, multi=FALSE)
    {
        reportStart <- function()
        {
            multipleNodeSets <- length(x$nodeSets) > 1
            if (multipleNodeSets)
            {
                Report("Dependent variables   Type      NodeSet(s) (R, C)\n",
                       outf)
                Report("-------------------   ----      -----------------\n",
                       outf)
                for (i in 1:length(x$depvars))
                {
                    atts <- attributes(x$depvars[[i]])
                    Report(c(format(atts$name, width=20),
                             format(atts$type, width=10)), outf)
                    for (j in 1:length(atts$nodeSets))
                    {
                        Report(c(format(atts$nodeSets, width=10),
                                 "(", atts$dim[j], ")"), sep="", outf)
                    }
                    Report("\n", outf)
                }
            }
            else
            {
                Report(c(x$observations, "observations,\n"), outf)
                Report(c(length(x$nodeSets[[1]]), "actors,\n"), outf)
                Report(c(sum(types=="oneMode"),
                         "dependent network variables,\n"),
                       outf)
                Report(c(sum(types=="bipartite"),
                         "dependent bipartite variables,\n"), outf)
                Report(c(sum(types=="behavior"),
                         "dependent behavior variables,\n"),
                       outf)
            }
            Report(c(length(x$cCovars), "constant actor covariates,\n"), outf)
            Report(c(length(x$vCovars),
                     "exogenous changing actor covariates,\n"), outf)
            Report(c(length(x$dycCovars), "constant dyadic covariates,\n"),
                   outf)
            Report(c(length(x$dyvCovars),
                     "exogenous changing dyadic covariates,\n"), outf)
            Report(c(c('no files',
                       'file')[1 + as.numeric(length(x$compositionChange))],
                     "with times of composition change.\n"), outf)
            if ((length(x$cCovars) > 0 || length(x$dycCovars) > 0) && multi)
            {
                Report(c("For multi-group projects, constant covariates are",
                         "treated as changing covariates.\n"), outf)
                if (length(x$dycCovars) > 0)
                {
                    Report(c("Note that missings in changing dyadic",
                             "covariates are not (yet) supported!\n"), outf)
                }
            }
         Report("\n", outf)
        }

        reportNetworks <- function()
        {
            Heading(2, outf, "Reading network variables.")
            anymissings <- FALSE
            for (i in 1:length(x$depvars))
            {
                depvar <- x$depvars[[i]]
                atts <- attributes(depvar)
                netname <- atts$name
                type <- atts$type
                if (type != "behavior")
                {
                    Report("Name of ", outf)
                    if (nNetworks > 1)
                    {
                        Report("this ", outf)
                    }
                    Report(c("network variable: ", netname, '.\n'),
                           sep="", outf)
                    Report(c(type, "network.\n"), outf)
                    if (type == "bipartite")
                    {
                        Report("This is a two-mode network.\n", outf)
                        Report(c("The number of units in the second mode is ",
                                 atts$dim[2], ".\n"), sep="", outf)
                    }
                    for (k in 1:x$observations)
                    {
                        if (!is.null(session))
                        {
                            filename <-
                                session$Filename[session$Name == netname
                                                 &session$Period == k]
                            Report(c("Observation moment ", k + periodFromStart,
                                     " was read from file ", filename, '. \n'),
                                   sep='', outf)
                        }
                        Report(c("For observation moment ", k + periodFromStart,
                                 ", degree distributions are as ",
                                 "follows:\nNodes\n"),
                               sep="", outf)
                        ## remove structurals ? NA or 0/1
                        if (attr(depvar, "sparse"))
                        {
                            require(Matrix)
                            tmpdepvar <- depvar[[k]]
                            tmpx1 <- tmpdepvar@x
                            use <- tmpx1 %in% c(10, 11)
                            tmpx1[use] <- tmpx1[use] - 10
                            tmpdepvar@x <- tmpx1
                            outdeg <- rowSums(tmpdepvar, na.rm=TRUE)
                            indeg <- colSums(tmpdepvar, na.rm=TRUE)
                            diag(tmpdepvar) <- 0
                            missrow <- rowSums(is.na(depvar[[k]]))
                            misscol <- colSums(is.na(depvar[[k]]))
                        }
                        else
                        {
                            tmpdepvar <- depvar[, , k]
                            use <- tmpdepvar %in% c(10, 11)
                            tmpdepvar[use] <- tmpdepvar[use] - 10
                            outdeg <- rowSums(tmpdepvar, na.rm=TRUE)
                            indeg <- colSums(tmpdepvar, na.rm=TRUE)
                            diag(tmpdepvar) <- 0
                            missrow <- rowSums(is.na(tmpdepvar))
                            misscol <- colSums(is.na(tmpdepvar))
                        }
                        tmp <- format(cbind(1:atts$netdims[1], outdeg, indeg))
                        Report(tmp[, 1], fill=60, outf)
                        Report("out-degrees\n", outf)
                        Report(tmp[, 2], fill=60, outf)
                        Report("in-degrees\n", outf)
                        Report(tmp[, 3], fill=60, outf)
                        ## report structural values
                        if (attr(depvar, "structural"))
                        {
                            if (attr(depvar, "sparse"))
                            {
                                nstruct0 <- sum(depvar[[k]] %in% c(10))
                                nstruct1 <- sum(depvar[[k]] %in% c(11))
                            }
                            else
                            {
                                nstruct0 <- sum(depvar[, , k] %in% c(10))
                                nstruct1 <- sum(depvar[, , k] %in% c(11))
                            }
                            if (nstruct0 + nstruct1 > 0)
                            {
                                Report(c("\nThe input file contains codes for ",
                                         "structurally determined values:\n"),
                                       sep="", outf );
                                if (attr(depvar, "sparse"))
                                {
                                    nstruct0 <- sum(depvar[[k]] %in% c(10))
                                    nstruct1 <- sum(depvar[[k]] %in% c(11))
                                }
                                else
                                {
                                    nstruct0 <- sum(depvar[, , k] %in% c(10))
                                    nstruct1 <- sum(depvar[, , k] %in% c(11))
                                }
                                Report(c('  ', nstruct0, ' structural zero'),
                                       sep='', outf)
                                Report(ifelse(nstruct0 > 1,
                                              "s were found (code 10).\n",
                                              " was found (code 10).\n"), outf)
                                Report(c('  ', nstruct1, ' structural one'),
                                       sep='', outf)
                                Report(ifelse(nstruct1 > 1,
                                              "s were found (code 11).\n",
                                              " was found (code 11).\n"),
                                       outf)
                                ##
                                if (attr(depvar, 'sparse'))
                                {
                                    nnonactive <-
                                        rowSums(depvar[[k]] == 10 |
                                                depvar[[k]] == 11, na.rm=TRUE)
                                    nnonactive <- nnonactive >= nrow(depvar[[k]])
                                }
                                else
                                {
                                    nnonactive <-
                                        rowSums(depvar[, , k] == 10 |
                                                depvar[, , k] == 11, na.rm=TRUE)
                                    nnonactive <- nnonactive >= nrow(depvar[, , k])
                                }
                                if (sum(nnonactive)  == 1)
                                {
                                    Report(c("Actor ", which(nnonactive),
                                             " is inactive at this ",
                                             "observation.\n"), sep='', outf)
                                }
                                else if (sum(nnonactive) > 1)
                                {
                                    Report(c("Actors ", which(nnonactive),
                                             " are inactive at this ",
                                             "observation.\n"), sep='', outf)
                                }
                            }
                        }
                        if (attr(depvar, "sparse"))
                        {
                            anymissings <- any(is.na(depvar[[k]]))
                        }
                        else
                        {
                            anymissings <- any(is.na(depvar[, , k]))
                        }
                        if (anymissings)
                        {
                            Report(c("\nFor observation moment ",
                                     k + periodFromStart,
                                     ", number of missing values ",
                                     "are:\nNodes\n"),
                                   sep="", outf)
                            tmp <- format(cbind(1:atts$netdims[1],
                                                missrow, misscol))
                            Report(tmp[, 1], fill=60, outf)
                            Report("missing in rows\n", outf)
                            Report(tmp[, 2], fill=60, outf)
                            Report("missing in columns\n", outf)
                            Report(tmp[, 3], fill=60, outf)
                            Report(c("Total number of missing data: ",
                                     sum(missrow),
                                     ", corresponding to a fraction of ",
                                     round(sum(missrow)/atts$netdims[1] /
                                           (atts$netdims[1] - 1), 3),
                                     ".\n"), sep="", outf)
                            if (k > 1)
                                Report(c("In reported in- and outdegrees,",
                                         "missings are not counted.\n"), outf)
                            Report("\n", outf)
                        }
                        else
                        {
                            Report(c("\nNo missing data for observation ",
                                     k + periodFromStart, ".\n\n"),
                                   sep= "", outf)
                        }
                    }
                    if (anymissings)
                    {
                        Report(c("There are missing data for this",
                               "network variable,\n"), outf)
                        Report(c("and the <<carry missings forward>>",
                               "option is active.\n"), outf)
                        Report("This means that for each tie variable,\n", outf)
                        Report(c("the last previous nonmissing value (if any)",
                                 "is imputed.\n"), outf)
                        Report(c("If there is no previous nonmissing value,",
                                 "the value 0 is imputed.\n"), outf)
                    }
                }
            }
            Report("\n", outf)
        }
        reportBehaviors <- function()
        {
            Heading(2, outf, "Reading dependent actor variables.")
            anymissings <- FALSE
            iBehav <- 0
            for (i in 1:length(x$depvars))
            {
                if (types[i] == "behavior")
                {
                    depvar <- x$depvars[[i]]
                    atts <- attributes(depvar)
                    netname <- atts$name
                    type <- atts$type

                    iBehav <- iBehav + 1
                    mymat <- depvar[, 1, ]
                    mystr <- paste(iBehav, switch(as.character(iBehav),
                                                  "1"=, "21"=, "31"= "st",
                                                  "2"=, "22"=, "32"= "nd",
                                                  "3"=, "23"=, "33"= "rd",
                                                  "th"), sep="")
                    Report(c(mystr, " dependent actor variable named ",
                             netname), sep="", outf)
                    if (!is.null(session))
                    {
                        filename <-
                            session$Filename[session$Name == netname]

                            Report(c(" was read from file ", filename, ".\n"),
                                   sep="", outf)
                    }
                    else
                    {
                        Report(".\n", outf)
                    }
                    ranged <- atts$range2
                    Report(c("Maximum and minimum rounded values are ",
                             round(ranged[1]), " and ",
                             round(ranged[2]), ".\n"), sep="", outf)
                    if (ranged[1] < 0 )
                        stop("Negative minima not allowed for dependent actor ",
                             "variables.\n")
                    if (ranged[2] > 255 )
                        stop("Maxima more than 255 not allowed for dependent",
                             "actor ", "variables.\n")
                    if (ranged[1] >= ranged[2] )
                        stop("Dependent actor variables must not be",
                             " constant.\n")
                    if (any(is.na(depvar)))
                    {
                        Report(c("Missing values in this actor variable are",
                                 "imputed",
                                 "by the mode per observation.\n"), outf)
                        Report(c("But if there is a previous nonmissing",
                                 "value,",
                                 "this is used as the imputed value.\n"), outf)
                        Report("Modal values:\nObservation  ", outf)
                        Report(c(format(1:x$observations+periodFromStart,
                                        width=4), '\n'), outf)
                        Report(c(format("Modes", width=12),
                                 format(atts$modes, width=4)), outf)
                        Report("\n", outf)
                    }
                    Report('\n', outf)
                }
            }
            Report(c("\nA total of",
                     nBehavs, "dependent actor variable"), outf)
            Report(ifelse(nBehavs > 1, "s.\n\n", ".\n\n"), outf)
            Report("Number of missing cases per observation:\n", outf)
            Report(c(" observation", format(1:x$observations+periodFromStart,
                                            width=10),
                     "      overall\n"), sep="", outf)
            for (i in 1:length(x$depvars))
            {
                if (types[i] == "behavior")
                {
                    depvar <- x$depvars[[i]][, 1, ]
                    netname <- atts$name
                    missings <- colSums(is.na(depvar))
                    Report(c(format(netname, width=12),
                             format(c(missings, sum(missings)),
                                    width=10), "      (",
                             format(round(sum(missings)/
                                          nrow(depvar)/ncol(depvar), 1),
                                    nsmall=1, width=4), ' %)\n'), sep="", outf)
                }
            }
            Report("\nMeans per observation:\n", outf)
            Report(c(" observation", format(1:x$observations+periodFromStart,
                                            width=10),
                     "      overall\n"), sep="", outf)
            for (i in 1:length(x$depvars))
            {
                if (types[i] == "behavior")
                {
                    depvar <- x$depvars[[i]][, 1, ]
                    netname <- atts$name
                    means <- colMeans(depvar, na.rm=TRUE)
                    Report(c(format(netname, width=14),
                             format(round(means, 3), nsmall=3,
                                    width=10), format(round(mean(means),
                                    3), width=10), '\n\n'), sep="", outf)
                }
            }
        }
        reportConstantCovariates <- function()
        {
            nCovars <- length(x$cCovars)
            covars <- names(x$cCovars)
            Heading(2, outf, "Reading constant actor covariates.")
            Report(c(nCovars, "variable"),outf)
            Report(ifelse(nCovars == 1, ", named:\n", "s, named:\n"), outf)
            for (i in seq(along=covars))
            {
                Report(c(format(covars[i], width=15), '\n'), outf)
            }
            Report(c("\nA total of", nCovars,
                     "non-changing individual covariate"), outf)
            Report(ifelse(nCovars == 1, ".\n\n", "s.\n\n"), outf)
            Report("Number of missing cases:\n", outf)
            for (i in seq(along=covars))
            {
                Report(c(format(covars[i], width=15),
                         sum(is.na(x$cCovars[[i]])), "  (",
                         format(round(sum(is.na(x$cCovars[[i]]))/
                                      length(x$cCovars[[i]]), 1),
                                width=3, nsmall=1), '%)\n'), outf)
            }
            Report("\nInformation about covariates:\n", outf)
            Report(c(format("minimum  maximum     mean", width=38,
                            justify="right"), "\n"), outf)
            for (i in seq(along=covars))
            {
                atts <- attributes(x$cCovars[[i]])
                Report(c(format(covars[i], width=10),
                         format(round(atts$range2[1], 1),
                                nsmall=1, width=8),
                         format(round(atts$range2[2], 1),
                                nsmall=1, width=7),
                         format(round(atts$mean, 3),
                                nsmall=3, width=10), "\n"), outf)
            }
            Report(c("The mean value", ifelse(nCovars == 1, " is", "s are"),
                     " subtracted from the covariate",
                     ifelse(nCovars == 1, ".\n\n", "s.\n\n")), sep="", outf)
        }
        reportChangingCovariates <- function()
        {
            nCovars <- length(x$vCovars)
            covars <- names(x$vCovars)
            use <- ! covars %in% names(x$cCovars)
            nCovars <- length(x$vCovars[use])
            Heading(2, outf, "Reading exogenous changing actor covariates.")
            Report(c(nCovars, "variable"),outf)
            Report(ifelse(nCovars == 1, ", named:\n", "s, named:\n"), outf)
            for (i in seq(along=covars[use]))
            {
                Report(c(format(covars[use][i], width=15), '\n'), outf)
            }
            Report(c("\nA total of", nCovars,
                     "exogenous changing actor covariate"), outf)
                Report(ifelse(nCovars == 1, ".\n\n", "s.\n\n"), outf)
            Report("Number of missing cases per period:\n", outf)
            Report(c(" period   ", format(1:(x$observations - 1) +
                                         periodFromStart, width=9),
                     "       overall\n"), sep="", outf)
            for (i in seq(along=covars))
            {
                if (use[i])
                {
                    thiscovar <- x$vCovars[[i]] ## matrix
                    misscols <- colSums(is.na(thiscovar))
                    Report(c(format(covars[i], width=10),
                             format(misscols, width=8),
                             format(sum(misscols), width=9), "     (",
                             format(round(sum(misscols)/nrow(thiscovar)/
                                          ncol(thiscovar), 1), nsmall=1,
                                    width=3), '%)\n'), outf)
                }
            }
            Report("\nInformation about changing covariates:\n", outf)
            Report(c(format("minimum  maximum     mean", width=38,
                            justify="right"), "\n"), outf)
            for (i in seq(along=covars))
            {
                if (use[i])
                {
                    atts <- attributes(x$vCovars[[i]])
                    Report(c(covars[i], '\n'), outf) # name
                    for (j in 1:(ncol(x$vCovars[[i]])))
                    {
                        Report(c("  period", format(j + periodFromStart,
                                                   width=3),
                                 format(round(atts$rangep[1, j], 1),
                                        nsmall=1, width=7),
                                 format(round(atts$rangep[2, j], 1),
                                        nsmall=1, width=7),
                                 format(round(atts$meanp[j], 3),
                                        nsmall=3, width=10), "\n"), outf)
                    }
                    Report(c(format("Overall", width=28),
                             format(round(atts$mean, 3), width=10, nsmall=3),
                             "\n"), outf)

                }
            }
            Report("\nThe overall mean value", outf)
            Report(c(ifelse(nCovars  == 1, " is", "s are"),
                     "subtracted from the covariate.\n\n"),  outf)
        }
        reportConstantDyadicCovariates <- function()
        {
            nCovars <- length(x$dycCovars)
            covars <- names(x$dycCovars)
            Heading(2, outf, "Reading constant dyadic covariates.")
            for (i in seq(along=covars))
            {
                Report(c("Dyadic covariate named ", covars[i], '.\n'),
                       sep="", outf)
            }
            Report(c("\nA total of", nCovars,
                     "dyadic individual covariate"), outf)
            Report(ifelse(nCovars == 1, ".\n\n", "s.\n\n"), outf)
            Report("Number of tie variables with missing data:\n", outf)
            for (i in seq(along=covars))
            {
                myvar <- x$dycCovars[[i]]
                diag(myvar) <- 0
                Report(c(format(covars[i], width=15),
                         sum(is.na(myvar)), "  (",
                         format(round(sum(is.na(myvar))/
                                      (length(myvar) - nrow(myvar)), 1),
                                width=3, nsmall=1), '%)\n'), outf)
            }
            Report("\nMeans of  covariates:\n", outf)
            Report(c(format("minimum  maximum     mean", width=38,
                            justify="right"), "\n"), outf)
            for (i in seq(along=covars))
            {
                atts <- attributes(x$dycCovars[[i]])
                Report(c(format(covars[i], width=10),
                         format(round(atts$range2[1], 1),
                                nsmall=1, width=8),
                         format(round(atts$range2[2], 1),
                                nsmall=1, width=7),
                         format(round(atts$mean, 3),
                                nsmall=3, width=10), "\n"), outf)
            }
            Report('\n', outf)
            Report(c(ifelse(nCovars == 1,
                            "This mean value is ", "These mean values are "),
                     "subtracted from the dyadic covariate",
                     ifelse(nCovars == 1, ".\n\n", "s.\n\n")), sep="", outf)
        }
        reportChangingDyadicCovariates <- function()
        {
            covars <- names(x$dyvCovars)
            use <- ! covars %in% names(x$dycCovars) ## need an attributes to say
            nCovars <- length(x$dyvCovars[use])
            Heading(2, outf, "Reading exogenous dyadic covariates.")
            Report(c("Note that no missing values are considered yet for",
                     "changing dyadic covariates.\n"), outf)
            for (i in seq(along=covars))
            {
                Report(c("Exogenous dyadic covariate named ", covars[i], '.\n'),
                       sep="", outf)
            }
            Report("\nMeans of  covariates:\n", outf)
            for (i in seq(along=covars))
            {
                atts <- attributes(x$dyvCovars[[i]])
                Report(c(format(covars[i], width=10),
                         format(round(atts$mean, 3),
                                nsmall=3, width=10), "\n"), outf)
            }
            Report('\n', outf)
            Report(c(ifelse(nCovars == 1, "This ", "These"),
                     "global mean value",
                     ifelse(nCovars == 1, " is ", "s are "),
                     "subtracted from the changing dyadic covariate",
                     ifelse(nCovars ==1, ".\n\n", "s.\n\n")), sep="", outf)
        }
        reportCompositionChange <- function()
        {
            comps <- x$compositionChange
            nComps <- length(comps)
            Heading(2, outf, "Reading files with times of composition change.")
            for (i in seq(along=comps))
            {
                nodeSet <- attr(comps[[i]], "nodeSet")
                Report(c("Composition changes for nodeSet ", nodeSet, '.\n'),
                       sep="", outf)
                events <- attr(comps[[i]], "events")
                for (j in 1:nrow(events))
                {
                    x <- events[j, ]
                    Report(c("Actor ", format(x$actor, width=2),
                             ifelse(x$event=="join", " joins ", " leaves"),
                             " network at time ",
                             format(round(x$period + x$time, 4), nsmall=4),
                             ".\n"), sep="", outf)
                }
                pertab <- table(events$period, events$event)
                for (period in row.names(pertab))
                {
                    joiners <- pertab[period, "join"]
                    leavers <- pertab[period, "leave"]
                    Report(c("\nIn period ", period, ", ", joiners,
                             ifelse(joiners == 1, " actor", " actors"),
                             " joined and ", leavers,
                             ifelse(leavers == 1, " actor", " actors"),
                             " left the network.\n"), sep="", outf)
                }
            }
        }
        types <- lapply(x$depvars, function(z) attr(z, "type"))
        reportStart()
        nNetworks <- sum(types != "behavior")
        nBehavs <- sum(types == "behavior")
        if (nNetworks > 0)
        {
            reportNetworks()
        }
        if (nBehavs > 0)
        {
            reportBehaviors()
        }
        if (length(x$cCovars) > 0)
        {
            reportConstantCovariates()
        }
        if (nData > 1 && length(x$vCovars) > length(x$cCovars) ||
            (nData ==1  && length(x$vCovars) > 0))
        {
            reportChangingCovariates()
        }
        if (length(x$dycCovars) > 0)
        {
            reportConstantDyadicCovariates()
        }
        if (nData > 1 && length(x$dyvCovars) > length(x$dycCovars) ||
            (nData ==1  && length(x$dyvCovars) > 0))
        {
            reportChangingDyadicCovariates()
        }
        if (length(x$compositionChange) > 0)
        {
            reportCompositionChange()
        }
        Report("\n\n", outf) ## end of reportDataObject
    }
    ## create output file. ## start of print01Report proper
    Report(open=TRUE, type="w", projname=modelname)
    Report("                            ************************\n", outf)
    Report(c("                                   ", modelname, ".out\n"),
           sep='', outf)
    Report("                            ************************\n\n", outf)
    Report(c("Filename is ", modelname, ".out.\n\n"), sep="", outf)
    Report(c("This file contains primary output for SIENA project <<",
        modelname, ">>.\n\n"), sep="", outf)
    Report(c("Date and time:", format(Sys.time(), "%d/%m/%Y %X"), "\n\n"), outf)
    packageValues <- packageDescription("RSiena", fields=c("Version", "Date"))
    Report(c("SIENA version ", packageValues[[1]], " (",
        format(as.Date(packageValues[[2]]), "%d %m %Y"), ")\n\n"), sep="", outf)

    if (!inherits(data, 'sienaGroup'))
    {
        nData <- 1
        data <- sienaGroupCreate(list(data), singleOK=TRUE)
    }
    else
    {
        nData <- length(data)
    }
    if (nData > 1)
    {
        Report("Multi-group input detected\n\n", outf)
        for (i in 1:nData)
        {
            Report(c("Subproject ", i, ": <", names(data)[i], ">\n"), sep="",
                   outf)
            reportDataObject1(data[[i]])
        }
        Report(c("Multi-group project", modelname, "contains", nData,
                 "subprojects.\n\n"), outf)
        periodFromStart <- 0
        for (i in 1:nData)
        {
            Heading(1, outf,
                    paste("Subproject ", i, ": <", names(data)[i], ">",
                          sep="", collapse="")
                    )
            reportDataObject(data[[i]], periodFromStart, multi=TRUE)
            periodFromStart <- periodFromStart + data[[i]]$observations
       }
    }
    else
    {
        Heading(1, outf, "Data input.")
        reportDataObject(data[[1]], 0, multi=FALSE)
    }
    atts <- attributes(data)
    nets <- atts$types != "behavior"
    if (length(data) > 1)
    {
        Report("Series of observations for the multi-group project:\n", outf)
        periodFromStart <- 0
        for (i in seq(along=data))
        {
            Report(c(format(1:data[[i]]$observations + periodFromStart), '\n'),
                   outf)
            periodFromStart <- periodFromStart + data[[i]]$observations
        }
        Report("\n", outf)
   }
    if (sum(nets) > 0)
    {
        if (any(atts$anyMissing[nets]))
        {
            netnames <- atts$netnames[nets]
            missings <- atts$anyMissing[nets]
            for (i in seq(along=netnames[missings]))
            {
                Report(c("There are missing data for network variable ",
                         netnames[i], ".\n"), sep = "", outf)
            }
        }
        Report("\n", outf)
    }
  if (sum(atts$types == 'oneMode') > 0)
    {
        balmean <- atts$"balmean"

        Report(c("The mean structural dissimilarity value subtracted",
                 "in the\n"), outf)
        Report("balance calculations is ", outf)
        for (i in seq(along=atts$types))
        {
            if (atts$types[i] == "oneMode")
            {
                if (sum(atts$types == "oneMode") > 1)
                {
                    Report(c("Network name:", netnames[i],
                             format(round(balmean[i], 4), nsmall=4, width=14),
                             '.\n'),
                           sep="", outf)
                }
                else
                {
                    Report(c(format(round(balmean[i], 4), nsmall=4, width=14),
                             '.\n'),
                           sep="", outf)
                }
            }
        }
    }
    if (sum(atts$types == "behavior") > 0 ||
        (nData ==1 && length(atts$cCovars) > 0) ||
        length(atts$vCovars) > 0)
    {
        Report(c("\nFor the similarity variable calculated from each actor",
                 "covariate,\nthe mean is subtracted.\nThese means are:\n"),
               outf)
        if (nData == 1)
        {
            for (i in seq(along=atts$cCovars))
            {
                if (atts$cCovarPoszvar[i])
                {
                    Report(c("Similarity", format(atts$cCovars[i], width=12),
                             ':', format(round(atts$cCovarSim[i], 4), width=12,
                                         nsmall=4),
                             '\n'), outf)
                }
            }
        }
        for (i in seq(along=atts$netnames))
        {
            if (atts$types[i] == "behavior" && atts$bPoszvar[i])
            {
                Report(c("Similarity", format(atts$netnames[i], width=12),
                          ':', format(round(atts$bSim[i], 4), nsmall=4,
                                      width=12),
                         '\n'), outf)
            }
        }
        for (i in seq(along=atts$vCovars))
        {
            if (atts$vCovarPoszvar[i])
            {
                Report(c("Similarity", format(atts$vCovars[i], width=12),
                         ':', format(round(atts$vCovarSim[i], 4), width=12,
                                     nsmall=4),
                         '\n'), outf)
            }
        }
    }
    printInitialDescription(data, myeff, modelname)
    ##close the files
    Report(close=TRUE)
}













