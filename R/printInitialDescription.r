printInitialDescription <- function(data, effects, modelName="Siena")
{
    initialNetworks <- function()
    {
        Heading(2, outf,'Change in networks:')
        difficult <- FALSE
        for (net in which(types=="oneMode"))
        {
            if (nOneMode > 1)
            {
                Heading(3, outf, c("Network", netnames[net]))
            }

            Report(c("For the following statistics, missing values (if any)",
                     "are not counted.\n"), outf)
            if (any(attr(data, "structural")))
            {
                Report(c("Further, structurally determined entries are",
                         "treated as observed entries,\n"), outf)
                Report(c("except for calculating the distances, where missings",
                         "and structurally determined entries\n"), outf)
                Report(c("at the beginning or end of a period",
                     "are disregarded.\n"), outf)
                Report(c("This may imply that the columns 0 => 1 and 1 => 0",
                         "do not add up to the distance column.\n"), outf)
            }
            ## is this network symmetric?
            if (gpatts$symmetric[net])
            {
                Report(c(ifelse(nobs == 2, "Both", c("All ", nobs)),
                        " observed networks are symmetric.\n"), outf)
                Report(c("Therefore, it is assumed that this is an",
                         "analysis of an non-directed relation.\n"), outf)
            }
            Report("\nNetwork density indicators:\n", outf)
            density <- rep(NA, nobs)
            degree <- rep(NA, nobs)
            nties <- rep(NA, nobs)
            missings <- rep(NA, nobs)
            periodFromStart <- 0
            for (group in 1:nData)
            {
                j <- match(netnames[net], names(data[[group]]$depvars))
                if (is.na(j))
                    stop("network names not consistent")
                depvar <- data[[group]]$depvars[[j]]
                atts <- attributes(depvar)
                subs <- 1:data[[group]]$observations + periodFromStart
                ones <- sapply(atts$vals, function(x){
                    if (is.na(x["11"]))
                    {
                        x["1"]
                    }
                    else
                    {
                        x["1"] + x["11"]
                    }
                })
                density[subs] <- ones / atts$nval
                if (any(ones >= atts$nval))
                {
                    difficult <- TRUE
                }
                degree[subs] <- (atts$netdims[1] - 1) * ones / atts$nval
                nties[subs] <- ones
                if (gpatts$symmetric[net])
                {
                    nties <- nties / 2
                }
                missings[subs] <- 1 - atts$nval/ atts$netdims[1] /
                    (atts$netdims[1] - 1)
                periodFromStart <- data[[group]]$observations
            }
            ## now do the format
            tmp <- rbind(format(round(density, 3), nsmall=3, width=7),
                         format(round(degree, 3), nsmall=3, width=7),
                         format(nties, width=7),
                         format(round(missings, 3), nsmall=3, width=7))
            startCol <- 1
            text <- c("observation time", "density", "average degree",
                      "number of ties", "missing fraction")
            text <- c(format(text[1], width=24), format(text[-1], width=25))
            eols <- rep("\n", length(text))
         #   browser()
            repeat
            {
                endCol <- startCol + 6
                endCol <- min(endCol, nobs)
                tmp3 <- rbind(format(startCol:endCol, width=7),
                              tmp[, startCol:endCol])
                tmp2 <- cbind(text, tmp3, eols)
                Report(t(tmp2), sep="", outf)
                startCol <- startCol + 7
                if (startCol > nobs)
                    break
            }
            Report("\n\n", outf)
            Report(c(ifelse(gpatts$symmetric[net], "Edge", "Tie"),
                   "changes between subsequent observations:\n"), outf)
            valmin <- gpatts$netRanges[1, net]
            valmax <- gpatts$netRanges[2, net]
            tmp <- expand.grid(valmin:valmax, valmin:valmax)
            heads <- paste(tmp[, 2], "=> ", tmp[,1])
           # browser()
            Report(c(format(" periods", width=16),
                     format(heads, width=10),
                     ifelse(valmin == 0 && valmax == 1,
                            "Distance Jaccard   Missing\n",
                            "Distance     Missing\n")), sep = "", outf)
            periodFromStart <- 0

            for (group in 1:nData)
            {
                j <- match(netnames[net], names(data[[group]]$depvars))
                if (is.na(j))
                    stop("network names not consistent")
                depvar <- data[[group]]$depvars[[j]]
                tmp <- getNetworkStartingVals(depvar, structValid=TRUE)
                atts <- attributes(depvar)
                matchange <- tmp$tmp[-c(1:2), , drop=FALSE]
                ## if symmetric: divide by 2
                if (gpatts$symmetric[net])
                {
                    tmp$tmp <- tmp$tmp %/% 2
                    matchange <- matchange %/% 2
                }
                for (per in 1:(atts$netdims[3] - 1))
                {
                    ntot <- tmp$tmp["matcnt", per]
                    if (gpatts$symmetric[net])
                    {
                        misd <- atts$netdims[1] * (atts$netdims[1] - 1) / 2 -
                            ntot
                    }
                    else
                    {
                        misd <- atts$netdims[1] * (atts$netdims[1] - 1) - ntot
                    }
                    if (valmin == 0 && valmax == 1)
                    {
                        jaccard <- format(round(matchange[4, per] /
                                         (matchange[4, per] +
                                          matchange[3, per] +
                                          matchange[2, per]), 3), nsmall=3,
                                          width=10)

                        Report(c(format(per + periodFromStart, width=3),
                                 " ==> ",
                                 format(per + 1 + periodFromStart, width=2),
                                 format(matchange[, per], width=10),
                                 format(attr(depvar, "distance")[per],
                                        width=10),
                                 jaccard, format(misd, width=6), " (",
                                 round(100 * misd/(ntot + misd)), "%)\n"),
                               sep="", outf)
                     }
                    else
                    {
                          Report(c(per + periodFromStart, " ==> ",
                                 format(per + 1 + periodFromStart, width=2),
                                 format(matchange[, per], width=10),
                                 format(attr(depvar,"distance")[per], width=10),
                                 format(misd, width=7), " (",
                                   round(100 * misd/(ntot + misd)), "%)\n"),
                                 sep="", outf)

                      }
                    if (valmin == 0 && valmax == 1)
                    {
                        if (matchange[4, per] / (matchange[3, per] +
                                                 matchange[4, per]) <
                            matchange[2, per] / (matchange[2, per] +
                                                 matchange[1, per]))
                        {
                            Report(c("\nThis means that in period ", per,
                                     ", proportionately less 1-ties stayed 1,\n",
                                     " than 0-ties became 1. A great reversal",
                                     " of the network pattern!\n",
                                     "For some model specifications this may",
                                     " lead to problems in estimation.\n"), outf)
                        }
                    }
                }
                periodFromStart <- periodFromStart + atts$netdims[3]
                Report("\n", outf)
            }
            if (gpatts$symmetric[net])
            {
                Report(c("The distances reported in the output file",
                         "for conditional estimation\n",
                         "for the network variable refer to the total",
                         "symmetric adjacency matrix,\n",
                         "and therefore are double the distance",
                         "reported above.)\n\n"), outf)
            }
            Report("Dyad Counts:\n", outf)
            if (valmin == 0 & valmax == 1)
            {
                Report(" observation    total    mutual    asymm.     null\n",
                       outf)
            }
            else
            {
                Report(" observation    total \n", outf)
                heads <- expand.grid(valmin:valmax, valmin:valmax)
                heads <- heads[, c(2, 1)] ## reverse the column order
                heads <- heads[heads[,2] >= heads[,1], ] ## remove col2 < col1
                heads <- paste("(", heads[, 1], ",", heads[, 2], ")", sep="")
                Report(heads, outf)
                Report("\n", outf)
                ## this one needs more work!
            }
            periodFromStart <- 0
            for (group in 1:nData)
            {
                j <- match(netnames[net], names(data[[group]]$depvars))
                if (is.na(j))
                    stop("network names not consistent")
                depvar <- data[[group]]$depvars[[j]]
                atts <- attributes(depvar)

                for (per in 1:(atts$netdims[3]))
                {
                    ## find the dyad tables
                    if (atts$sparse)
                    {
                      #  require(Matrix)
                        mymat <- depvar[[per]]
                        mymat1 <- mymat@i
                        mymat2 <- mymat@j
                        mymat3 <- mymat@x
                        if (any(duplicated(paste(mymat1, mymat2))))
                            stop("sparse matrix has duplicate triples")
                        mymat3[mymat3==10] <- 0
                        mymat3[mymat3==11] <- 1
                        mymat1 <- mymat1[is.na(mymat3) | mymat3 != 0]
                        mymat2 <- mymat2[is.na(mymat3) | mymat3 != 0]
                        mymat3 <- mymat3[is.na(mymat3) | mymat3 != 0]
                        ## for not 0/1 need to tabulate the equivs to mutuals
                        ## and asymms by
                        ## the non zero one.
                        if (valmin == 0 && valmax ==1)
                        {
                            ij <- paste(mymat1[!is.na(mymat3)],
                                        mymat2[!is.na(mymat3)])
                            ji <- paste(mymat2[!is.na(mymat3)],
                                        mymat1[!is.na(mymat3)])
                            missij <- paste(mymat1[is.na(mymat3)],
                                            mymat2[is.na(mymat3)])
                            missji <- paste(mymat2[is.na(mymat3)],
                                            mymat1[is.na(mymat3)])
                            mutual <- sum(ij %in% ji) / 2
                            nondyads <- sum(ji %in% missij)
                            asymm <- length(ij) - nondyads - mutual * 2
                            missdyads <- sum(!missij %in% missji) +
                                sum(missij %in% missji) / 2
                            nulls <- atts$netdims[1] *
                                (atts$netdims[2] - 1) / 2 -
                                    missdyads - mutual - asymm
                            totDyad <- nulls + mutual + asymm
                        }
                    }
                    else
                    {
                        mymat <- depvar[, , per]
                        mymat[mymat == 10] <- 0
                        mymat[mymat == 11] <- 1
                        diag(mymat) <- NA
                        dyadTable <- table(mymat, t(mymat))
                        diag(dyadTable) <- diag(dyadTable) / 2
                        if (valmin == 0 && valmax ==1)
                        {
                            mutual <- dyadTable[2, 2]
                            asymm <- dyadTable[2, 1]
                            nulls <- dyadTable[1, 1]
                            totDyad <- nulls + asymm + mutual

                        }
                    }
                    if (valmin == 0 && valmax == 1)
                    {
                        Report(c(format(per + periodFromStart, width=6),
                                 ".", format(totDyad, width=14),
                                 format(mutual, width=9),
                                 format(asymm, width=10),
                                 format(nulls, width=10), "\n"),
                               sep="", outf)
                    }
                }
                periodFromStart <- periodFromStart + atts$netdims[3]
            }
            if (difficult)
            {
                Report(c("There is a density equal to 0.0 or 1.0.\n",
                         "This may lead to difficulties.\n"), outf)
            }
            Report("\n", outf)
            Report(c("Standard values for initial parameter values\n",
                     "-------------------------------------------------\n\n"),
                   sep = "", outf)
            myeff <- effects[effects$name == netnames[net],]
            myrate <- myeff[myeff$shortName == "Rate",]
            for (i in 1:nrow(myrate))
            {
                Report(c(format(myrate$effectName[i], width=35),
                     format(round(myrate$initialValue[i], 4), nsmall=4,
                            width=10), "\n"), outf)
            }
            myobj <- myeff[myeff$shortName == "density" &
                           myeff$type == "eval",]
            untrimmed <- myobj$untrimmedValue
            Report(c(format(myobj$effectName, width=46),
                     format(round(untrimmed, 4), nsmall=4, width=10),
                     "\n"), outf)
            if (untrimmed > 3)
            {
                Report(c("The initial parameter value is very low.",
                         "It is truncated to -3. \n"), outf)
            }
            else
                if (untrimmed < -3)
                {
                    Report(c("The initial parameter value is very high.",
                             "It is truncated to 3.\n"), outf)
                }
        }
        Report("\n", outf)
    }

    initialBehaviors <- function()
    {
        Report("\n", outf)
        Heading(2, outf,"Dependent actor variables:")
        if (nBehav == 1)
        {
            Report(c(netnames[types == "behavior"], "\n\n"), outf)
            Heading(3, outf, "Marginal distribution")
        }
        else
        {
            Heading(3, outf, "Marginal distributions")
        }
        for (net in which(types=="behavior"))
        {
            if (nBehav > 1)
            {
                Report(c("Dependent actor variable ", net, ": ",
                         netnames[net], "\n"), sep="", outf)
            }
            Report(c(rep(" ", 3*nobs + 7), "Observations\n"), sep="", outf)
            Report(c(format("values", width=16),
                     format(1 : nobs, width=6), "\n"), sep="", outf)
            Report(c(rep(" ", 16), rep("------", nobs), "\n"), sep="", outf)
            periodFromStart <- 0
            bRange <- gpatts$bRange[net]
            minval <- gpatts$behRange[, net][1]
            maxval <- gpatts$behRange[, net][2]
            vals <- matrix(0, ncol=nobs, nrow=bRange + 1)
            minsub <- max(1, minval)
            missings <- rep(0, nobs)
            missingN <- rep(NA, nobs)
            for (group in 1:nData)
            {
                j <- match(netnames[net], names(data[[group]]$depvars))
                if (is.na(j))
                    stop("network names not consistent")
                depvar <- data[[group]]$depvars[[j]]
                atts <- attributes(depvar)
                for (i in 1: atts$netdims[3])
                {
                   # browser()
                    mytab <- table(depvar[, 1, i])
                    vals[as.numeric(names(mytab)) + 1 - minval,
                         periodFromStart + i] <- mytab
                   # vals[factor(names(mytab)), periodFromStart + i] <- mytab
                    missings[periodFromStart + i] <- sum(is.na(depvar[, 1, i]))
                    missingN[periodFromStart + i] <- atts$netdims[1]
                }
                periodFromStart <- periodFromStart + atts$netdims[3]
            }
            for (i in 1:nrow(vals))
            {
                Report(c(format(i - 1 + minval, width=3), rep(" ", 13),
                         format(vals[i, ], width=6), "\n"), sep="", outf)
            }
            if (gpatts$anyMissing[net])
            {
                Report(c(format("missing", width=16),
                          format(missings, width=6), "\n\n"), sep="", outf)
                Report(c(format("(fraction missing"),
                         format(round(missings / missingN, 2), nsmall=2,
                                width=6), "  )\n"), sep="", outf)
            }
            else
            {
                Report("No missings\n", outf)
            }

        }
        Report("\n\n", outf)
        Heading(3, outf, "Changes")
        for (net in which(types=="behavior"))
        {
            if (nBehav > 1)
            {
                Report(c("Dependent actor variable ", net, ": ",
                         netnames[net], "\n"), sep="", outf)
            }
            Report(c(" periods    actors:  down   up   constant  missing  ;",
                     "   steps:   down    up  total\n"), sep="", outf)
            periodFromStart <- 0
            for (group in 1:nData)
            {
                j <- match(netnames[net], names(data[[group]]$depvars))
                if (is.na(j))
                    stop("network names not consistent")
                depvar <- data[[group]]$depvars[[j]]
                atts <- attributes(depvar)
                for (i in 1:(atts$netdims[3] - 1))
                {
                    mytab <- table(depvar[, 1, i + 1] - depvar[, 1, i])
                    numvals <- as.numeric(names(mytab))
                    ups <- mytab[numvals > 0]
                    downs <- mytab[numvals < 0]
                    constants <- mytab[numvals == 0]
                    stepsup <- sum(ups * numvals[numvals > 0])
                    stepsdown <- sum(-1 * downs * numvals[numvals < 0])
                    Report(c(format(i + periodFromStart, width=3),
                             "  =>", format(i + 1 + periodFromStart, width=3),
                             format(sum(downs), width=14),
                             format(sum(ups), width=6),
                             format(sum(constants), width=8),
                             format(sum(is.na(depvar[, 1, i + 1]) |
                                        is.na(depvar[, 1, i])), width=10),
                             format(stepsdown, width=21),
                             format(stepsup, width=6),
                             format(stepsup + stepsdown, width=6), "\n"),
                           sep="", outf)
                }
                Report("\n", outf)
                periodFromStart <- periodFromStart + atts$netdims[3]
            }
            myeff <- effects[effects$name == netnames[net],]
            myobj <- myeff[myeff$shortName == "linear" &
                           myeff$type == "eval",]
            Report(c("For this variable, the standard initial ",
                     "behavioral tendency parameter is ",
                     format(round(myobj$initialValue,4), nsmall=4,width=8),
                     "\n"), sep="",
                   outf)
        }

        Report("\n", outf)
    }

    if (!inherits(data, "sienaGroup"))
    {
        nData <- 1
        data <- sienaGroupCreate(list(data), singleOK=TRUE)
    }
    else
    {
        nData <- length(data)
    }
    Report("\n\n", outf)
    Heading(1, outf, "Initial data description.")
    gpatts <- attributes(data)
    types <- gpatts$types
    netnames <- gpatts$netnames
    nobs <- gpatts$observations + length(data)
    nOneMode <- sum(types == "oneMode")
    nBehav <- sum(types == "behavior")
    if (nOneMode> 0)
    {
        initialNetworks()
    }
    if (nBehav > 0)
    {
        initialBehaviors()
    }
    Report(c("Initialisation of project <<", modelName,
             ">> executed succesfully.\n"), sep="", outf)
}

