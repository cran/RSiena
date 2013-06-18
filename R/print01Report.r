#/******************************************************************************
# * SIENA: Simulation Investigation for Empirical Network Analysis
# *
# * Web: http://www.stats.ox.ac.uk/~snidjers/siena
# *
# * File: print01Report.r
# *
# * Description: This module contains the function to print the initial report
# *****************************************************************************/
##@print01Report Reporting
print01Report <- function(data, myeff, modelname="Siena", session=NULL,
						  getDocumentation=FALSE)
{
	##@reportDataObject1 internal print01Report
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
	##@reportDataObject internal print01Report
	reportDataObject <- function(x, periodFromStart=0, multi=FALSE,
								 session=session)
	{
		##@reportStart internal print01Report
		reportStart <- function()
		{
			multipleNodeSets <- length(x$nodeSets) > 1
			if (multipleNodeSets)
			{
				Report("Dependent variables	  Type		 NodeSet(s) (R, C)\n",
					   outf)
				Report("-------------------	  ----		 -----------------\n",
					   outf)
				for (i in 1:length(x$depvars))
				{
					atts <- attributes(x$depvars[[i]])
					Report(c(format(atts$name, width=20),
							 format(atts$type, width=12)), outf)
					for (j in 1:length(atts$nodeSet))
					{
						if (j > 1)
						{
							Report(', ', outf)
						}
						Report(c(format(atts$nodeSet[j]),
								 " (", atts$netdims[j], ")"), sep="", outf)
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

		##@reportNetworks internal print01Report
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
								 atts$netdims[2], ".\n"), sep="", outf)
					}
					for (k in 1:x$observations)
					{
						if (!is.null(session))
						{
							namesession <-
								session[session$Name == netname, ]
							filename <- namesession$Filename
							if (length(filename) > 1)
							{
								if (namesession$Format[1] == "Siena net")
								{
									period <-
										strsplit(namesession$Period, " ")
									sub <- sapply(period, function(x) k %in% x)
								}
								else
								{
								   period <-
									unlist(strsplit(namesession$Period, " "))
								 sub <- match(k, period)
								}
								filename <- filename[sub]
							}
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
						   # require(Matrix)
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
						if (attr(depvar, "type") == "bipartite")
						{
							tmp <- format(cbind(1:atts$netdims[1], outdeg))
						}
						else
						{
							tmp <- format(cbind(1:atts$netdims[1], outdeg,
												indeg))
						}

						Report(tmp[, 1], fill=60, outf)
						Report("out-degrees\n", outf)
						Report(tmp[, 2], fill=60, outf)
						if (attr(depvar, "type") != "bipartite")
						{
							Report("in-degrees\n", outf)
							Report(tmp[, 3], fill=60, outf)
						}
						## report structural values
						if (attr(depvar, "structural"))
						{
							if (attr(depvar, "sparse"))
							{
								nstruct0 <- sum(depvar[[k]]@x %in% c(10))
								nstruct1 <- sum(depvar[[k]]@x %in% c(11))
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
									nstruct0 <- sum(depvar[[k]]@x %in% c(10))
									nstruct1 <- sum(depvar[[k]]@x %in% c(11))
								}
								else
								{
									nstruct0 <- sum(depvar[, , k] %in% c(10))
									nstruct1 <- sum(depvar[, , k] %in% c(11))
								}
								Report(c('	', nstruct0, ' structural zero'),
									   sep='', outf)
								Report(ifelse(nstruct0 > 1,
											  "s were found (code 10).\n",
											  " was found (code 10).\n"), outf)
								Report(c('	', nstruct1, ' structural one'),
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
									nnonactive <- nnonactive >=
										nrow(depvar[, , k])
								}
								if (sum(nnonactive)	 == 1)
								{
									Report(c("Actor ", which(nnonactive),
											 " is inactive at this ",
											 "observation.\n"), sep='', outf)
								}
								else if (sum(nnonactive) > 1)
								{
									Report(c("Actors", which(nnonactive),
											 "are inactive at this",
											 "observation.\n"), fill=80, outf)
								}
							}
						}
						if (attr(depvar, "sparse"))
						{
							depvark <- depvar[[k]]
							diag(depvark) <- 0
							anymissings <- any(is.na(depvark))
						}
						else
						{
							depvark <- depvar[, , k]
							diag(depvark) <- 0
							anymissings <- any(is.na(depvark))
						}
						if (anymissings)
						{
							Report(c("\nFor observation moment ",
									 k + periodFromStart,
									 ", number of missing values ",
									 "are:\n"),
								   sep="", outf)
							if (attr(depvar, "type") != "bipartite")
							{
								Report("Nodes\n", outf)
								tmp <- format(cbind(1:atts$netdims[1],
													missrow, misscol))
								Report(tmp[, 1], fill=60, outf)
								Report("missing in rows\n", outf)
								Report(tmp[, 2], fill=60, outf)
								Report("missing in columns\n", outf)
								Report(tmp[, 3], fill=60, outf)
								mult <- atts$netdims[1] - 1
							}
							else
							{
								Report("Senders\n", outf)
								tmp <- format(cbind(1:atts$netdims[1],
													missrow))
								Report(tmp[, 1], fill=60, outf)
								Report("missing in rows\n", outf)
								Report(tmp[, 2], fill=60, outf)
								tmp <- format(cbind(1:atts$netdims[2],
													misscol))
								Report("Receivers\n", outf)
								Report(tmp[, 1], fill=60, outf)
								Report("missing in columns\n", outf)
								Report(tmp[, 2], fill=60, outf)
								mult <- atts$netdims[2]
							}
							Report(c("Total number of missing data: ",
									 sum(missrow),
									 ", corresponding to a fraction of ",
									 format(round(sum(missrow)/
												  atts$netdims[1] /
												  mult, 3),
											nsmall=3),
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
		##@reportBehaviors internal print01Report
		reportBehaviors <- function()
		{
			Heading(2, outf, "Reading dependent actor variables.")
			iBehav <- 0
			for (i in 1:length(x$depvars))
			{
				if (types[i] == "behavior")
				{
					depvar <- x$depvars[[i]]
					atts <- attributes(depvar)
					netname <- atts$name

					iBehav <- iBehav + 1
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
						Report(c("But if there is a previous (or later)",
								 "nonmissing value,",
								 "this is used as the imputed value.\n"), outf)
						Report("Modal values:\nObservation	", outf)
						Report(c(format(1:x$observations+periodFromStart,
										width=4), '\n'), outf)
						Report(c(format("Modes", width=12),
								 format(atts$modes, width=4)), outf)
						Report("\n", outf)
					}
					depvar2 <- depvar
					depvar2[is.na(depvar2)] <- 0
					if (!isTRUE(all.equal(as.vector(depvar2) -
										  round(as.vector(depvar2)),
										  rep(0, length(depvar2)))))
					{
						 Report(c("Non-integer values noted in this behavior",
								  "variable: they will be truncated.\n")
								 , outf)
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
					 "		overall\n"), sep="", outf)
			for (i in 1:length(x$depvars))
			{
				if (types[i] == "behavior")
				{
					depvar <- x$depvars[[i]][, 1, ]
					atts <- attributes(x$depvars[[i]])
					netname <- atts$name
					missings <- colSums(is.na(depvar))
					Report(c(format(netname, width=12),
							 format(c(missings, sum(missings)),
									width=10), "	  (",
							 format(round(100 * sum(missings)/
										  nrow(depvar)/ncol(depvar), 1),
									nsmall=1, width=4), ' %)\n'), sep="", outf)
				}
			}
			Report("\nMeans per observation:\n", outf)
			Report(c(" observation", format(1:x$observations+periodFromStart,
											width=10),
					 "		overall\n"), sep="", outf)
			for (i in 1:length(x$depvars))
			{
				if (types[i] == "behavior")
				{
					depvar <- x$depvars[[i]][, 1, ]
					atts <- attributes(x$depvars[[i]])
					netname <- atts$name
					means <- colMeans(depvar, na.rm=TRUE)
					Report(c(format(netname, width=14),
							 format(round(means, 3), nsmall=3,
									width=10), format(round(mean(means),
									3), width=10), '\n'), sep="", outf)
				}
			}
		}
		##@reportConstantCovariates internal print01Report
		reportConstantCovariates <- function()
		{
			nCovars <- length(x$cCovars)
			covars <- names(x$cCovars)
			Heading(2, outf, "Reading constant actor covariates.")
			if (!is.null(session))
			{
				covarssession <- session[session$Type == "constant covariate", ]
				for (i in 1:nrow(covarssession))
				{
					names <- strsplit(covarssession$Name[i],
									  " ", fixed=TRUE)[[1]]
					ncases <- length(x$cCovars[[match(names[1], covars)]])
					Report(c("Covariate data file", covarssession$Filename[i],
						   "with", length(names), "variables,", ncases,
							 "cases, named:\n"), outf)
					Report(paste(names, "\n"), outf, sep="")
				}
			}
			else
			{
				Report(c(nCovars, "variable"),outf)
				Report(ifelse(nCovars == 1, ", named:\n", "s, named:\n"), outf)
				for (i in seq(along=covars))
				{
					Report(c(format(covars[i], width=15), '\n'), outf)
				}
			}
			Report(c("\nA total of", nCovars,
					 "non-changing individual covariate"), outf)
			Report(ifelse(nCovars == 1, ".\n\n", "s.\n\n"), outf)
			Report("Number of missing cases:\n", outf)
			for (i in seq(along=covars))
			{
				Report(c(format(covars[i], width=15),
						 sum(is.na(x$cCovars[[i]])), "	(",
						 format(round(100 * sum(is.na(x$cCovars[[i]]))/
									  length(x$cCovars[[i]]), 1),
								width=3, nsmall=1), '%)\n'), outf)
			}
			Report("\nInformation about covariates:\n", outf)
			Report(c(format("minimum  maximum	  mean", width=38,
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
			if (nData <= 1)
			{
				Report(c("The mean value", ifelse(nCovars == 1, " is", "s are"),
					 " subtracted from the covariate",
					 ifelse(nCovars == 1, ".\n\n", "s.\n\n")), sep="", outf)
			}
		}
		##@reportChangingCovariates internal print01Report
		reportChangingCovariates <- function()
		{
			nCovars <- length(x$vCovars)
			covars <- names(x$vCovars)
			use <- ! covars %in% names(x$cCovars)
			nCovars <- length(x$vCovars[use])
			Heading(2, outf, "Reading exogenous changing actor covariates.")
			if (!is.null(session))
			{
				if (nData > 1)
				{
					covarssession <-
						session[session$Type == "constant covariate", ]
					for (i in 1:nrow(covarssession))
					{
						names <- strsplit(covarssession$Name[i],
										  " ", fixed=TRUE)[[1]]
						ncases <- length(x$vCovars[[match(names[1], covars)]])
						Report(c("Covariate data file",
								 covarssession$Filename[i],
								 "with", length(names), "variables,", ncases,
								 "cases, named:\n"), outf)
						Report(paste(names, "\n"), outf, sep="")
					}
				}
				covarssession <- session[session$Type == "changing covariate", ]
				for (i in seq(along=covarssession[,1]))
				{
					ncases <- nrow(x$vCovars[[match(covarssession$Name[i],
													covars)]])
					Report(c("Exogenous changing covariate ",
							 covarssession$name[i], " read from file ",
							 covarssession$Filename[i], ".\n"), sep="", outf)
					Report(c("Number of cases is ", ncases, ".\n"), sep="",
						   outf)
				}
			}
			else
			{
				Report(c(nCovars, "variable"),outf)
				Report(ifelse(nCovars == 1, ", named:\n", "s, named:\n"), outf)
				for (i in seq(along=covars[use]))
				{
					Report(c(format(covars[use][i], width=15), '\n'), outf)
				}
			}
			Report(c("\nA total of", nCovars,
					 "exogenous changing actor covariate"), outf)
				Report(ifelse(nCovars == 1, ".\n\n", "s.\n\n"), outf)
			Report("Number of missing cases per period:\n", outf)
			Report(c(" period	", format(1:(x$observations - 1) +
										 periodFromStart, width=9),
					 "		 overall\n"), sep="", outf)
			for (i in seq(along=covars))
			{
				if (use[i])
				{
					thiscovar <- x$vCovars[[i]] ## matrix
					misscols <- colSums(is.na(thiscovar))
					Report(c(format(covars[i], width=10),
							 format(misscols, width=8),
							 format(sum(misscols), width=9), "	   (",
							 format(round(100 * sum(misscols)/nrow(thiscovar)/
										  ncol(thiscovar), 1), nsmall=1,
									width=3), '%)\n'), outf)
				}
			}
			Report("\nInformation about changing covariates:\n", outf)
			Report(c(format("minimum  maximum	  mean", width=38,
							justify="right"), "\n"), outf)
			for (i in seq(along=covars))
			{
				if (use[i])
				{
					atts <- attributes(x$vCovars[[i]])
					Report(c(covars[i], '\n'), outf) # name
					for (j in 1:(ncol(x$vCovars[[i]])))
					{
						Report(c("	period", format(j + periodFromStart,
												   width=3),
								 format(round(atts$rangep[1, j], 1),
										nsmall=1, width=7),
								 format(round(atts$rangep[2, j], 1),
										nsmall=1, width=7),
								 format(round(atts$meanp[j], 3),
										nsmall=3, width=10), "\n"), outf)
					}
					Report(c(format("Overall", width=29),
							 format(round(atts$mean, 3), width=10, nsmall=3),
							 "\n"), outf)

				}
			}
			if (nData <= 1)
			{
				Report("\nThe overall mean value", outf)
				Report(c(ifelse(nCovars	 == 1, " is", "s are"),
					 "subtracted from the covariate.\n\n"),	 outf)
			}
		}
		##@reportConstantDyadicCovariates internal print01Report
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
				if (attr(x$dycCovars[[i]], "sparse"))
				{
					myvar <- x$dycCovars[[i]][[1]]
				}
				else
				{
					myvar <- x$dycCovars[[i]]
				}
				diag(myvar) <- 0
				Report(c(format(covars[i], width=15),
						 sum(is.na(myvar)), "  (",
						 format(round(100 * sum(is.na(myvar))/
									  (length(myvar) - nrow(myvar)), 1),
								width=3, nsmall=1), '%)\n'), outf)
			}
			Report("\nMeans of	covariates:\n", outf)
			Report(c(format("minimum  maximum	  mean", width=38,
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
		##@reportChangingDyadicCovariates internal print01Report
		reportChangingDyadicCovariates <- function()
		{
			covars <- names(x$dyvCovars)
			use <- ! covars %in% names(x$dycCovars) ## need an attributes to say
			nCovars <- length(x$dyvCovars[use])
			Heading(2, outf, "Reading exogenous dyadic covariates.")
			## Report(c("Note that no missing values are considered yet for",
			##			"changing dyadic covariates.\n"), outf)
			for (i in seq(along=covars))
			{
				Report(c("Exogenous dyadic covariate named ", covars[i], '.\n'),
					   sep="", outf)
			}
			Report("Number of tie variables with missing data per period:\n",
				   outf)
			Report(c(" period	", format(1:(x$observations - 1) +
										  periodFromStart, width=9),
					 "		 overall\n"), sep="", outf)
			for (i in seq(along=covars))
			{
				if (use[i])
				{
					sparse <- attr(x$dyvCovars[[i]], "sparse")
					vardims <- attr(x$dyvCovars[[i]], "vardims")
					thiscovar <- x$dyvCovars[[i]] ## array/list of sparse mats
					if (!sparse)
					{
						missvals <- colSums(is.na(thiscovar), dims=2)
					}
					else
					{
						missvals <- sapply(thiscovar, function(x)sum(is.na(x)))
					}
					Report(c(format(covars[i], width=10),
							 format(missvals, width=8),
							 format(sum(missvals), width=9), "	   (",
							 format(round(100 * sum(missvals)/vardims[1]/
										  vardims[2]), nsmall=1,
										  width=3), '%)\n'), outf)
				}
			}
			Report("\nMeans of	covariates:\n", outf)
			for (i in seq(along=covars))
			{
				atts <- attributes(x$dyvCovars[[i]])
				Report(c(format(covars[i], width=10),
						 format(round(atts$mean, 3),
								nsmall=3, width=10), "\n"), outf)
			}
			Report('\n', outf)
			Report(c(ifelse(nCovars == 1, "This ", "These "),
					 "global mean value",
					 ifelse(nCovars == 1, " is ", "s are "),
					 "subtracted from the changing dyadic covariate",
					 ifelse(nCovars ==1, ".\n\n", "s.\n\n")), sep="", outf)
		}
		##@reportCompositionChange internal print01Report
		reportCompositionChange <- function()
		{
			comps <- x$compositionChange
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
			(nData ==1	&& length(x$vCovars) > 0))
		{
			reportChangingCovariates()
		}
		if (length(x$dycCovars) > 0)
		{
			reportConstantDyadicCovariates()
		}
		if (nData > 1 && length(x$dyvCovars) > length(x$dycCovars) ||
			(nData ==1	&& length(x$dyvCovars) > 0))
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
	if (getDocumentation)
	{
		tt <- getInternals()
		return(tt)
	}
	Report(openfiles=TRUE, type="w", projname=modelname)
	Report("							************************\n", outf)
	Report(c("									 ", modelname, ".out\n"),
		   sep='', outf)
	Report("							************************\n\n", outf)
	Report(c("Filename is ", modelname, ".out.\n\n"), sep="", outf)
	Report(c("This file contains primary output for SIENA project <<",
		modelname, ">>.\n\n"), sep="", outf)
	Report(c("Date and time:", format(Sys.time(), "%d/%m/%Y %X"), "\n\n"), outf)
	packageValues <- packageDescription(pkgname, fields=c("Version", "Date"))
	rforgeRevision <-  packageDescription(pkgname,
										  fields="Repository/R-Forge/Revision")
	if (is.na(rforgeRevision))
	{
		revision <- ""
	}
	else
	{
		revision <- paste(" R-forge revision: ", rforgeRevision, " ", sep="")
	}
	Report(c("SIENA version ", packageValues[[1]], " (",
		format(as.Date(packageValues[[2]]), "%d %m %Y"), ")",
			 revision, "\n\n"), sep="", outf)

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
			thisSession <- session[session$Group == names(data)[i],]
			reportDataObject(data[[i]], periodFromStart, multi=TRUE,
							 session=thisSession)
			periodFromStart <- periodFromStart + data[[i]]$observations
	   }
	}
	else
	{
		Heading(1, outf, "Data input.")
		reportDataObject(data[[1]], 0, multi=FALSE, session=session)
	}
	atts <- attributes(data)
	nets <- atts$types != "behavior"
	if (length(data) > 1)
	{
		Heading(1, outf, "Further processing of multi-group data.")
		Report("Series of observations for the multi-group project:\n", outf)
		periodFromStart <- 0
		for (i in seq(along=data))
		{
			Report(c(format(1:data[[i]]$observations + periodFromStart), '\n'),
				   outf)
			periodFromStart <- periodFromStart + data[[i]]$observations
		}
		Report("\n", outf)
		if (length(atts$vCovars) == 1) 
		{
			Report(c("The overall mean value ",
				format(round(atts$vCovarMean, 4), nsmall=3, width=12), 
					 " is subtracted from covariate ", atts$vCovars, 
					  ".\n\n"), sep="", outf)
		}
		else if (length(atts$vCovars) >= 2) 
		{
			Report(c("The mean values are subtracted from the covariates:\n"), outf)
			for (i in seq(along=atts$vCovars))
			{
				Report(c(format(atts$vCovars[i], width=15), 
				format(round(atts$vCovarMean[i], 4), nsmall=3, width=12), '\n'), outf)
			}
		}		
	}
	periodNos <- attr(data, "periodNos")
	if (any(atts$anyUpOnly[nets]))
	{
		netnames <- atts$netnames[nets]
		upOnly <- atts$anyUpOnly[nets]
		for (i in which(upOnly))
		{
			if (sum(nets) > 1)
			{
				Report(c("Network ", netnames[i], ":\n"), sep = "", outf)
			}
			if (atts$allUpOnly[i])
			{
				Report("All network changes are upward.\n", outf)
				Report("This will be respected in the simulations.\n", outf)
				Report("Therefore, there is no outdegree parameter.\n\n", outf)
			}
			else
			{
				Report(c("All network changes are upward for the following",
						 "periods:\n"), outf)
				periodsUp <- unlist(lapply(data, function(x)
					{
						attr(x$depvars[[match(netnames[i], names(x$depvars))]],
							"uponly")
					}))
				periods <- periodNos[c(1:length(periodsUp))[periodsUp]]
				Report(paste(periods, " => ", periods + 1, ";",
							 sep=""), fill=80, outf)
				Report("This will be respected in the simulations.\n\n", outf)
			}
		}
	}
	if (any(atts$anyDownOnly[nets]))
	{
		netnames <- atts$netnames[nets]
		downOnly <- atts$anyDownOnly[nets]
		for (i in which(downOnly))
		{
			if (sum(nets) > 1)
			{
				Report(c("Network ", netnames[i], "\n"), sep = "", outf)
			}
			if (atts$allDownOnly[i])
			{
				Report("All network changes are downward.\n", outf)
				Report("This will be respected in the simulations.\n", outf)
				Report("Therefore, there is no outdegree parameter.\n\n", outf)
			}
			else
			{
				periodsDown <-
					unlist(lapply(data, function(x)
					   {
						   attr(x$depvars[[match(netnames[i],
												 names(x$depvars))]],
								"downonly")
					   }))
				Report(c("All network changes are downward for the",
						 "following periods:\n"), outf)
				periods <- periodNos[c(1:length(periodsDown))[periodsDown]]
				Report(paste(periods, " => ", periods + 1, ";",
							 sep=""), fill=80, outf)
				Report("This will be respected in the simulations.\n\n", outf)
			}
		}
	}
	if (any(atts$anyUpOnly[!nets]))
	{
		netnames <- atts$netnames[!nets]
		upOnly <- atts$anyUpOnly[!nets]
		for (i in which(upOnly))
		{
			Report(c("\nBehavior variable ", netnames[i], ":\n"), sep = "",
				   outf)
			if (atts$allUpOnly[i])
			{
				Report("All behavior changes are upward.\n", outf)
				Report("This will be respected in the simulations.\n", outf)
				Report("Therefore, there is no linear shape parameter.\n\n",
						outf)
			}
			else
			{
				Report(c("All behavior changes are upward for the following",
						 "periods:\n"), outf)
				periodsUp <-
					sapply(data, function(x)
					   {
						   attr(x$depvars[[match(netnames[i],
												 names(x$depvars))]],
								"uponly")
					   })
				periods <- periodNos[c(1:length(periodsUp))[periodsUp]]
				Report(paste(periods, " => ", periods + 1, ";",
							 sep=""), fill=80, outf)
				Report("This will be respected in the simulations.\n\n", outf)
			}
		}
	}
	if (any(atts$anyDownOnly[!nets]))
	{
		netnames <- atts$netnames[!nets]
		downOnly <- atts$anyDownOnly[!nets]
		for (i in which(downOnly))
		{
			Report(c("\nBehavior ", netnames[i], ":\n"), sep = "", outf)
			if (atts$allDownOnly[i])
			{
				Report("All behavior changes are downward.\n", outf)
				Report("This will be respected in the simulations.\n", outf)
				Report("Therefore, there is no linear shape parameter.\n\n",
						outf)
			}
			else
			{
				periodsDown <-
					sapply(data, function(x)
					   {
						   attr(x$depvars[[match(netnames[i],
												 names(x$depvars))]],
								"downonly")
					   })
				Report(c("All behavior changes are downward for the",
						 "following periods:\n"), outf)
				periods <- periodNos[c(1:length(periodsDown))[periodsDown]]
				Report(paste(periods, " => ", periods + 1, ";",
							 sep=""), fill=80, outf)
				Report("This will be respected in the simulations.\n\n", outf)
			}
	   }
	}
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
	 if (any(atts$anyMissing[!nets]))
	{
		netnames <- atts$netnames[!nets]
		missings <- atts$anyMissing[!nets]
		for (i in seq(along=netnames[missings]))
		{
			Report(c("There are missing data for behavior variable ",
					 netnames[i], ".\n"), sep = "", outf)
		}
	}
   Report("\n", outf)

	if (sum(atts$types == 'oneMode') > 0)
	{
		netnames <- atts$netnames[nets]
		if (nData > 1)
		{
			balmean <-
				lapply(data, function(x)
					   sapply(x$depvars, function(y) attr(y, "balmean")))
		}
		else
		{
			balmean <- atts$"balmean"
		}
		if (nData > 1 || sum(atts$types == "oneMode") > 1)
		{
			Report(c("The mean structural dissimilarity values subtracted",
					 "in the\n"), outf)
			Report("balance calculations are\n", outf)
		}
		else
		{
			Report(c("The mean structural dissimilarity value subtracted",
					 "in the\n"), outf)
			Report("balance calculations is ", outf)
		}
		for (i in seq(along=atts$types))
		{
			if (atts$types[i] == "oneMode")
			{
				if (nData > 1)
				{
					thisbalmean <- sapply(balmean, function(x)x[[netnames[i]]])
					##	if (sum(atts$types == "oneMode") > 1)
					if (sum(atts$types != "behavior") > 1)
					{
						Report(c("for network ", netnames[i],":"), sep="",
							   outf)
					}
					Report("\n", outf)
					mystr <- format(paste("Subproject ", 1:nData, " <",
								  atts$names, "> ", sep=""))
					for (j in seq(along=thisbalmean))
					{
						Report(c(mystr[j], ": ",
								 format(round(thisbalmean[j], 4), nsmall=4,
										width=14), "\n"), sep="", outf)
					}
				}
				else
				{
					##	if (sum(atts$types == "oneMode") > 1)
					if (sum(atts$types != "behavior") > 1)
					{
						Report(c("for network ", format(netnames[i], width=12),
								 format(round(balmean[i], 4),
										nsmall=4, width=14), '.\n'),
							   sep="", outf)
					}
					else
					{
						Report(c(format(round(balmean[i], 4), nsmall=4,
										width=14), '.\n'), sep="", outf)
					}
				}
			}
		}
	}
	if (sum(atts$types == "behavior") > 0 ||
		(nData ==1 && length(atts$cCovars) > 0) ||
		length(atts$vCovars) > 0)
	{
		netnames <- atts$netnames
		vCovarSim2 <-
			lapply(data, function(x)
				   lapply(x$vCovars, function(y) attr(y, "simMeans")))
		behSim2 <-
			lapply(data, function(x)
				   lapply(x$depvars, function(y) attr(y, "simMeans")))
		cCovarSim2 <-
			lapply(data, function(x)
				   lapply(x$cCovars, function(y) attr(y, "simMeans")))
		if (nData > 1)
		{
			vCovarSim <-
				lapply(data, function(x)
					   sapply(x$vCovars, function(y) attr(y, "simMean")))
			behSim <-
				lapply(data, function(x)
					   sapply(x$depvars, function(y) attr(y, "simMean")))
		}
		else
		{
			vCovarSim <- atts$"vCovarSim"
			behSim <- atts$"bSim"
		}
		Report(c("\nFor the similarity variable calculated from each actor",
				 "covariate,\nthe mean is subtracted.\nThese means are:\n"),
			   outf)
		if (nData == 1) ## ie we may have constant covariates
		{
			for (i in seq(along=atts$cCovars))
			{
				if (atts$cCovarPoszvar[i])
				{
					Report(c("Similarity", format(atts$cCovars[i], width=12),
							 ':', format(round(atts$cCovarSim[i], 4), width=12,
										 nsmall=4), '\n'), outf)
				}
			}
		}
		for (i in seq(along=atts$netnames))
		{
			if (atts$types[i] == "behavior" && atts$bPoszvar[i])
			{
				if (nData > 1)
				{
					thisSim <- sapply(behSim, function(x)x[[netnames[i]]])
					Report(c("Similarity ", format(atts$netnames[i], width=12),
							 ":\n"), sep="", outf)
					mystr <- format(paste("	 Subproject ", 1:nData, " <",
								  atts$names, "> ", sep=""))
					for (j in seq(along=thisSim))
					{
						Report(c(mystr[j], format(round(thisSim[j], 4),
												  nsmall=4, width=12), "\n"),
							   sep="", outf)
					}
					Report("\n", outf)
				}
				else
				{
					Report(c("Similarity", format(atts$netnames[i], width=12),
							 ':', format(round(atts$bSim[i], 4), nsmall=4,
										 width=12), '\n'), outf)
				}
			}
		}
		for (i in seq(along=atts$vCovars))
		{
			covarnames <- atts$vCovars
			if (atts$vCovarPoszvar[i])
			{
				if (nData > 1)
				{
					thisSim <- sapply(vCovarSim, function(x)x[[covarnames[i]]])
					Report(c("Similarity ", format(covarnames[i], width=12),
							 ":\n"), sep="", outf)
					mystr <- format(paste("	 Subproject ", 1:nData, " <",
										  atts$names, "> ", sep=""))
					for (j in seq(along=thisSim))
					{
						Report(c(mystr[j], format(round(thisSim[j], 4),
												  nsmall=4, width=12), "\n"),
							   sep="", outf)
					}
					Report("\n", outf)
				}
				else
				{
					Report(c("Similarity", format(atts$vCovars[i], width=12),
							 ':', format(round(atts$vCovarSim[i], 4), width=12,
										 nsmall=4), '\n'), outf)
				}
			}
		}
		## report on alter similarity means
		Report(c("\nFor the distance-2 similarity variable calculated from",
				 "each actor \ncovariate, the mean is also subtracted.\n",
				 "These means are:\n"), outf)
		if (nData == 1) ## ie we may have constant covariates
		{
			sims <- cCovarSim2[[1]]
			for (i in seq(along=sims))
			{
				if (atts$cCovarPoszvar[i])
				{
					Report(c("Distance-2 Similarity ",
							 format(atts$cCovars[i], width=12),
							 ':', format(names(sims[[i]]), width=12,
										 justify="right"), '\n'), sep="", outf)
					Report(c(rep(" ",  35),
							 format(round(sims[[i]], 4),  width=12, nsmall=4),
							 '\n'), sep="", outf)
				}
			}
		}
		for (i in seq(along=atts$netnames))
		{
			if (atts$types[i] == "behavior" && atts$bPoszvar[i])
			{
				if (nData > 1)
				{
					thisSim <- lapply(behSim2, function(x)x[[netnames[i]]])
					Report(c("Distance-2 Similarity ",
							 format(atts$netnames[i], width=12), ":\n"),
						   sep="", outf)
					Report(c(rep(" ", 24),
							 format(names(thisSim[[1]]),
									width=12, justify="right"), "\n"), sep="",
						   outf)
					mystr <- format(paste("	 Subproject ", 1:nData, " <",
								  atts$names, "> ", sep=""))
					for (j in seq(along=thisSim))
					{
						Report(c(mystr[j], format(round(thisSim[[j]], 4),
												  nsmall=4, width=12), "\n"),
							   sep="", outf)
					}
					Report("\n", outf)
				}
				else
				{
					sims <- behSim2[[1]]
					Report(c("Distance-2 Similarity ", format(atts$netnames[i],
															 width=12),
							 ':', format(names(sims[[i]]), width=12), '\n'),
						   sep="", outf)
					Report(c(rep(" ",  35), format(round(sims[[i]], 4), width=12,
										 nsmall=4), '\n'), sep="", outf)
				}
			}
		}
		for (i in seq(along=atts$vCovars))
		{
			covarnames <- atts$vCovars
			if (atts$vCovarPoszvar[i])
			{
				if (nData > 1)
				{
					thisSim <- lapply(vCovarSim2, function(x)x[[covarnames[i]]])
					Report(c("Distance-2 Similarity ", format(covarnames[i],
															  width=12),
							 ":\n"), sep="", outf)
					Report(c(rep(" ", 24), format(names(thisSim[[1]]),
												  width=12, justify="right"),
							 "\n"),	 sep="", outf)
					mystr <- format(paste("	 Subproject ", 1:nData, " <",
										  atts$names, "> ", sep=""))
					for (j in seq(along=thisSim))
					{
						Report(c(mystr[j], format(round(thisSim[[j]], 4),
												  nsmall=4, width=12), "\n"),
							   sep="", outf)
					}
					Report("\n", outf)
				}
				else
				{
					sims <- vCovarSim2[[1]]
					Report(c("Distance-2 Similarity", format(atts$vCovars[i],
															 width=12),
							 ':', format(names(sims[[i]]), width=12,
										 justify="right"),
							 '\n'), outf)
					Report(c(rep(" ",  35), format(round(sims[[i]], 4), width=12,
												   nsmall=4), '\n'), sep="",
						   outf)
				}
			}
		}
	}
	## report on constraints
	if (any(atts$anyHigher) || any(atts$anyDisjoint) || any(atts$anyAtLeastOne))
	{
		Report("\n", outf)
		highers <- atts[["anyHigher"]]
		disjoints <- atts[["anyDisjoint"]]
		atleastones <- atts[["anyAtLeastOne"]]
		if (any(highers))
		{
			higherSplit <- strsplit(names(highers)[highers], ",")
			lapply(higherSplit, function(x)
			   {
				   Report(c("Network ", x[1], " is higher than network ", x[2],
							".\n"), sep="", outf)
				   Report("This will be respected in the simulations.\n\n",
						  outf)
			  })
		}
		if (any(disjoints))
		{
			disjointSplit <- strsplit(names(disjoints)[disjoints],',')
			lapply(disjointSplit, function(x)
			   {
				   Report(c("Network ", x[1], " is disjoint from network ",
							x[2], ".\n"), sep="", outf)
				   Report("This will be respected in the simulations.\n\n",
						  outf)
			  })
		}
		if (any(atleastones))
		{
			atLeastOneSplit <- strsplit(names(atleastones)[atleastones],',')
			lapply(atLeastOneSplit, function(x)
			   {
				   Report(c("A link in at least one of networks ",
							x[1], " and", x[2],
						   " always exists.\n"), sep="", outf)
				   Report("This will be respected in the simulations.\n\n",
						  outf)
			  })
		}
	}
	printInitialDescription(data, myeff, modelname)
	##close the files
	Report(closefiles=TRUE)
}













