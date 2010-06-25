##*****************************************************************************
## * SIENA: Simulation Investigation for Empirical Network Analysis
## *
## * Web: http://www.stats.ox.ac.uk/~snidjers/siena
## *
## * File: sienaTimeTest.r
## *
## * Description: This file contains the sienaTimeTest code for testing the
## * significance of additional time dummies interacted with effects, and
## * sienaTimeFix, which is called to set up time dummy interacted effects.
## ****************************************************************************/
##@sienaTimeTest siena07 Does test for time homogeneity of effects
sienaTimeTest <- function (sienaFit, effects=NULL, condition=FALSE)
{
	observations <- sienaFit$f$observations
	# There must be more than 2 observations to do a time test!
	if (observations <=2)
	{
		stop("You must have at least three time periods to test
			 for non-heterogeneity across time.")
	}
	## Screen out the undesired effects
	if (!is.null(effects)) {
		escreen = setdiff(1:nrow(sienaFit$effects), effects)
	} else {
		escreen = 99999
	}
	# Identify which effects are rate parameters
	indRateEffects <- which(sienaFit$effects$shortName[-escreen]=="Rate")
	# Identify which effects are estimated dummy terms
	indDummiedEffects <- grep("Dummy", sienaFit$effects$effectName[-escreen])
	# The effects which will be tested are stored here. Take all of the
	# effects, and take out the rate and dummy effects. These indices
	# will have to be changed for the moments, scores, etc. after we
	# screen the sienaFit ingredients.
	indBaseEffects <- setdiff(1:nrow(sienaFit$effects[-escreen,]), c(indRateEffects,
														  indDummiedEffects))
	baseNames=sienaFit$effects[-escreen,]$effectName[indBaseEffects]
	# toTest will hold booleans for which of the time dummies have not
	# been estimated and thus are candidates for time tests.
	toTest <- array(TRUE, dim=c(length(indBaseEffects),
								sienaFit$f$observations - 2))
	rownames(toTest) <- sienaFit$effects[-escreen,]$effectNumber[indBaseEffects]
	colnames(toTest) <- 2:(sienaFit$f$observations - 1)
	# dummyByEffect gets passed from sienaTimeTest so that other functions
	# know which dummies belong to which base effects.
	dummyByEffect <- array(0, dim=c(length(indBaseEffects),
									sienaFit$f$observations - 2))
	dimnames(dummyByEffect) <- dimnames(toTest)
	# dscreen is the first important screening vector, which will determine
	# which egoX dummies are fixed, so that we do not consider them as included
	dscreen <- which(sienaFit$effects[-escreen,]$shortName=='egoX' &
					sienaFit$effects[-escreen,]$fix
					 & length(grep("Dummy",
					sienaFit$effects[-escreen,]$effectName)) > 0)
	if (length(dscreen)==0)
	{
		dscreen <- 99999
	}
	## If the estimation was unconditional, the rate parameters will have scores
	## and moments which must also be screened out. Ruth, is there a simple way to
	## check conditioning? I tried $conditional and it doesnt seem to do what I
	## intuitively expected. For now, I just check the dimensionality of the scores,
	## as it will match the number of included "effects" on dimension 3 if uncond.
	## estimation was used.
	if (dim(sienaFit$sf2[,,-escreen])[3] == dim(sienaFit$effects[,,-escreen])[1]) {
		rscreen <- indRateEffects
	} else {
		rscreen <- 99999
	}
	## Go through each effect which had a time dummy included, and incorporate this
	## information into the toTest vector. i.e. if a time dummy was estimated, set
	## its element in toTest equal to FALSE so that we do not time test it
	for (i in sienaFit$effects[-escreen,]$effectNumber
		[sienaFit$effects[-escreen,]$timeDummy != ',']){
		tmp <- toString(sienaFit$effects[-escreen,]$timeDummy[
					   sienaFit$effects[-escreen,]$effectNumber == i])
		tmp <- strsplit(tmp, split=",", fixed=TRUE)[[1]]
		if (length(which(!tmp == '')) > 0)
		{
			## The effect we are looking at is a time dummy.
			if (tmp[1]=='isDummy' & !(i %in% sienaFit$effects[-escreen,]$
									  effectNumber[dscreen]))
			{
				## Dont test this dummy...
				toTest[rownames(toTest)==as.numeric(tmp[3]),
					colnames(toTest)==as.numeric(tmp[2])] <- FALSE
				## We want to be able to reference this effect given an
				## index for the base effect and a time period, so store
				## this information in dummyByEffect -- this is used
				## extensively in plot.sienaTimeTest
				dummyByEffect[rownames(toTest)==as.numeric(tmp[3]),
					colnames(toTest)==as.numeric(tmp[2])]  <-
					which(sienaFit$effects[-escreen,]$
					effectNumber[-c(rscreen,dscreen)]==i)
			}

		}
		else
		{
			## The effect we are looking at had a time dummy,
			## nothing required for now.
			next
		}
	}
	##  nEffects, nSims, nameslist, nDummies convert commonly used ingredients
	##  from sienaFit into an easily accessed form based on the screens
	##  set up above
	nEffects <- length(indBaseEffects) + sum(!toTest)
	## With the use of multiple nodes, sometimes the sienaFit object comes back
	## with the wrong number of iterations!! Fixing it by looking elsewhere:
	## Used to be: nSims <- sienaFit$n3
	nSims <- dim(sienaFit$sf2[,,-escreen])[1]
	nameslist <- list(
					Iteration=paste("it", 1:nSims, sep=""),
					Wave=paste("Wave", 1:(observations - 1), sep=""),
					Effect=sienaFit$effects[-escreen,]$effectName[-c(dscreen,rscreen)]
					)
	nDummies <- sum(toTest)
	nTotalEffects <- nDummies + nEffects
	## obsStats, moment, scores are the crucial ingredients from sienaFit which
	## screen for the base effects and make the rest of the code clean
	obsStats <- t(sienaFit$targets2[-c(dscreen,rscreen,escreen), ])
	moment <- sienaFit$sf2[, , -c(dscreen,rscreen,escreen)] - rep(obsStats, each=nSims)
	scores <- sienaFit$ssc[ , , -c(dscreen,rscreen,escreen)]
	## Because the sienaFit object does not have a strict class definition,
	## the $sf2 and $targets2 arrays cannot be expected to always have the
	## proper format. The best we can do is therefore to die gracefully if
	## the arrays do not line up:
	G <- array(0, dim=c(nSims, observations - 1, nEffects + nDummies))
	SF <- array(0, dim=c(nSims, observations - 1, nEffects + nDummies))
	if (sum(dim(G[, , 1:nEffects]) != dim(moment))+
		sum(dim(SF[, , 1:nEffects]) != dim(scores))>0) {
		stop("The moments and scores in your sienaFit have unexpected dimensions.\n
			It is possible that your model specifications are not yet implemented\n
			in sienaTimeTest. Please contact the developers.\n\nDid you include
			the base effect?\n")
	}
	## Will be used to construct the dummy names for output
	dummyNames <- rep("", nDummies)
	## Set the base effects G equal to the moments from sienaFit
	G[, , 1:nEffects] <- moment
	## inc used for incrementing through the dummies
	inc <- nEffects
	for (i in 1:nrow(toTest))
	{
		for (j in 1:ncol(toTest))
		{
			## Go through each dummy to be tested
			if (toTest[i, j])
			{
				inc <- inc + 1
				## And add scores and moments for the specific time period j+1
				G[, j + 1, inc] <- moment[, j + 1, i]
				dummyNames[inc-nEffects] <- paste("(*)Dummy", j + 1, ":",
												  nameslist$Effect[i], sep="")
			}
		}
	}
	## Put names onto G for easy reference
	dimnames(G) <- list(nameslist$Iteration, nameslist$Wave,
						c(nameslist$Effect, dummyNames))
	## Make the covariance matrix for the new moments
	sigma <- cov(apply(G, c(1, 3), sum))
	## Basically repeat this process for the scores:
	SF[, , 1:nEffects] <- scores
	inc <- nEffects
	dummyProps <- list()
	for (i in 1:nrow(toTest))
	{
		for (j in 1:ncol(toTest))
		{
			if (toTest[i, j])
			{
				inc <- inc + 1
				SF[, j + 1, inc] <- scores[, j + 1, i]
				## Save some information on these dummies for later;
				## these operations dont relate directly to the scores
				dummyByEffect[i, j]=inc
				dummyProps$shortName[inc] <- sienaFit$effects[-escreen,]$shortName[i]
				dummyProps$interaction1[inc] <- sienaFit$effects[-escreen,]$interaction1[i]
				dummyProps$type[inc] <- sienaFit$effects[-escreen,]$type[i]
				dummyProps$period[inc] <- j + 1
			}
		}
	}
	## Copy the dimnames for G
	dimnames(SF) <- dimnames(G)
	## We have now set up all of the ingredients properly, so we may proceed with
	## the score type test of Schweinberger (2007)
	D <- derivativeFromScoresAndDeviations(SF, G)
	fra <- apply(G, 3, sum) / nSims
	doTests <- c(rep(FALSE, nEffects), rep(TRUE, nDummies))
	jointTest <- ScoreTest(nTotalEffects, D, sigma, fra, doTests, maxlike=FALSE)
	jointTestP <- 1 - pchisq(jointTest$testresOverall, nDummies)
	if (! condition) {
		individualTest <- jointTest$testresulto[1:nDummies]
	} else {
		individualTest <- sapply(1:nDummies, function (i)
			{ doTests <- rep(FALSE, nEffects + nDummies)
				doTests[nDummies+i] <- TRUE
				test <- ScoreTest(nTotalEffects, D, sigma, fra, doTests, FALSE)
				test$testresulto[1]
			})
	}
	individualTestP <- 2 * (1-pnorm(abs(individualTest))[1:nDummies])
	rownames(jointTestP) <- c("Joint Significant Test")
	colnames(jointTestP) <- c("p-Val")
	thetaOneStep <- c(sienaFit$theta[-c(dscreen,rscreen,escreen)], rep(0, nDummies)) +
			jointTest$oneStep
	effectTest <- sapply(1:length(indBaseEffects), function (i)
					{
						 doTests <- rep(FALSE, nEffects + nDummies)
						 tmp <- which(dummyProps$shortName ==
									  sienaFit$effects[-escreen,]$shortName[i] &
									  dummyProps$interaction1 ==
									  sienaFit$effects[-escreen,]$interaction1[i])
						 if (length(tmp) > 0)
						 {
							doTests[tmp] <- TRUE
							test <- ScoreTest(nTotalEffects, D, sigma, fra,
										   doTests, FALSE)
							test$testresOverall
						 }
						 else
						 {
							NA
						 }
					})

	dim(effectTest) <- c(length(indBaseEffects), 1)
	effectTestP <- round(1 - pchisq(effectTest, apply(toTest, 1, sum)), 5)
	rownames(effectTestP) <- baseNames
	colnames(effectTestP) <- c("p-Val")
	thetaStar <- cbind(c(sienaFit$theta[-c(dscreen,rscreen,escreen)], rep(0, nDummies)),
					   thetaOneStep,
					   round(c(2-2 * pnorm(abs(sienaFit$theta[-c(dscreen,rscreen,escreen)]/
											 sqrt(diag(sienaFit$covtheta)[-c(dscreen,
											rscreen,escreen)]))),
							   individualTestP), 5))
	colnames(thetaStar) <- c("Initial Est.", "One Step Est.", "p-Value")
	rownames(thetaStar) <- dimnames(G)[[3]]
	returnObj <- list(
					  JointTest=jointTestP,
					  EffectTest=effectTestP,
					  IndividualTest=thetaStar,
					  JointTestStatistics=jointTest,
					  EffectTestStatistics=effectTest,
					  IndividualTestStatistics=individualTest,
					  CovDummyEst=jointTest$covMatrix,
					  Moments=G,
					  NonRateIndices=indBaseEffects,
					  Waves=dim(G)[2],
					  Sims=dim(G)[1],
					  Effects=dim(G)[3],
					  DummyIndexByEffect=dummyByEffect,
					  DummyStdErr=sqrt(diag(jointTest$covMatrix)),
					  OriginalEffects=nEffects,
					  OriginalThetaStderr=sqrt(diag(sienaFit$covtheta))[-c(dscreen,
									  		rscreen,escreen)],
					  SienaFit=sienaFit,
					  DummyProps=dummyProps,
					  ToTest=toTest,
					  ScreenedEffects=setdiff(c(rscreen,escreen),99999)
					  )
	class(returnObj) <- "sienaTimeTest"
	returnObj
}
summary.sienaTimeTest <- function(object, ...)
{
	if (!inherits(object, "sienaTimeTest"))
	{
		stop("not a legitimate Siena time test object")
	}
	class(object) <- c("summary.sienaTimeTest", class(object))
	object
}
print.summary.sienaTimeTest <- function(x, ...)
{
	if (!inherits(x, "summary.sienaTimeTest"))
	{
		stop("not a legitimate Siena time test summary object")
	}
	print.sienaTimeTest(x)
## Additional output to the print will go in here:
	cat("\nIndividual significance tests and one-step estimators:\n")
	print(x$IndividualTest)
	cat("\nParameter-wise joint significance tests (i.e. each
		parameter across all dummies):\n")
	print(x$EffectTest)
	if (x$Waves <=2)
	{
		cat("\n\nNote that these parameter-wise tests have a different
			form than the individual tests, thus testing with 3 observations
			may yield different individual and parameter-wise values.\n\n")
	}
	tmp <- paste(" (", 1:length(rownames(x$IndividualTest)), ") ",
				 rownames(x$IndividualTest), "\n", sep="")
	cat("\nUse the following indices for plotting:\n", tmp)
	tmp <- paste(" (", 1:length(x$NonRateIndices), ") ",
				 rownames(x$IndividualTest)[x$NonRateIndices], "\n", sep="")
	cat("\nIf you would like to fit time dummies to your model, use the
		timeDummy column in your effects object.")
	cat("\nType \"?sienaTimeTest\" for more information on this output.\n")
	invisible(x)
}
print.sienaTimeTest <- function(x, ...)
{
	if (!inherits(x, "sienaTimeTest"))
	{
		stop("not a legitimate Siena time test object")
	}
	effectNames <- rownames(x$IndividualTest)
	dummies <- grepl("Dummy", effectNames)
	dummyIndex <- paste(" (", 1:sum(dummies), ") ", effectNames[dummies],
						"\n", sep="")
	cat("Joint significance test of the dummy parameters:\np-Val = ",
		x$JointTest,
		", \nWhere H0: The following parameters are zero:\n",
		dummyIndex
		)
	invisible(x)
}
plot.sienaTimeTest <- function(x, pairwise=FALSE, effects=1:2,
	dims=c(2, 1), scale=.2, plevels=c(.1, .05, .025),
	multiplot=FALSE, ...)
{
	require(lattice)
	timetest <- x
	if (pairwise)
	{
## On a "pairwise" call, print a pairwise plot of moments
		if (length(intersect(effects, 1:timetest$OriginalEffects))!=
			length(effects))
		{
			cat("Detected an error with the effects included. For a
				parameter-plot, use the following indices:")
			tmp <- paste(" (", 1:length(rownames(timetest$IndividualTest)),
						 ") ", rownames(timetest$IndividualTest), "\n", sep="")
			cat("\nUse the following indices for plotting the pairwise
				moment correlations:\n", tmp)
			stop(" ")
		}
		if (length(effects)==0)
		{
			x <- timetest$Moments

		}
		else
		{
			if (class(effects)!="integer")
			{
				stop("Effects is not a vector of integers.")
			}
			x <- timetest$Moments[, , effects]
		}
		panel.cor <- function(x, y, digits=2, prefix="", cex.cor, ...)
		{
			usr <- par("usr"); on.exit(par(usr))
			par(usr = c(0, 1, 0, 1))
			r <- abs(cor(x, y))
			txt <- format(c(r, 0.123456789), digits=digits)[1]
			txt <- paste(prefix, txt, sep="")
			if(missing(cex.cor)) cex.cor <- 0.8 / strwidth(txt)
			text(0.5, 0.5, txt, cex = cex.cor * r)
		}
		pairs(apply(x, c(1, 3), sum),
			  lower.panel=panel.smooth,
			  upper.panel=panel.cor,
			  pch=5, ...)

	}
	else
	{
## Otherwise, make the parameter plots:
		if (length(intersect(effects, 1:1:nrow(x$ToTest)))!=length(effects))
		{
			cat("Detected an error with the effects included. For a parameter
				-plot, use the following indices:")
			tmp <- paste(" (", 1:length(timetest$NonRateIndices), ") ",
						 rownames(timetest$IndividualTest)
						 [timetest$NonRateIndices], "\n", sep="")
			cat("\nUse the following indices for plotting the effect-wise
				fitted parameters:\n", tmp)
			stop(" ")
		}
		if (multiplot)
		{
			dims=c(1,1)
			nplots=length(effects)
		}
		else
		{
		nplots=(dims[1] * dims[2])
		}
		if (length(effects) > nplots)
		{
			stop("You have included space for ", nplots, " plots, but have
				 requested ",
				 length(effects), " effects to be plotted. Use dims=c(x, y).")
		}
		xaxis <- 1:timetest$Waves
		if (length(effects)==1)
		{
			yaxis <- timetest$IndividualTest[as.vector(c(effects,
									timetest$DummyIndexByEffect[effects, ])), 2]
			dim(yaxis) <- c(1, timetest$Waves)

		}
		else
		{
			yaxis <- timetest$IndividualTest[as.vector(t(cbind(effects,
									timetest$DummyIndexByEffect[effects, ]))), 2]
			yaxis <- matrix(yaxis, nrow=length(effects), ncol=timetest$Waves,
							byrow=TRUE)
		}
		rownames(yaxis) <- rownames(timetest$IndividualTest)[effects]
		colnames(yaxis) <- 1:timetest$Waves
		vals <- yaxis
		basevals <- array(yaxis[, 1], dim=dim(yaxis))
		basevals[, 1] <- 0
		yaxis <- yaxis + basevals
		pvals <- timetest$IndividualTest[c(effects, as.vector(
									t(timetest$DummyIndexByEffect[effects, ]))), 3]
		dummysd <- abs(c(vals / qnorm(1 - pvals / 2)))
		dummysd[effects] <- timetest$OriginalThetaStderr[effects]
		dim(dummysd) <- c(length(effects), timetest$Waves)
		dim(pvals) <- dim(dummysd)
		dim(vals) <- dim(dummysd)
		rownames(dummysd) <- rownames(timetest$IndividualTest)[effects]
		colnames(dummysd) <- 1:timetest$Waves
##Function to print a panel:
		makeplot <- function (i)
		{
			ymin=min(yaxis[i, ] - scale * abs(yaxis[i, ]))
			ymax=max(yaxis[i, ] + scale * abs(yaxis[i, ]))
			xyplot(yaxis[i, ] ~ xaxis,
				   type = "p", main = rownames(timetest$EffectTest)[effects[i]] ,
				   sub=paste("p=", timetest$EffectTest[effects[i]]), bty="n",
				   xlab="Wave", ylab="Parameter Value", auto.key=TRUE,
				   ylim=c(ymin, ymax), xlim=c(0, length(xaxis) + 1),
				   panel=function(x, y){
                       for (j in 1:length(x))
                       {
                           if ( c(FALSE, timetest$ToTest[effects[i], ])[j] )
                           {
                               tmp="red"

                           }
                           else
                           {
                               tmp="gray"
                           }
                           l <- yaxis[i, j] - abs(qnorm(plevels[1]  /  2,
                                                        sd=dummysd[i, j]))
                           u <- yaxis[i, j] + abs(qnorm(plevels[1] / 2,
                                                        sd=dummysd[i, j]))
                           panel.xyplot(c(x[j], x[j]), c(l, u), reference=TRUE,
                                        col=tmp, alpha=.50,
                                        type="l", lend=1, lwd=10)
                           panel.xyplot(c(x[j], x[j]), c(l, u), reference=TRUE,
                                        col=tmp, alpha=.75,
                                        type="p", pch=45, cex=3)
                           l <- yaxis[i, j] - abs(qnorm(plevels[2] / 2,
                                                        sd=dummysd[i, j]))
                           u <- yaxis[i, j] + abs(qnorm(plevels[2] / 2,
                                                        sd=dummysd[i, j]))
                           panel.xyplot(c(x[j], x[j]), c(l, u), reference=TRUE,
                                        col=tmp, alpha=.50,
                                        type="l", lend=1, lwd=10)
                           panel.xyplot(c(x[j], x[j]), c(l, u), reference=TRUE,
                                        col=tmp, alpha=.75,
                                        type="p", pch=45, cex=3)
                           l <- yaxis[i, j] - abs(qnorm(plevels[3] / 2,
                                                        sd=dummysd[i, j]))
                           u <- yaxis[i, j] + abs(qnorm(plevels[3] / 2,
                                                        sd=dummysd[i, j]))
                           panel.xyplot(c(x[j], x[j]), c(l, u), reference=TRUE,
                                        col=tmp, alpha=.25,
                                        type="l", lend=1, lwd=10)
                       }
                       panel.xyplot(x, y, type="s", reference=TRUE,
                                    col="black", alpha=.75, pch=2)
                       panel.xyplot(x, y, type="p", reference=TRUE, pch=20,
                                    col=1)
                       panel.abline(a=timetest$IndividualTest[effects[i], 1],
                                    reference=TRUE,
                                    col="black", lwd=2, alpha=.75)

				   }, ...)
		}

		if (length(effects) > 1 & !multiplot)
		{
            print(makeplot(1), newpage=TRUE, more=TRUE, split=c(1, 1, dims[1],
                                                        dims[2]))
			if(dims[1] > dims[2])
			{
				col=1
				row=2

			}
			else
			{
				row=1
				col=2
			}
			if (length(effects) > 2)
			{
				for (i in 2:(length(effects)-1))
				{
					print(makeplot(i), more=TRUE, split=c(row, col, dims[1],
														  dims[2]))
					col <- col + 1
					if(col > dims[2])
					{
						col=1
						row <- row + 1
					}
				}
			}
			print(makeplot(length(effects)), split=c(row, col, dims[1],
                                             dims[2]))
		}

		else if (length(effects) > 1 & multiplot)
		{
			for (i in 1:(length(effects)))
			{
				dev.new()
				print(makeplot(i), newpage=TRUE, more=FALSE,
                      split=c(1, 1, 1, 1))
			}
		}
		else
		{
			print(makeplot(1), newpage=TRUE, more=FALSE, split=c(1, 1, 1, 1))
		}
	}
}
##@sienaTimeFix siena07 Adds time dummy terms to the effects object
sienaTimeFix <- function(effects, data)
{
    if (inherits(data, "sienaGroup"))
    {
        warning("Time dummy not implemented for multi-group projects")
        effects$timeDummy <- ","
    }
    else
    {
        observations <- data$observations - 1
        if (observations < 2 && any(effects$timeDummy != ","))
        {
            warning("Time dummies not relevant with only 2 periods")
            effects$timeDummy <- ","
        }
    }
    use <- effects$name == effects$name[1]
    if (any(effects$timeDummy[!use] != ","))
    {
        warning("Time dummy only implemented for first dependent variable")
        effects$timeDummy[!use] <- ","
    }
    if (length(unique(effects$groupName)) > 1)
    {
        warning("Time dummy not implemented for multi-group projects")
        effects$timeDummy <- ","
    }
 # Josh tested these covariate effects, they work as-is for sienaTimeFix.
 #   covar <- effects$interaction1 != ""
 #   if (any(effects$timeDummy[covar] != ","))
 #   {
 #       warning("Time dummy not implemented for covariate effects")
 #       effects$timeDummy[covar] <- ","
 #   }
    implemented <- (effects$type == "eval" | effects$shortName == "RateX")
	if (any(effects$timeDummy[!implemented] !=','))
	{
		warning("Time dummy effects are only implemented",
                " for one mode network effects of type eval or for RateX.")
        effects$timeDummy[!implemented] <- ","
	}
	if (all(effects$timeDummy==',') )
	{
##		No time dummy interactions to add, so kick the inputs back.
		return(list (effects=effects, data=data))
	}
	else
	{
## 	One mode, eval effects, or RateX effects:
		alreadyDummied <- grep("isDummy", effects$timeDummy)
		effects$timeDummy[effects$timeDummy=="all"]  <-
            paste(2:(data$observations-1), collapse = ",")
		if (length(alreadyDummied)  >  0)
		{
## Just remove those effects that have already been dummied so as to
## not messy things up. The assumption is that the user will retain
## all of the previous dummied effects within the column.
			effects <- effects[-alreadyDummied, ]
		}
		dummiedEffects <- effects$effectNumber[effects$timeDummy != ',' & (effects$type=='eval' | effects$shortName=='RateX')]
		covToAdd <- NULL
		rateCovToAdd <- NULL
		dummyCombos <- list()
		ctr=1
## This might need to be changed for sienaGroup:
		nact=dim(data$depvars[[1]])[1]
		nper=dim(data$depvars[[1]])[3]
		for (i in dummiedEffects)
		{
## Get the time periods that we want dummied for effect i:
			tmp <- toString(effects$timeDummy[effects$effectNumber == i])
			tmp <- strsplit(tmp, split=",", fixed=TRUE)[[1]]
			if (length(which(!tmp == '')) > 0)
			{
				tmp=as.numeric(tmp)
				tmp=tmp[tmp > 1 & tmp  <  nper]

			}
			else
			{
				next
			}
			if (length(which(!is.numeric(tmp))) > 0)
			{
				stop("Invalid input for time dummy column of effects object:", tmp)
			}
			if (length(tmp) > 0)
			{
				if (effects$type[effects$effectNumber==i]=='eval') {
					dummyCombos[[ctr]]=list(effectNumber=i, periods=tmp)
					ctr=ctr + 1
					covToAdd <- unique(c(covToAdd, tmp))
				} else if (effects$shortName[effects$effectNumber==i]=='RateX') {
					## RateX effect, has to be dealt with differently. Just add them now:
					for (p in tmp) {
						dname <- paste(effects$interaction1[effects$effectNumber==i],
								"Dummy",p,sep="")
						base <- matrix(0,nact,nper-1)
						## Figure out the base values:
						dvind <- which(names(data$cCovars) ==
							effects$interaction1[effects$effectNumber==i])
						if ( length(dvind) == 0) {
						## It is a varCovar, not a coCovar 
							dvind <- which(names(data$vCovars) ==
											effects$interaction1[effects$effectNumber==i])
							if (length(dvind)==0) {
								stop("Having trouble finding the covariate for your rate effect. Please
									 contact the developers.")
							}
							base[,p] <- data$vCovars[[dvind]][,p]
						} else {
							## Stick them into the right time spot
							base[,p] <- data$cCovars[[dvind]]
							## make a new varCovar:
						}
						base <- varCovar(base)
						base <- addAttributes.varCovar(base, name=dname)
						data$vCovars[[length(data$vCovars)+1]] <- base
						names(data$vCovars)[length(data$vCovars)] <- dname
						## Now add the rate term:
						tmprow <- allEffects[allEffects$functionName==
										'Amount of change x xxxxxx' & allEffects$type=='rate'
										& allEffects$effectGroup=='covarNonSymmetricRate', ]
						tmprow$name <- effects$name[effects$shortName=='RateX' &
										effects$type=='rate'][1]
						tmprow$effectFn <- 'NULL'
						tmprow$statisticFn <- 'NULL'
						tmprow$netType <- 'oneMode'
						tmprow$groupName <- 'Group1'
						tmprow$group <- 1
						tmprow$fix <- FALSE
						tmprow$include <- TRUE
						tmprow$effectNumber <- max(effects$effectNumber) + 1
						tmprow <- tmprow[, colnames(effects)]
						tmprow$effectName <- gsub('xxxxxx', dname, tmprow$effectName)
						tmprow$functionName <- gsub('xxxxxx', dname, tmprow$functionName)
						tmprow$interaction1 <- dname
						tmprow$timeDummy <- paste('isDummy', p, i, sep=',')
						rownames(tmprow) <- dname
						effects <- rbind(effects, tmprow)
					}
				}
			}
		}
## Add the required covariate effects to the effect objects
		ctr <- length(data$vCovars) + 1
		for (i in covToAdd)
		{
			dname <- paste("Dummy", i, sep='')
			tmp <- array(0, c(nact, nper-1))
			tmp[, i]=1
			tmp <- varCovar(tmp)
			tmp <- addAttributes.varCovar(tmp, name=dname)
			data$vCovars[[ctr]] <- tmp
			names(data$vCovars)[ctr] <- dname
			ctr <- ctr + 1
			tmprow <- allEffects[allEffects$functionName==
			'Sum of outdegrees x xxxxxx' & allEffects$type=='eval'
			& allEffects$effectGroup=='covarNonSymmetricObjective', ]
			tmprow$name <- effects$name[effects$shortName=='density' &
			effects$type=='eval'][1]
			tmprow$effectFn <- 'NULL'
			tmprow$statisticFn <- 'NULL'
			tmprow$netType <- 'oneMode'
			tmprow$groupName <- 'Group1'
			tmprow$group <- 1
			tmprow$fix <- TRUE
			tmprow$include <- TRUE
			tmprow$effectNumber <- max(effects$effectNumber) + 1
			tmprow <- tmprow[, colnames(effects)]
			tmprow$effectName <- gsub('xxxxxx', dname, tmprow$effectName)
			tmprow$functionName <- gsub('xxxxxx', dname, tmprow$functionName)
			tmprow$interaction1 <- dname
			tmprow$timeDummy <- paste('isDummy', i,
									  effects$effectNumber[effects$shortName=='density' &
									  effects$type=='eval'], sep=',')
			rownames(tmprow) <- dname
			effects <- rbind(effects, tmprow)
		}
		for (i in seq(along=dummyCombos))
		{
			baseNum=dummyCombos[[i]]$effectNumber
			for (j in seq(along=dummyCombos[[i]]$periods))
			{
				dname <- paste("Dummy", dummyCombos[[i]]$periods[j], sep="")
				dummyNum <- effects$effectNumber[rownames(effects)==dname]
				if (effects$shortName[baseNum] != 'density')
				{
## Make a user specified interaction
## for the time dummy interacted effect
					tmprow <- allEffects[allEffects$shortName=='unspInt'
					& allEffects$type=='eval'
					& allEffects$effectGroup=='unspecifiedNetInteraction', ]
					tmprow$name <- effects$name[effects$effectNumber==baseNum]
					tmprow$effectFn <- 'NULL'
					tmprow$statisticFn <- 'NULL'
					tmprow$netType <- 'oneMode'
					tmprow$groupName <- 'Group1'
					tmprow$group <- 1
					tmprow$fix <- FALSE
					tmprow$include <- TRUE
					tmprow$effectNumber <- max(effects$effectNumber) + 1
					tmprow <- tmprow[, colnames(effects)]
					tmprow$effectName <- 'unspecified interaction effect'
					tmprow$functionName <- 'unspecified interaction statistic'
					rownames(tmprow) <- paste(dname, baseNum, sep='.')
					tmprow$effect1 <- baseNum
					tmprow$effect2 <- dummyNum
					tmprow$timeDummy <- paste('isDummy',
											  dummyCombos[[i]]$periods[j], baseNum, sep=',')
					effects <- rbind(effects, tmprow)

				}
				else
				{
					effects$fix[effects$effectNumber==dummyNum] <- FALSE
				}
			}
		}
		list(effects=effects, data=data)
	}
}
##@includeTimeDummy DataCreate
includeTimeDummy <- function(myeff, ..., timeDummy="all", name=myeff$name[1],
		type="eval", interaction1="", interaction2="", include=TRUE,
		character=FALSE)
{

	if (character)
	{
		dots <- sapply(list(...), function(x)x)
	}
	else
	{
		dots <- substitute(list(...))[-1] ##first entry is the word 'list'
	}
	if (length(dots) == 0)
	{
		stop("need some effect short names")
	}
	if (!character)
	{
		effectNames <- sapply(dots, function(x)deparse(x))
	}
	else
	{
		effectNames <- dots
	}
	use <- myeff$shortName %in% effectNames &
			myeff$type==type &
			myeff$name==name &
			myeff$interaction1 == interaction1 &
			myeff$interaction2 == interaction2
	myeff[use, "timeDummy"] <- timeDummy
    myeff[use, "include"] <- include
    print.data.frame(myeff[use,])
	myeff
}
