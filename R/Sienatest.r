#/******************************************************************************
# * SIENA: Simulation Investigation for Empirical Network Analysis
# *
# * Web: http://www.stats.ox.ac.uk/~snidjers/siena
# *
# * File: Sienatest.r
# *
# * Description: This module contains the function for instability analysis and
# * score tests.
# *
# *****************************************************************************/
##@InstabilityAnalysis siena07 Not currently used
InstabilityAnalysis<- function(z)
{
	##I think this is not correct, because of scaling. cond number of var matrix of X
	## can be obtained via svd(data) (which is stored in z$sf). Square of ratio
	## of smallest to largest singular value.
	Report('Instability Analysis\n')
	constant <- z$diver
	test <- z$test
	covtheta <- z$covtheta
	covZ <- z$msf
	covth <- covtheta[!(test|constant), !(test|constant)]
	covth <- MatrixNorm(covth)
	eigenv <- eigen(covth,symmetric=TRUE)$values
	ma <- max(eigenv)
	mi <- min(eigenv)
	if (mi!=0)
	{
		cond.n <- ma / mi
	}
	Report('Instability analysis\n', lf)
	Report('--------------------\n\n', lf)
	Report('Variance-covariance matrix of parameter estimates', lf)
	##if (global boolean1 )
	## Report(' (without coordinates that are kept constant):\n',lf)
	##else
	Report(c(':\n\nCondition number = ',format(cond.n, width=4, nsmall=4,
				digits=1),
			' \n\n'),sep='',lf)
	Report(c('Eigen Values  ', format(eigenv, width=6, nsmall=6, digits=1)), lf)
	Report('\n\n',lf)
	covZ <- MatrixNorm(covZ)
	eigenvZ <-eigen(covZ, symmetric=TRUE)$values
	ma <- max(eigenvZ)
	mi <- min(eigenvZ)
	if (mi!=0)
	{
		cond.n <- ma / mi
	}
	Report('Variance-covariance matrix of X',lf)
	Report(c(':\n\nCondition number = ', format(cond.n, width=4, nsmall=4,
				digits=1),
			' \n\n'), sep='', lf)
	Report(c('Eigen Values  ', format(eigenvZ, width=6, nsmall=6, digits=1)), lf)
	Report(c('\n\n',date(),'\n'),sep='',lf)
	mysvd <- svd(z$sf)$d
	ma <- max(mysvd)
	mi <- min(mysvd)
	cond.n <- (ma/mi)^2
	Report(c(':\n\nCondition number2 = ',
			format(cond.n, width=4, nsmall=4, digits=1),
			' \n\n'), sep='', lf)
	Report(c('Singular Values  ', format(mysvd, width=6, nsmall=6, digits=1)), lf)
	Report(c('\n\n', date(), '\n'), sep='', lf)
}

##@MatrixNorm siena07 Not currently used. May be incorrect.
MatrixNorm <- function(mat)
{
	tmp <- apply(mat, 2, function(x) x / sqrt(crossprod(x)))
	##or  sweep(mat,2,apply(mat,2,function(x)x/sqrt(crossprod(x))
	tmp
}
##@TestOutput siena07 Print report
TestOutput <- function(z, x)
{
	testn <- sum(z$test)
	## browser()
	if (testn)
	{
		if (x$maxlike)
		{
			Heading(2, outf,'Score test <c>')
		}
		else
		{
			Heading(2, outf, 'Generalised score test <c>')
		}
		Report('Testing the goodness-of-fit of the model restricted by\n', outf)
		j <- 0
		for (k in 1:z$pp)
		{
			if (z$test[k])
			{
				j <- j + 1
				Report(c(" (", j, ")   ",
						format(paste(z$requestedEffects$type[k], ":  ",
								z$requestedEffects$effectName[k],
								sep=""),
							width=50), " = ",
						sprintf("%8.4f", z$theta[k]), "\n"),
					sep = "", outf)
			}
		}
		Report('_________________________________________________\n',outf)
		Report('                ',outf)
		Report('   \n',outf)
		if (testn > 1)
			Report('Joint test:\n-----------\n',outf)
		Report(c('   c = ',sprintf("%8.4f", z$testresOverall),
				'   d.f. = ',j,'   p-value '),sep='',outf)
		pvalue <- 1-pchisq(z$testresOverall,j)
		if (!is.na(pvalue))
		{
			if (pvalue < 0.0001)
			{
				Report('< 0.0001',outf)
			}
			else
			{
				Report(c('= ',sprintf("%8.4f",pvalue)), sep = '', outf)
			}
		}
		else
		{
			Report('  NA  ',outf)
		}
		if (testn==1)
			Report(c('\n   one-sided (normal variate): ',
					sprintf("%8.4f",z$testresulto[1])), sep = '', outf)
		if (testn> 1)
		{
			Report('\n\n',outf)
			for (k in 1:j)
			{
				Report(c('(',k,') tested separately:\n'),sep='',outf)
				Report('-----------------------\n',outf)
				Report(' - two-sided:\n',outf)
				Report(c('  c = ', sprintf("%8.4f", z$testresult[k]),
						'   d.f. = 1  p-value '), sep = '', outf)
				pvalue<- 1-pchisq(z$testresult[k],1)
				if (!is.na(pvalue))
				{
					if (pvalue < 0.0001)
					{
						Report('< 0.0001\n',outf)
					}
					else
					{
						Report(c('= ', sprintf("%8.4f", pvalue), '\n'), sep = '',
							outf)
					}
				}
				else
				{
					Report('  NA  ',outf)
				}
				Report(c(' - one-sided (normal variate): ',
						sprintf("%8.4f", z$testresulto[k])), sep = '', outf)
				if (k < j)
				{
					Report('\n\n',outf)
				}
			}
		}
		Report('    \n_________________________________________________\n\n',
			outf)
		Report('One-step estimates: \n\n', outf)
		for (i in 1 : z$pp)
		{
			onestepest <- z$oneStep[i] + z$theta[i]
			Report(c(format(paste(z$requestedEffects$type[i], ':  ',
							z$requestedEffects$effectName[i], sep = ''),
						width=50),
					sprintf("%8.4f", onestepest), '\n'), sep = '', outf)
		}
		Report('\n',outf)
	}
}

##@ScoreTest siena07 Do score tests
ScoreTest<- function(pp, dfra, msf, fra, test, redundant, maxlike)
{
	testresult <- rep(NA, pp) ##for chisq per parameter
	testresulto <- rep(NA, pp) ##for one-sided tests per parameter
	##first the general one
	ans <- EvaluateTestStatistic(maxlike, test, redundant, dfra, msf, fra)
	testresOverall <- ans$cvalue
	covMatrix <- ans$covMatrix
	if (sum(test) == 1)
	{
		testresulto[1] <- ans$oneSided
	}
	else
	{
		## single df tests
		use <- !test
		k <- 0
		for (i in 1:pp)
		{
			if (test[i])
			{
				k <- k + 1
				use[i] <- TRUE
				ans <- EvaluateTestStatistic(maxlike, test[use], redundant[use],
					dfra[use, use], msf[use, use], fra[use])
				testresult[k] <- ans$cvalue
				testresulto[k] <- ans$oneSided
				use[i] <- FALSE
			}
		}
	}
	##onestep estimator
	if (maxlike)
		dfra2 <- dfra + msf
	else
		dfra2 <- dfra
	if (inherits(try(dinv2 <- solve(dfra2), silent=TRUE), "try-error"))
	{
		Report("Error message for inversion to get onestep estimator: \n", cf)
		dinv2 <- dfra2
		dinv2[] <- NA
		oneStep <- rep(NA, nrow(dfra2))
	}
	else
	{
		oneStep<- -dinv2 %*% fra
	}
	list(testresult=testresult, testresulto=testresulto,
		testresOverall=testresOverall, covMatrix=covMatrix,
		oneStep=oneStep, dinv2= dinv2, dfra2=dfra2)
}
##@EvaluateTestStatistic siena07 Calculate score test statistics
EvaluateTestStatistic<- function(maxlike, test, redundant, dfra, msf, fra)
{
	##uses local arrays set up in the calling procedure
	d11 <- dfra[!(test|redundant), !(test|redundant), drop=FALSE]
	d22 <- dfra[test, test, drop=FALSE]
	d21 <- dfra[test, !(test|redundant), drop=FALSE]
	d12 <- t(d21)
	sigma11 <- msf[!(test|redundant), !(test|redundant), drop=FALSE]
	sigma22<- msf[test, test,drop=FALSE]
	sigma12 <- msf[!(test|redundant), test, drop=FALSE]
	sigma21<- t(sigma12)
	z1 <- fra[!(test|redundant)]
	z2 <- fra[test]
	if (inherits(try(id11 <- solve(d11), silent=TRUE), "try-error"))
	{
		cat('Score test: Error for inversion of d11 \n')
		oneSided <- NA
		v9 <- d22
		v9[] <- NA
		cvalue <- matrix(NA, 1, 1)
	}
	else
	{
		rg <- d21 %*% id11
		if (!maxlike)
		{
			##orthogonalise deviation vector
			ov <- z2 - rg %*% z1
			##compute var(ov) = sigma22 - (d21 %*% id11) %*% sigma12 -
			##      sigma21 %*% t(id11)%*% t(d21) +
			##      d21%*%id11 %*% sigma11 %*% t(id11) %*% t(d21)
			v2 <- sigma21 - rg %*% sigma11
			v6 <- v2 %*% t(id11) %*% t(d21)
			v9 <- sigma22 -  rg %*% sigma12 - v6
		}
		else
		{
			ov <- -z2
			v9 <- d22 - rg %*% d12
		}
		if (inherits(try(vav <- solve(v9), silent=TRUE), "try-error"))
			## vav is the inverse variance matrix of ov
		{
			cat('Score test: Error for inversion of v9\n')
			vav <- v9
			vav[] <- NA
			cvalue <- NA
			oneSided <- NA
		}
		else
		{
			cvalue <- t(ov) %*% vav %*% ov
			if (cvalue < 0) cvalue <- 0
			if (sum(test) == 1)
			{
				if (vav > 0)
					oneSided <- ov * sqrt(vav)
				else
					oneSided <- 0
				if (!maxlike) oneSided <- - oneSided
				## change the sign for intuition for users
			}
			else
			{
				oneSided <- 0
			}
		}
	}
	list(cvalue=cvalue, oneSided=oneSided, covMatrix=v9)
}


##@scoreTest Calculate score test
score.Test <- function(ans, test=ans$test)
	# use: ans must be a sienaFit object;
	# test must be a boolean vector with length equal to the number of parameters of ans,
	# or a vector of integer numbers between 1 and ans$pp.
{
	if ((is.numeric(test)) || (is.integer(test)))
	{
		if (max(test) > ans$pp)
		{
			stop(paste('The maximum requested coordinate is too high:',
					max(test)))
		}
		test <- (1:ans$pp) %in% test
	}
	if (sum(test) <= 0) stop(paste('Something should be tested, but the total requested is',
			sum(test)))
	if (length(test) != ans$pp) stop('Dimensions of must agree')
	if (any(test & (!ans$fix))) cat('Warning: some tested parameters were not fixed; do you know what you are doing??? \n')	
	fra <- colMeans(ans$sf, na.rm=TRUE)
	redundant <- (ans$fix & (!test))
	teststat <- EvaluateTestStatistic(ans$maxlike, test, redundant, ans$dfra, ans$msf, fra)$cvalue
	df <- sum(test)
	pval <- 1 - pchisq(teststat, df)
	list(chisquare = teststat, df = df, pvalue = pval)
}

##@Wald.RSiena  Calculate Wald test statistics
Wald.RSiena <- function(A, ans)
{
	if (is.vector(A)){A <- matrix(A, nrow=1)}
	if (dim(A)[2] != ans$pp){stop(paste('A must have', ans$pp, 'columns.'))}
	sigma <- ans$covtheta
	if (any(is.na(sigma))) {
		# happens when some coordinates were fixed
		# in the call of siena07 leading to ans;
		# then the non-used part of sigma,
		# which partially consists of NA,
		# is replaced by the identity matrix.
		zero.cols <- apply(A, 2, function(colum){all(colum==0)})
		sigma[zero.cols, ] <- 0
		sigma[, zero.cols] <- 0
		diag(sigma)[zero.cols] <- 1
	}
	th <- A %*% ans$theta
	covmat <- A %*% sigma %*% t(A)
	chisq <- drop(t(th) %*% solve(covmat) %*% th)
	d.f. <- nrow(A)
	pval <- 1 - pchisq(chisq, d.f.)
	list(chisquare = chisq, df = d.f., pvalue = pval)
}

##@Multipar.RSiena  Calculate Wald test statistic for hypothesis that subvector = 0.
Multipar.RSiena <- function(ans, ...)
{
	p <- length(ans$theta)
	k <- length(c(...))
	A <- matrix(0, nrow=k, ncol=p)
	A[cbind(1:k,c(...))] <- 1
	Wald.RSiena(A, ans)
}
