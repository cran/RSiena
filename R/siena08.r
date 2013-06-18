#/******************************************************************************
# * SIENA: Simulation Investigation for Empirical Network Analysis
# *
# * Web: http://www.stats.ox.ac.uk/~snijders/siena
# *
# * File: siena08.r
# *
# * Description: This module contains the code for the meta analysis of a
# * collection of Siena fits.
# *****************************************************************************/
##@siena08 siena08
siena08 <- function(..., projname="sienaMeta", bound=5, alpha=0.05, maxit=20)
{
    dots <- as.list(substitute(list(...)))[-1] ##first entry is the word 'list'
    if (length(dots) == 0)
    {
        stop('need some sienafits')
    }
    nm <- names(dots)
    if (is.null(nm))
    {
        fixup <- seq(along=dots)
    }
    else
    {
        fixup <- nm == ''
    }
    dep <- sapply(dots[fixup], function(x) deparse(x)[1])
    if (is.null(nm))
    {
        nm <- dep
    }
    else if (length(dep) > 0)
    {
        nm[fixup] <- dep
    }
    dots <- list(...)
    names(dots) <- nm
    if (any(duplicated(nm)))
    {
        stop('names must be unique')
    }
    ex <- dots
    projnames <- sapply(ex, function(x)x$x$projname)
    if (any(duplicated(projnames)))
    {
        projnames <- paste(projnames, "(", names(ex), ")", sep="")
    }

    ##first make a data frame
    mydf <- do.call(rbind, lapply(1:length(ex), function(i, y)
                              {
                                  x <- y[[i]]
                                  n <- length(x$theta)
                                  scoretests <- rep(NA, length(x$theta))
                                  if (!is.null(x$testresulto))
                                  {
                                      scoretests[x$test] <-
                                          x$testresulto[!is.na(x$testresulto)]
                                  }
                                  data.frame(objname=rep(names(y)[i], n),
                                             projname=rep(projnames[i], n),
                                             theta=x$theta,
                                             effects=
                                             paste(format(x$requestedEffects$type,
                                                          width=5),
                                                   x$requestedEffects$effectName,
                                                   sep=": "),
                                             tconv=x$tconv,
                                             version=rep(x$version, n),
                                             revision=rep(x$revision, n),
                                             scoretests=scoretests,
                                             ## add anything more needed
                                             ## for the report
                                             se=sqrt(diag(x$covtheta)))

                              }, y=ex
                                  )
                    )
    ## make sure the effects are in the right order.
    mydf$effects <- factor(mydf$effects, levels=unique(mydf$effects))
    ##count the score tests

    ## produce meta analysis object
    dometa <- function(x)
    {
        ##tidy up the data frame x first so theta and sj2 match in length
        ##not actually necessary for iwlsm, as model.frame deals with it.
        ## but easier for the remainder...
        x1 <- x[!is.na(x$theta) & !is.na(x$se) & x$se < bound,]
        if (any(x1$theta != 0))
        {
            if (sum((x1$se < bound)) >= 3)
            {
                suppressWarnings(check.correl <- cor.test(x1$theta, x1$se,
                                                          method="spearman"))
                ## warnings will be given in case of ties, not important here
            }
            else
            {
                check.correl <- data.frame(estimate=NA, p.value=NA,
                                           method="no correlation test")
            }
            regfit <- iwlsm(theta ~ 1, psi=psi.iwlsm, data=x1,
                            ses=x1$se^2, maxit=maxit)
            regfit$terms <- NA
            regfit$model <- NULL
            regfit$psi <- NULL
            ## symbols ttilde, Qstat, Tsq as in Snijders & Baerveldt (2003),
            ##(18), (17), (15)
            Tsq <- sum((x1$theta / x1$se)^2)
            regsummary <- summary(regfit)
            tratio <- regsummary$coef[1, 3]
            ttilde <- sum(x1$theta / x1$se^2) / sqrt(sum(1 / x1$se^2))
            Qstat <- Tsq - ttilde^2
            cjplus <- -2 * sum(pnorm(x1$theta / x1$se, lower.tail=FALSE,
                                     log.p=TRUE))
            cjminus <- -2 * sum(pnorm(x1$theta / x1$se, log.p=TRUE))
            cjplusp <- 1 - pchisq(cjplus, 2 * nrow(x1))
            cjminusp <- 1 - pchisq(cjminus, 2 * nrow(x1))
            ## ML estimates and confidence intervals
            maxxlik <- maxlik(x1$theta, x1$se)
            cmu  <- confint.mu(x1$theta, x1$se, alpha)
            csig <- confint.sig(x1$theta, x1$se, alpha)
            ret1 <- list(cor.est=check.correl$estimate,
                         cor.pval=check.correl$p.value,
                         cor.meth=check.correl$method,
                         regfit=regfit, regsummary=regsummary,
                         Tsq=Tsq, pTsq=1 - pchisq(Tsq, nrow(x1) - 1),
                         tratio=tratio,
                         ptratio=2 * pnorm(abs(tratio), lower.tail=FALSE),
                         Qstat=Qstat,
                         pttilde=1 - pchisq(Qstat, nrow(x1) - 1),
                         cjplus=cjplus, cjminus=cjminus,
                         cjplusp=cjplusp, cjminusp=cjminusp, n1=nrow(x1),
                         mu.ml=maxxlik$mu, sigma.ml=maxxlik$sigma,
                         mu.ml.se=maxxlik$se.mu,
                         mu.confint=cmu, sigma.confint=csig)
        }
        else
        {
            ret1 <- NULL
            ret1$n1 <- 0
        }
        if (any(!is.na(x$scoretests)))
        {
            n <- sum(!is.na(x$scoretests))
            cjplus <- -2 * sum(pnorm(x$scoretests, lower.tail=FALSE,
                                     log.p=TRUE), na.rm=TRUE)
            cjminus <- -2 * sum(pnorm(x$scoretests, log.p=TRUE),
                                na.rm=TRUE)
            cjplusp <- 1 - pchisq(cjplus, 2 * n)
            cjminusp <- 1 - pchisq(cjminus, 2 * n)
            ret1$scoreplus <- cjplus
            ret1$scoreminus <- cjminus
            ret1$scoreplusp <- cjplusp
            ret1$scoreminusp <- cjminusp
            ret1$ns <- n
        }
        else
        {
            ret1$ns <- 0
        }
        ret1
    }

    meta <- by(mydf, mydf$effects, dometa)
    meta$thetadf <- mydf
    class(meta) <- "sienaMeta"
    meta$projname <- projname
    meta$bound <- bound
    ## count the score tests
    meta$scores <- by(mydf, mydf$effects, function(x)
                      any(!is.na(x$scoretests)))
    meta
}


## methods

##@print.sienaMeta Methods
print.sienaMeta <- function(x, file=FALSE, ...)
{
    exitfn <- function()
    {
        if (file)
        {
            ## close the report file
            Report(closefiles=TRUE)
        }

    }
    on.exit(exitfn())
    projname <- x$projname
    if (file)
    {
        Report(openfiles=TRUE, type="a", projname=projname) # initialise a file
    }
    else
    {
        Report(openfiles=TRUE, type="n") #initialise with no file
    }
    ## projnames <- unique(x$thetadf$projname)
    ## nProjects <- length(projnames)
    effects <- unique(x$thetadf$effects)
    ## nEffects <- length(effects)
    ## results

    ## estimates

    Report(c("\nUpper bound used for standard error is",
             format(round(x$bound, 4), width=9, nsmall=2), ".\n"), sep="", outf)

    dashes <- paste(rep("-", 80), collapse="")
    x$thetadf$excl <- ifelse(x$thetadf$se > x$bound,
                             " EXCLUDED from meta-analysis", "")
    by(x$thetadf, x$thetadf$effects, function(x, y)
   {
       i <- match(x$effects[1], effects)
       y <- y[[effects[i]]]
       Report(c("\n", dashes, "\nParameter ", i, ": ",
                as.character(x$effects[1]), "\n", dashes, "\n"),
              sep="", outf)
       tmp <- paste("Data set ", 1:nrow(x), ", ", format(x$projname),
                    " :  Estimate ",
                    format(round(x$theta, 4), width=12),
                    " (standard error ",
                    format(round(x$se, 4), nsmall=4,
                           width=12), ")", x$excl, "\n", sep="")
       Report(c(tmp, "\n"), sep="", outf)
       Report(c(" ", y$n1, " datasets used.\n\n"), sep="", outf)
       if (y$n1 > 0)
       {
           Report("Test that estimates and standard errors are uncorrelated",
                  outf)
           if (is.na(y$cor.est))
           {
               Report("\ncannot be performed.\n\n", outf)
           }
           else
           {
               Report(c(": \n", y$cor.meth, " =", format(round(y$cor.est, 4),
                                                         width=9),
                        ", two-sided ",reportp(y$cor.pval,3), "\n\n"), sep="",
                      outf)
           }
           Report(c("Estimates and test based on IWLS modification of",
                    "Snijders & Baerveldt (2003)\n"),
                  outf)
           Report(c("---------------------------------------------------",
                    "-------------------------\n"), sep="",
                  outf)
           Report(c("Test that all parameters are 0 : \n"), outf)
           Report(c("chi-squared =", format(round(y$Tsq, 4), width=9),
                    ", d.f. = ", y$n1, ", ",
                    reportp(y$pTsq, 3), "\n\n"), sep="", outf)
           Report(c("Estimated mean parameter",
                    format(round(y$regsummary$coefficients[1, 1], 4), width=9),
                    " (s.e.", format(round(y$regsummary$coefficients[1, 2], 4),
                                     width=9), "), two-sided ",
                    reportp(2 * pt(-abs(y$regsummary$coefficients[1, 3]),
                                   y$n1 - 1), 3), "\n\n"), sep="", outf)
           Report(c("Estimated standard deviation",
                    format(round(y$regsummary$stddev, 4), width=9)), outf)
           Report("\nTest that variance of parameter is 0 :\n", outf)
           Report(c("Chi-squared = ", format(round(y$Qstat, 4), width=9),
                    " (d.f. = ", y$n1-1, "), ", reportp(y$pttilde, 3),
                    "\n\n"), sep="", outf)
           Report(c("Estimates and confidence intervals under normality",
                    "assumptions\n"),
                  outf)
           Report(c("-------------------------------------------------------",
                    "-------\n"), outf)
           Report(c("Estimated mean parameter",
                    format(round(y$mu.ml, 4), width=9),
                    " (s.e.",format(round(y$mu.ml.se, 4), width=9),
                    "), two-sided ",
                    reportp(2 * pt(-abs(y$mu.ml/y$mu.ml.se),
                                   y$n1 - 1), 3), "\n"), sep="", outf)
           Report(c(format(round(y$mu.confint[3], 2), width=4),
                    "level confidence interval [",
                    format(round(y$mu.confint[1], 4), width=7),
                    ",",
                    format(round(y$mu.confint[2], 4), width=7), "]\n"), outf)
           Report(c("Estimated standard deviation",
                    ifelse((y$sigma.ml > 0.0001) | (y$sigma.ml < 0.0000001),
                           format(round(y$sigma.ml, 4), width=9), " < 0.0001"),
                     "\n"), outf)
           Report(c(format(round(y$sigma.confint[3], 2), width=4),
                    "level confidence interval [",
                    format(round(y$sigma.confint[1], 4), width=7),
                    ",",
                    format(round(y$sigma.confint[2], 4), width=7), "]\n\n"), outf)
           Report("Fisher's combination of one-sided tests\n", outf)
           Report("----------------------------------------\n", outf)
           Report("Combination of right one-sided p-values:\n", outf)
           Report(c("Chi-squared = ", format(round(y$cjplus, 4), width=9),
                    " (d.f. = ", 2 * y$n1, "), ", reportp(y$cjplusp, 3),
                    "\n"), sep="", outf)
           Report("Combination of left one-sided p-values:\n", outf)
           Report(c("Chi-squared = ", format(round(y$cjminus, 4), width=9),
                    " (d.f. = ", 2 * y$n1, "), ", reportp(y$cjminusp, 3),
                    "\n"), sep="", outf)
       }
       else
       {
           Report(c("There were no data sets satisfying the bounds for",
                    "this parameter.\n No combined output is given.\n"), outf)
       }
   }, y=x)
    ##score tests
    if (any(x$scores))
    {
        Report(c("\n\n", paste(rep("-", 65), collapse=""),
                 "\nScore tests:\nFisher combination\n",
                 paste(rep("-", 65), collapse=""), "\n"), sep="", outf)

        by(x$thetadf, x$thetadf$effects, function(x, y)
       {
           i <- match(x$effects[1], effects)
           y <- y[[effects[i]]]
           if (y$ns > 0)
           {
               Report(c("\n", "(", i, ")   ",
                        as.character(x$effects[1]), "\n"), sep="", outf)
               tmp <- paste("Data set ", 1:nrow(x), ", ", format(x$projname),
                            " : z = ", ifelse(is.na(x$scoretests), "NA",
                                              format(round(x$scoretests, 4),
                                                     width=12)),
                            "\n", sep="")
               Report(c(tmp, "\n"), sep="", outf)
               Report("Combination of right one-sided p-values:\n", outf)
               Report(c("Chi-squared = ", format(round(y$scoreplus, 4),
                                                 width=9),
                        " (d.f. = ", 2 * y$ns, "), ",
                        reportp(y$scoreplusp, 3), "\n"), sep="", outf)
               Report("Combination of left one-sided p-values:\n", outf)
               Report(c("Chi-squared = ",
                        format(round(y$scoreminus, 4), width=9),
                        " (d.f. = ", 2 * y$ns, "), ",
                        reportp(y$scoreminusp, 3), "\n"), sep="", outf)
           }
       }, y=x)
    }
    invisible(x)
}

##@reportp Miscellaneous
reportp <- function(p, ndec)
{
    thr <- exp(-ndec * log(10.0));
    if (is.na(p))
    {
        "p undefined"
    }
    else if (p > thr)
    {
        c("p = ", format(round(p, ndec), width=ndec+2, scientific=FALSE,
                         nsmall=ndec))
    }
    else
    {
        c("p < ", format(round(thr, ndec), width=ndec+2, scientific=FALSE,
                         nsmall=ndec))
    }

}

##@plot.sienaMeta Methods
plot.sienaMeta <- function(x, ..., layout = c(2,2))
{
    library(lattice)
    tmp <- xyplot(theta ~ se|effects,
                  data=x$thetadf[is.na(x$thetadf$scoretests),],
                  ylab="estimates",
                  xlab="standard errors", layout=layout,
                  panel=function(x, y)
              {
                  panel.xyplot(x, y)
                  panel.abline(0, qnorm(0.025))
                  panel.abline(0, qnorm(0.975))
              },
                  prepanel=function(x,y)
              {   list(xlim=c(min(0,min(x)),max(0,max(x))),
                       ylim=c(min(0,min(y)),max(0,max(y))))
              },
                scales="free")
    tmp[!sapply(tmp$y.limits, function(x)all(is.na(x)))]
}

##@summary.sienaMeta Methods
summary.sienaMeta <- function(object, file=FALSE, extra=TRUE, ...)
{
    object$file <- file
    object$extra <- extra
    class(object) <- c("summary.sienaMeta", class(object))
    object
}

##@print.summary.sienaMeta Methods
print.summary.sienaMeta <- function(x, file=FALSE, extra=TRUE, ...)
{
    exitfn <- function()
    {
        if (file)
        {
            ## close the report file
            Report(closefiles=TRUE)
        }

    }
    on.exit(exitfn())
    if (!is.null(x$file))
    {
        file <- x$file
    }
    if (!is.null(x$extra))
    {
        extra <- x$extra
    }
    projname <- x$projname
    if (file)
    {
        Report(openfiles=TRUE, type="a", projname=projname) # initialise a file
    }
    else
    {
        Report(openfiles=TRUE, type="n") #initialise with no file
    }
    if (file)
    {
        ## do some sums for the heading

        namelen <- nchar(projname) + 4
        ##linelen <- 80
        astlen <- min(namelen + 10, 80)
        nBlanks <- if (astlen < 80)
        {
            (80 - astlen) %/% 2
        }
        else
        {
            0
        }
        nBlanks2 <- if (namelen < 80)
        {
            (80 - namelen) %/% 2
        }
        else
        {
            0
        }
        Report(c(rep(" ", nBlanks), rep("*", astlen), "\n"), sep="", outf)
        Report(c(rep(" ", nBlanks2), projname, ".out\n"), sep="", outf)
        Report(c(rep(" ", nBlanks), rep("*", astlen), "\n"), sep="", outf)
        Report(c("Filename is ", projname, ".out.\n\n"), sep="", outf)
        Report(c("This file contains primary output for SIENA project <<",
                 projname, ">>.\n\n"), sep="", outf)
        Report(c("Date and time:", format(Sys.time(),
                                          "%d/%m/%Y %X"), "\n\n"), outf)
        packageValues <- packageDescription(pkgname,
                                            fields=c("Version", "Date"))
        rforgeRevision <-  packageDescription(pkgname,
                                              fields=
                                              "Repository/R-Forge/Revision")
        if (is.na(rforgeRevision))
        {
            revision <- ""
        }
        else
        {
            revision <- paste(" R-forge revision: ", rforgeRevision,
                              " ", sep="")
        }
        Report(c("SIENA version ", packageValues[[1]], " (",
                 format(as.Date(packageValues[[2]]), "%d %m %Y"), ")",
                 revision, "\n\n"), sep="", outf)
    }
    Report(c("================================= SIENA08 ",
             "================================================\n",
             "Multilevel use of Siena algorithms according to ",
             "Snijders & Baerveldt (2003) with extension\n",
             "=================================================",
             "=========================================\n\n"), sep="", outf)
    projnames <- unique(x$thetadf$projname)
    nProjects <- length(projnames)
    effects <- unique(x$thetadf$effects)
    nEffects <- length(effects)
    Report(c("Number of projects in the list is " ,
             length(projnames), ".\n"), sep="", outf)
    Report("The names of these projects are :\n", outf)
    tmp <- paste("project", 1:nProjects, ": <", projnames,
                 ">\n", sep="")
    Report(tmp, sep="", outf)
    Report(c("\nOptions for running Siena08:\n",
             "-> Parameters are excluded from the meta-analysis when their ",
             "standard\n", "   error exceeds an upper bound of ",
             round(x$bound,digits=2), ".\n"), sep="", outf)
    if (extra)
    {
        Report("-> Extra output requested\n", outf)
    }
    else
    {
        Report("-> No extra output requested\n", outf)
    }
    ## RSiena version
    Report(c("\nThe RSiena Version of the first fit object is ",
             x$thetadf$version[1], ".\n\n"), sep="", outf)
    ## project names
    by(x$thetadf, x$thetadf$projname, function(x)
   {
       Report(c("Object <", x$projname[1], "> contains estimates of ",
                nrow(x), " parameters.\n"), sep="", outf)
       Report(c("The number of valid score tests found was ",
                sum(is.na(x$scoretests)), ".\n"),
              sep="", outf);
   })
    ##parameters:
    Report(c("\nA total of", nEffects, "parameters in", nProjects,
             "projects :\n"),   outf)
    Report(paste(format(1:length(effects)), ". " , effects, "\n", sep=""),
           sep="", outf)
    Report(c("\nThe projects contain the parameters as follows",
             "(1=present, 0=absent):\n\n"), outf)
    row1 <- c(1:nEffects)
    rows <- do.call(rbind, tapply(x$thetadf$effects,
                                  x$thetadf$projname, function(x)
                              {
                                  as.numeric(effects %in% x)
                              }
                                  )
                    )
    rows <- format(rbind(row1, rows), width=3)
    row2 <- rep(paste(rep("-", nchar(rows[1, 1])), collapse=""), nEffects)
    rows <- rbind(rows[1,], row2, rows[-1, ])
    col1 <- format(rbind("Project",
                         paste(rep("-", 7), collapse=""),
                         cbind(1:nProjects)), justify="centre")
    col2 <- format(rbind( "|",  "--+--",  cbind(rep("|", nProjects))),
                   justify="centre")
    rows <- cbind(col1, col2, rows)
    Report(t(rows), sep = c(rep.int("", ncol(rows) - 1), "\n"), outf)
    by(x$thetadf, x$thetadf$projname, function(x)
   {
       Report(c("\nProject", x$projname[1], "\n"), outf)
       tmp <- paste("par.", 1:nEffects, "estimate",
                    format(round(x$theta, 4), width=12), "(s.e. ",
                    format(round(x$se, 4), width=11), ") (conv_t",
                    format(round(x$tconv, 4), width=9), ")\n")
       Report(format(tmp), sep="", outf)
       Report(c("\nMaximal absolute convergence t-statistic = ",
                format(round(max(abs(x$tconv)), 4), width=9), "\n"), outf)
       ## score tests
   })
    ## results

    ## estimates
    Report(c("\n\n", paste(rep("=", 29), collapse=""),
             "\nResults of the meta-analysis:\n",
             paste(rep("=", 29), collapse=""), "\n"), sep="", outf)

    Report(c("\nUpper bound used for standard error is",
             format(round(x$bound, 4), width=9, nsmall=2), ".\n"), sep="", outf)

    dashes <- paste(rep("-", 80), collapse="")
    x$thetadf$excl <- ifelse(x$thetadf$se > x$bound,
                             " EXCLUDED from meta-analysis", "")
    by(x$thetadf, x$thetadf$effects, function(x, y)
   {
       i <- match(x$effects[1], effects)
       y <- y[[effects[i]]]
       Report(c("\n", dashes, "\nParameter ", i, ": ",
                as.character(x$effects[1]), "\n", dashes, "\n"),
              sep="", outf)
       tmp <- paste("Data set ", 1:nrow(x), ", ", format(x$projname),
                    " :  Estimate ",
                    format(round(x$theta, 4), width=12),
                    " (standard error ",
                    format(round(x$se, 2), nsmall=2,
                           width=11), ")", x$excl, "\n", sep="")
       Report(c(tmp, "\n"), sep="", outf)
       Report(c(" ", y$n1, " datasets used.\n\n"), sep="", outf)
       if (y$n1 > 0)
       {
           if (extra)
           {
               Report(c("IWLS modification of Snijders-Baerveldt (2003) method ",
                        "of combining estimates"), outf)
               Report(c("\n--------------------------------------------",
                        "---------------------------------\n"), sep="", outf)
               Report(c("This method assumes that true parameters and",
                        " standard errors are uncorrelated.\n",
                        "This can be checked by the plot method ",
                        "and the test below.\n\n"), sep="", outf)
           }
           Report("Test that estimates and standard errors are uncorrelated",
                  outf)
           if (is.na(y$cor.est))
           {
               Report("\ncannot be performed.\n\n", outf)
           }
           else
           {
               Report(c(": \n", y$cor.meth, " =", format(round(y$cor.est, 4),
                                                            width=9),
                        ", two-sided ",reportp(y$cor.pval,3), "\n\n"),
                      sep="", outf)
           }
           Report(c("Test that all parameters are 0 : \n"), outf)
           Report(c("chi-squared =", format(round(y$Tsq, 4), width=9),
                    ", d.f. = ", y$n1, ", ",
                    reportp(y$pTsq, 3), "\n\n"), sep="", outf)
           Report(c("Estimated mean parameter",
                    format(round(y$regsummary$coefficients[1, 1], 4), width=9),
                    " (s.e.", format(round(y$regsummary$coefficients[1, 2], 4),
                                     width=9), "), two-sided ",
                    reportp(pt(y$regsummary$coefficients[1, 3],
                               y$n1 - 1), 3), "\n"), sep="", outf)
           Report(c("based on IWLS modification of Snijders & Baerveldt (2003). ",
                  "\n\n"), sep="", outf)
           Report(c("Residual standard error",
                    format(round(y$regsummary$stddev, 4), width=9)), outf)
           Report("\nTest that variance of parameter is 0 :\n",outf)
           Report(c("Chi-squared = ", format(round(y$Qstat, 4), width=9),
	            " (d.f. = ", y$n1-1, "), ", reportp(y$pttilde, 3),
	            "\n"), sep="", outf)
           Report(c("based on IWLS modification of Snijders & Baerveldt (2003).",
                  "\n\n"), sep="", outf)
           Report("Fisher's combination of one-sided tests\n", outf)
           Report("----------------------------------------\n", outf)
           Report("Combination of right one-sided p-values:\n", outf)
           Report(c("Chi-squared = ", format(round(y$cjplus, 4), width=9),
                    " (d.f. = ", 2 * y$n1, "), ", reportp(y$cjplusp, 3),
                    "\n"), sep="", outf)
           Report("Combination of left one-sided p-values:\n", outf)
           Report(c("Chi-squared = ", format(round(y$cjminus, 4), width=9),
                    " (d.f. = ", 2 * y$n1, "), ", reportp(y$cjminusp, 3),
                    "\n"), sep="", outf)
       }
       else
       {
           Report(c("There were no data sets satisfying the bounds for",
                    "this parameter.\n No combined output is given.\n"), outf)
       }
   }, y=x)
    ##score tests
    if (any(x$scores))
    {
        Report(c("\n\n", paste(rep("-", 65), collapse=""),
                 "\nScore tests:\nFisher combination\n",
                 paste(rep("-", 65), collapse=""), "\n"), sep="", outf)

        invisible(by(x$thetadf, x$thetadf$effects, function(x, y)
       {
           i <- match(x$effects[1], effects)
           y <- y[[effects[i]]]
           if (y$ns > 0)
           {
               Report(c("\n", "(", i, ")   ",
                        as.character(x$effects[1]), "\n"), sep="", outf)
               tmp <- paste("Data set ", 1:nrow(x), ", ", format(x$projname),
                            " : z = ",
                            ifelse(is.na(x$scoretests), "NA",
                                   format(round(x$scoretests, 4), width=12)),
                            "\n", sep="")
               Report(c(tmp, "\n"), sep="", outf)
               Report("Combination of right one-sided p-values:\n", outf)
               Report(c("Chi-squared = ", format(round(y$scoreplus, 4),
                                                 width=9),
                        " (d.f. = ", 2 * y$ns, "), ",
                        reportp(y$scoreplusp, 3), "\n"), sep="", outf)
               Report("Combination of left one-sided p-values:\n", outf)
               Report(c("Chi-squared = ",
                        format(round(y$scoreminus, 4), width=9),
                        " (d.f. = ", 2 * y$ns, "), ",
                        reportp(y$scoreminusp, 3), "\n"), sep="", outf)
           }
       }, y=x))
    }
}

