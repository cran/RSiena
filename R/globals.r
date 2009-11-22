##/*****************************************************************************
## * SIENA: Simulation Investigation for Empirical Network Analysis
## *
## * Web: http://www.stats.ox.ac.uk/~snidjers/siena
## *
## * File: globals.r
## *
## * Description: This file contains the code to create and use global objects
## *
## ****************************************************************************/
##@outf Objects/File project .out file
outf <- NULL
##@lf Objects/File suppressed or to console
lf <- NULL
##@bof Objects/File suppressed or to console
bof <- NULL
##@cf Objects/File suppressed or to console
cf <- NULL

##@Reportfun Reporting Part of global mechanism
Reportfun<- function(x, verbose = FALSE)
{
    x <- x
    beverbose <- verbose
    function(txt, dest, fill=FALSE, sep=" ", hdest,
             open=FALSE, close=FALSE,
             type=c("a", "w"),  projname="Siena" , verbose=FALSE)
    {
        if (open)
        {
            type <- match.arg(type)
            beverbose <<- verbose
            if (type =='w')
            {
                x$outf <<- file(paste(projname, ".out", sep=""), open="w")
            }
            else
            {
                x$outf <<- file(paste(projname, ".out", sep=""), open="a")
            }

        }
        else if (close)
        {
            close(x[["outf"]])
        }
        else
        {
            if (missing(dest) && missing(hdest))
            {
                cat(txt, fill = fill, sep = sep)
            }
            else
            {
                if (missing(dest))
                {
                    if (hdest  %in% c("cf", "lf", "bof"))
                    {
                        if (beverbose)
                        {
                            cat(txt, fill=fill, sep=sep)
                        }
                    }
                    else
                    {
                        cat(txt, file = x[[hdest]], fill = fill, sep = sep)
                    }
                }
                else
                {
                    if (deparse(substitute(dest)) %in% c("cf", "lf", "bof"))
                    {
                        if (beverbose)
                        {
                            cat(txt, fill=fill, sep=sep)
                        }
                    }
                    else
                    {
                        cat(txt, file=x[[deparse(substitute(dest))]],
                            fill=fill, sep=sep)
                    }
                }
            }
       }
    }
}

##@Report Globals
Report <- local({verbose <-  NULL;
                 Reportfun(list(outf=outf, lf=lf, cf=cf, bof=bof), verbose)})
##@UserInterrupt Siena07/GlobalFunctions Global (within siena07)
UserInterrupt <- local({A <-  FALSE;function(x){if (!missing(x))A<<-x;A}})
##@EarlyEndPhase2 siena07/GlobalFunctions
EarlyEndPhase2 <- local({A <-  FALSE;function(x){if (!missing(x))A<<-x;A}})
##@UserRestart siena07/GlobalFunctions Global (within siena07)
UserRestart <- local({A <-  FALSE;function(x){if (!missing(x))A<<-x;A}})
##@UserInterruptFlag siena07/GlobalFunctions Global (within siena07)
UserInterruptFlag <- local({A <-  FALSE;function(x){if (!missing(x))A<<-x;A}})
##@EarlyEndPhase2Flag siena07/GlobalFunctions Global (within siena07)
EarlyEndPhase2Flag <- local({A <-  FALSE;function(x){if (!missing(x))A<<-x;A}})
##@UserRestartFlag siena07/GlobalFunctions Global (within siena07)
UserRestartFlag <- local({A <-  FALSE;function(x){if (!missing(x))A<<-x;A}})
##@is.batch siena07/GlobalFunctions Global (within siena07)
is.batch <- local({A <-  FALSE;function(x){if (!missing(x))A<<-x;A}})
##@DONE siena01/GlobalFunctions Used to communicate with siena.exe and sienaScript
DONE <- local({A <-  FALSE;function(x){if (!missing(x))A<<-x;invisible(A)}})
##@FRANstore siena07/GlobalFunctions Used to pass data to other processes
FRANstore <- local({A <-  NULL;function(x){if (!missing(x)) A<<-x;A}})

##@Heading Reporting Global function
Heading<- function(level=1, dest, text, fill=FALSE)
{
    ch <- c("=", "-", " ")[level]
    if (missing(dest))
    {
        Report(c("@", level, "\n", text, "\n"), sep="", fill=fill)
        Report(rep(ch, sum(nchar(text)) + 3), sep="", fill=fill)
        Report("\n\n")
    }
    else
    {
        dest <- deparse(substitute(dest))
        Report(c("@", level, "\n", text, "\n"), hdest=dest, sep="", fill=fill)
        Report(rep(ch, sum(nchar(text))), hdest=dest, sep="", fill=fill)
        if (level < 3)
            Report("\n\n", hdest = dest)
        else
            Report("\n", hdest = dest)
    }
}

##@PrtOutMat Reporting
PrtOutMat<- function(mat, dest)
{
    if (missing(dest))
        Report(format(t(mat)), sep=c(rep.int(" ", ncol(mat) - 1), "\n"))
    else
    {
        Report(format(t(mat)), sep=c(rep.int(" ", ncol(mat) - 1), "\n"),
               hdest=deparse(substitute(dest)))
        Report("\n", hdest=deparse(substitute(dest)))
    }
}
##@NullChecks siena07/GlobalFunctions Resets global flags
NullChecks <- function()
{
    UserInterrupt(FALSE)
    EarlyEndPhase2(FALSE)
    UserRestart(FALSE)
    UserInterruptFlag(FALSE)
    EarlyEndPhase2Flag(FALSE)
    UserRestartFlag(FALSE)
}

##@CheckBreaks siena07/GlobalFunctions Reads global flags

CheckBreaks <- function()
{
    UserInterruptFlag(UserInterrupt())
    EarlyEndPhase2Flag(EarlyEndPhase2())
    UserRestartFlag(UserRestart())
}
