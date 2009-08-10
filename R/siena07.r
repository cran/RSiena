##/******************************************************************************
## * SIENA: Simulation Investigation for Empirical Network Analysis
## *
## * Web: http://stat.gamma.rug.nl/siena.html
## *
## * File: siena07.r
## *
## * Description: This file contains the main controlling module for the model
## * fitting
## * Also contains utility functions
## *****************************************************************************/

outf <- NULL
lf <- NULL
bof <- NULL
cf <- NULL

siena07<- function(x, batch = FALSE, verbose = FALSE, useCluster = FALSE,
                   nbrNodes = 2, initC=FALSE,
                   clusterString=rep("localhost", nbrNodes), tt=NULL,
                   parallelTesting=FALSE, ...)
{
    exitfn <- function()
    {
       if (!is.batch())
       {
           tkdestroy(tkvars$tt)
       }
       ## close the report file
       Report(close=TRUE)
    }
    on.exit(exitfn())


    time0 <-  proc.time()['elapsed']
    z <- NULL ## z is the object for all control information which may change.
    ## x is designed to be readonly. Only z is returned.

    if (useCluster)
    {
        require(snow)
        require(rlecuyer)
        x$firstg <- x$firstg * nbrNodes
        z$int <- nbrNodes
    }
    else
    {
        z$int <- 1
    }
    if (parallelTesting)
    {
        set.seed(1, kind='Wich')
        ## randomseed2 is for second generator needed only for parallel testing
        randomseed2 <- .Random.seed
      #  .Random.seed[2:4] <- as.integer(c(1,2,3))
        randomseed2[2:4] <- as.integer(c(3,2,1))
        seed <- 1
        newseed <- 1
    }
    else
    {
        randomseed2 <-  NULL
        ## x$randomSeed is the user seed, if any
        if (!is.null(x$randomSeed))
        {
            set.seed(x$randomSeed)
            seed <- x$randomSeed
       }
        else
        {
            if (exists(".Random.seed"))
            {
                rm(.Random.seed, pos=1)
            }
            newseed <- trunc(runif(1) * 1000000)
            set.seed(newseed)  ## get R to create a random number seed for me.
            seed <- NULL
        }
    }
    z$randomseed2 <- randomseed2

    ## set the global is.batch
    is.batch(batch)

    ## open the output file
    Report(open=TRUE, projname=x$projname, verbose=verbose)
    InitReports(seed, newseed)

    ## reset the globals for interrupts
    NullChecks()

    ## create the screen
    if (!is.batch())
    {
        tkvars <- runtk(tt=tt)
        z$tkvars<- tkvars
        z$pb <- list(pb=tkvars$pb, pbval=0, pbmax=1)
    }
    else
    {
        z$pb <- list(pb=NULL, pbval=0, pbmax=1)
    }

    z <- robmon(z, x, useCluster, nbrNodes, initC, clusterString,...)

    time1 <-  proc.time()['elapsed']
    Report(c("Total computation time", round(time1 - time0, digits=2),
             "seconds.\n"), outf)

    if (useCluster)
        stopCluster(z$cl)

    class(z) <- "sienaFit"
    z
}
InitReports <- function(seed, newseed)
{
    Report("\n\n-----------------------------------\n", outf)
    Report("New Analysis started.\n", outf)
    Report(c("Date and time:", format(Sys.time(),"%d/%m/%Y %H:%M:%S")), outf)
    Report("\nNew results follow.\n", outf)
    Report("-----------------------------------\n", outf)
    rforgeRevision <-  packageDescription("RSiena",
                                          fields="Repository/R-Forge/Revision")
    if (is.na(rforgeRevision))
    {
        revision <- ""
    }
    else
    {
        revision <- paste(" R-forge revision: ", rforgeRevision, " ", sep="")
    }
    Report(c("\nSiena version ",
             packageDescription("RSiena", fields = "Version"), " (",
             format(as.Date(packageDescription("RSiena", fields = "Date")),
                    "%d %b %y"), ")",
             revision, "\n\n"), sep = '',  outf )
    Heading(1, outf, "Estimation by stochastic approximation algorithm.")
    if (is.null(seed))
    {
        Report("Random initialization of random number stream.\n", outf)
        Report(sprintf("Current random number seed is %d.\n", newseed), outf)
    }
    else
    {
        Report(sprintf("Current random number seed is %d.\n", seed), outf)
    }
}
WriteOutTheta <- function(z)
{
    if (!is.batch())
    {
        tkdelete(z$tkvars$current, "1.0", "end")
        tmp <- paste(c("", rep("\n", z$pp - 1)),
                    format(round(z$theta,4), width=12, sep=""),
                    collapse="")
        tkinsert(z$tkvars$current, "1.0", tmp)
    }
    else
    {
        Report(c("theta:", format(z$theta, digits=3), "\n"))
    }
    Report("Current parameter values:\n", cf)
    Report(format(z$theta), cf, fill=80)
}

DisplayTheta<- function(z)
{
    if ((z$Phase == 2 || z$nit == 1 ) && (z$nit <= 30))
    {
        if (!is.batch())
        {
            tkdelete(z$tkvars$current, "1.0", "end")
            tmp<- paste(c("", rep("\n", z$pp - 1)),
                        format(z$theta, width=12, sep=""),
                    collapse="")
            tkinsert(z$tkvars$current, "1.0", tmp)
        }
        else
        {
          Report(c("theta:", format(z$theta, digits=3), "\n"))
        }
    }

}

NullChecks <- function()
{
    UserInterrupt(FALSE)
    EarlyEndPhase2(FALSE)
    UserRestart(FALSE)
    UserInterruptFlag(FALSE)
    EarlyEndPhase2Flag(FALSE)
    UserRestartFlag(FALSE)
}

CheckBreaks <- function()
{
    UserInterruptFlag(UserInterrupt())
    EarlyEndPhase2Flag(EarlyEndPhase2())
    UserRestartFlag(UserRestart())
}

AnnouncePhase <- function(z, x, subphase=NULL)
{
    if (!is.batch())
    {
        tkdelete(z$tkvars$phase, 0, "end")
        tkinsert(z$tkvars$phase, 0, paste(" ", z$Phase))
        tkdelete(z$tkvars$subphase, 0, "end")
        tkdelete(z$tkvars$iteration, 0, "end")
        tkinsert(z$tkvars$iteration, 0, format(0, width=6))
    }
    if (missing(subphase))
    {
        Report(c("\nStart phase", z$Phase, "\n"), cf)
    }
    else
    {
        if (!is.batch())
        {
            tkinsert(z$tkvars$subphase, 0, paste(" ", subphase))
        }
        Report(c("\nStart phase ", z$Phase, ".", subphase, "\n"), sep="", cf)
    }
    if (z$Phase == 0)
    {
        if (!is.batch())
        {
            tkconfigure(z$tkvars$current, height=z$pp)
            tkconfigure(z$tkvars$deviation, height=z$pp)
            tkconfigure(z$tkvars$quasi, height=z$pp)
        }
        n1pos <- z$n1 * (z$pp + 1)
        z$n2min0 <- 7 + z$pp
        z$n2minimum<- rep(0, x$nsub)
        z$n2maximum<- rep(0, x$nsub)
    ## 2.5198421 = 2^(4/3); this gives a gain parameter of order n^(-3/4) ##
        if (x$nsub>0)
        {
            z$n2minimum[1] <- trunc(z$n2min0 * 2.52)
            z$n2maximum[1] <- z$n2minimum[1] + 200
            if (x$nsub > 1)
            {
                for (i in 2:x$nsub)
                {
                    z$n2minimum[i] <- trunc(z$n2minimum[i-1] * 2.52)
                    z$n2maximum[i] <- z$n2minimum[i] + 200
                }
            }
        }
        z$n2partsum <- c(0, cumsum(z$n2maximum))
        n2sum <- sum(z$n2maximum)
        ##Progress bar
        pbmax <- n1pos + n2sum + x$n3
        z$n1pos<- n1pos
        if (!x$maxlike && z$FinDiff.method)
            pbmax <- pbmax + x$n3 * z$pp
        z$pb$pbval<- 0
        z$pb <- createProgressBar(z$pb, maxvalue=pbmax)
        z$pb$pbmax<- pbmax
   }
    if (z$Phase==2)
    {
        propo <- z$n1pos + z$n2partsum[subphase]
        if (propo> getProgressBar(z$pb))
            z$pb<-setProgressBar(z$pb,propo)
    }
    if (z$Phase ==3)
    {
        propo <- z$n1pos + z$n2partsum[x$nsub + 1]
        if (!z$AllUserFixed)
            z$pb<- setProgressBar(z$pb,propo)
       else
        {
            max <- x$n3
            z$pb<-createProgressBar(z$pb,max)
       }
    }
    z
}

roundfreq<- function(w)
{
    vec1 <- c(1, 2, 3, 4, 31, 66, 101, 300, 500)
    vec2 <- c(1, 2, 3, 20, 50, 100, 200, 500)
    if (is.batch())
        w <- vec2[findInterval(w, vec1, all.inside=TRUE)]
    else
        w <- vec2[findInterval(w, vec1[1:7], all.inside=TRUE)]
    w
}

model.create<- function(fn=simstats0c, usesimstats0c=TRUE,
                        projname="Siena", MaxDegree=0, useStdInits=FALSE,
                        n3=1000, nsub=4, maxlike=FALSE, diag=TRUE,
                        condvarno=0, condname='',
                        firstg=0.2, cond=FALSE, findiff=FALSE,  seed=NULL)
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


Report <- local({verbose <-  NULL;
                 Reportfun(list(outf=outf, lf=lf, cf=cf, bof=bof), verbose)})

UserInterrupt <- local({A <-  FALSE;function(x){if (!missing(x))A<<-x;A}})
EarlyEndPhase2 <- local({A <-  FALSE;function(x){if (!missing(x))A<<-x;A}})
UserRestart <- local({A <-  FALSE;function(x){if (!missing(x))A<<-x;A}})
UserInterruptFlag <- local({A <-  FALSE;function(x){if (!missing(x))A<<-x;A}})
EarlyEndPhase2Flag <- local({A <-  FALSE;function(x){if (!missing(x))A<<-x;A}})
UserRestartFlag <- local({A <-  FALSE;function(x){if (!missing(x))A<<-x;A}})
is.batch <- local({A <-  FALSE;function(x){if (!missing(x))A<<-x;A}})
DONE <- local({A <-  FALSE;function(x){if (!missing(x))A<<-x;invisible(A)}})
FRANstore <- local({A <-  NULL;function(x){if (!missing(x)) A<<-x;A}})

DisplayThetaAutocor <- function(z)
{
    if (!is.batch())
    {
        tkdelete(z$tkvars$current, "1.0", "end")
        tmp<- paste(c("", rep("\n", z$pp - 1)),
                    format(round(z$theta, 4), width=12, sep=""),
                    collapse="")
        tkinsert(z$tkvars$current, "1.0", tmp)
        tkdelete(z$tkvars$quasi, "1.0", "end")
        tmp<- paste(c("", rep("\n", z$pp - 1)),
                    format(round(z$ac, 4), width=12, sep=""),
                    collapse="")
        tkinsert(z$tkvars$quasi, "1.0", tmp)
    }
    else
    {
        Report(c("theta", format(z$theta, digits=3),"\n"))
        Report(c("ac", format(z$ac, digits=3), "\n"))
  }

}
DisplayandWritetheta <- function(z)
{
    if (!is.batch())
    {
        tkdelete(z$tkvars$current, "1.0", "end")
        tmp<- paste(c("", rep("\n", z$pp - 1)),
                    format(round(z$theta, 4), width=12, nsmall=4, sep=""),
                    collapse="")
        tkinsert(z$tkvars$current, "1.0", tmp)
    }
    else
    {
        Report(c("theta", format(z$theta, digits=3), "\n"))
    }
}
DisplayTheta <- function(z)
{
        if (!is.batch())
        {
            tkdelete(z$tkvars$current, "1.0", "end")
            tmp<- paste(c("", rep("\n", z$pp - 1)),
                        format(round(z$theta, 4), width=12, sep="", nsmall=4),
                        collapse="")
            tkinsert(z$tkvars$current, "1.0", tmp)
        }

}
DisplayDeviations <- function(z, fra)
{
        if (!is.batch())
        {
            tkdelete(z$tkvars$deviations, "1.0", "end")
            tmp<- paste(c("", rep("\n", z$pp - 1)),
                        format(round(fra, 4), width=12, sep="", nsmall=4),
                        collapse="")
            tkinsert(z$tkvars$deviations, "1.0", tmp)
        }
}
DisplayIteration <- function(z)
{
    if (!is.batch())
    {
        tkdelete(z$tkvars$iteration, 0, "end")
        tkinsert(z$tkvars$iteration, 0, format(z$nit, width=6))
        tcl("update")
    }
}

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

PrtOutMat<- function(mat,dest)
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
Root<- function(x)
{
    ifelse(abs(x) > 1e-36, sqrt(abs(x)), 1e-18)
}

getProgressBar <- function(pb)
{
    if (is.batch())
        val <- pb$pbval
    else
        val <- as.numeric(tclvalue(tkcget(pb$pb, "-value")))
    val
}

setProgressBar <- function(pb, val)
{
    if (is.batch())
    {
        pb$pbval <- val
    }
    else
    {
        tkconfigure(pb$pb, value=val)
        tcl("update")
    }
    pb
}
createProgressBar <- function(pb, maxvalue)
{
    if (is.batch())
        pb$pbmax <- maxvalue
    else
        tkconfigure(pb$pb, maximum=maxvalue)
    pb
}

tkErrorMessage <- function()
{
    tkmessageBox(geterrmessage(), icon="error")
}

errorHandler <- function()
{
    opts <- options()
    if (!is.batch())
    {
        options(show.error.messages=FALSE)
        options(error=tkErrorMessage)
    }


}
