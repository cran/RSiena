DataReport <- function(z, x, f)
{
    ##f could be a group, but has attributes like a group even if not!
    oneMode <- attr(f, "types") == "oneMode"
    behavior <- attr(f, "types") == "behavior"
    nOneMode <- sum(oneMode)
    nBehavior <- sum(behavior)
    oneModeNames <- attr(f, "netnames")[oneMode]
    behaviorNames <- attr(f, "netnames")[behavior]
    symmetric <- attr(f, "symmetric")[oneMode]
    nDepVars <- nOneMode + nBehavior
    observations <- attr(f, "observations") ##note this is total number of
    ##  periods to process
    exogenous <- attr(f, 'compositionChange')
    exoOptions <- attr(f, 'exoptions')
    for (i in 1:nOneMode)
    {
        if (nOneMode > 1)
        {
            Report(sprintf("Network %d %s\n", i, oneModeNames[i]), outf)
        }
        if (symmetric[i])
            ModelTypeStrings <- c("Forcing model")
        else
            ModelTypeStrings <- c("Standard actor-oriented model")
        Report(sprintf("Model Type %d: %s\n",
                       x$ModelType, ModelTypeStrings[x$ModelType]), outf)
    }
    Report("Estimation method: ", outf)
    if (x$cconditional)
    {
        Report("conditional moment estimation\n", outf)
        Report(c('Conditioning variable is the total number of observed',
                 'changes ("distance") \n'), outf)
        if (z$condtype == 'oneMode')
        {
            if (nOneMode == 1)
            {
                Report('in the network variable.\n', outf)
            }
            else
            {
               Report(c("in network variable", x$condname, '\n'), outf)
            }
        }
        else
        {
          Report(c("in behavioral variable", x$condname, '\n'), outf)
        }
        if (observations == 1)
        {
            ## note can only be one group
            Report(c("Distance for simulations is",
                     format(f[[1]]$distances, width=4), '.\n'), outf)
        }
        else
        {
             Report("Distances for simulations are\n", outf)
 #            Report(c("period   :", format(1:observations, width=4),
               Report(c("period   :", sapply(f[1], function(x)
                                             sapply(x$depvars, function(y)
                                             format(1:observations, width=4))),
                    '\n'), sep='', outf)
             Report(c("distance : ", format(sapply(f[1],
                                                   function(x)x$distances),
                      width = 4), '.\n'),sep='',  outf)
        }
    }
    else
    {
        Report("unconditional moment estimation.\n", outf)

        if (exogenous)
        {
            Report("Changing composition: no conditional moment estimation.\n",
                   outf)
            ## TODO report separately for different node sets if appropriate
            Report(c("Joiners/leavers option:", exoOptions, "\n"), sep="",
                   outf)
        }
        Report('\nTime duration for simulations ', outf)
        if (observations >= 2)
        {
            Report('in each period ', outf)
        }
        Report('is 1.0.\n', outf)
    }
    if (z$FinDiff.method)
    {
        Report(c('Standard errors are estimated with the finite difference',
                 'method.\n'), outf)
        if (!x$FinDiff.method)
        {
            Report("Note that the option requested has been over-ridden\n",
                   outf)
        }
    }
    else
    {
        Report(c('Standard errors are estimated with the likelihood ratio',
               'method.\n'), outf)
    }
    ## any max degree restrictions?
    if (any(x$MaxDegree > 0))
    {
        for (i in 1:length(x$MaxDegree))
        {
            if (x$MaxDegree[i])
            {
                mdnet <- names(x$MaxDegree)[i]
                Report(c("Dependent network variable", mdnet, ':\n'), outf)
                if (attr(f, 'symmetric')[match(mdnet), attr(f, "netnames")])
                {
                    Report(c("All graphs are constrained to having degrees not",
                             "larger than", x$MaxDegree[i], '.\n'), outf)
                }
                else
                {
                    Report(c("All digraphs are constrained to having",
                             "outdegrees not larger than", x$MaxDegree[i],
                             '.\n'), outf)
                }

            }
        }
    }
    Report(sprintf("Initial value of gain parameter is %10.7f.\n", x$firstg),
           outf)
    Report(sprintf("Number of subphases in Phase 2 is %d.\n\n", x$nsub), outf)
    Report('Initial parameter values are \n', outf)
    if (x$cconditional)
    {
        if (observations == 1)
        {
            if (nDepVars == 1)
            {
                Report(format('  0. Rate parameter', width = 43), outf)
            }
            else
            {
                Report(format('  0. Rate parameter of conditioning variable',
                              width = 43), outf)
            }
            Report(sprintf("%9.4f\n", attr(f, "condEffects")$initialValue[1]),
                   outf)
        }
        else ## observations > 1
        {
            for (i in 1:observations)
            {
                if (nDepVars == 1)
                {
                    Report(format(paste('  0.', format(format(i), width=2),
                                        ' Rate parameter', sep=''),
                                  width = 43), outf)
                }
                else
                {
                    Report(format(paste('  0.', format(format(i), width=2),
                                        ' Rate parameter cond.',
                                  'variable period', i, sep=''),
                                  width = 43), outf)
                }
                Report(sprintf("%9.4f\n",
                               attr(f, "condEffects")$initialValue[i]), outf)
            }
        }
    }

    fixed <- ifelse(z$fixed, '  (fixed) ', '')
    tmp <- paste(sprintf("%3d",1:length(z$effects$effectName)), '. ',
                 format(paste(z$effects$type, ':  ', z$effects$effectName,
                              sep = ''), width = 52),
                 sprintf("%9.4f", z$effects$initialValue), fixed, '\n',
                 sep = '', collapse = '')
    Report(tmp, outf)
    ## targets:
    Report("\n\nObserved values of target statistics are\n", outf)
    tmp <- paste(sprintf("%3d",1:length(z$effects$effectName)), '. ',
                 format(z$effects$functionName, width = 66),
                 sprintf("%9.4f",
                         ifelse(z$effects$type=='endow', -z$targets,
                                z$targets)),
                 '\n', sep = '', collapse = '')
    Report(tmp, outf)
    Report(c('\n', nrow(z$effects), 'parameters,', nrow(z$effects),
             'statistics\n'),outf)
}
