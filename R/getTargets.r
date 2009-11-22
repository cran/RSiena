##@getTargets Miscellaneous Written for Krista. Use as RSiena:::getTargets
getTargets <- function(data, effects)
{
    f <- unpackData(data)
    effects <- effects[effects$include,]
    ##dyn.load(dllpath)
    ##
    ## dyn.load('d:/sienasvn/siena/src/RSiena.dll')
    pData <- .Call('setupData', PACKAGE="RSiena",
                   list(as.integer(f$observations)),
                   list(f$nodeSets))
    ## register a finalizer
    ans <- reg.finalizer(pData, clearData, onexit = FALSE)
    ans<- .Call('OneMode', PACKAGE="RSiena",
                pData, list(f$nets))
    ans<- .Call('Behavior', PACKAGE="RSiena", pData,
               list(f$behavs))
    ans<-.Call('ConstantCovariates', PACKAGE="RSiena",
               pData, list(f$cCovars))
    ans<-.Call('ChangingCovariates',PACKAGE="RSiena",
               pData,list(f$vCovars))
    ans<-.Call('DyadicCovariates',PACKAGE="RSiena",
               pData,list(f$dycCovars))
    ans<-.Call('ChangingDyadicCovariates',PACKAGE="RSiena",
               pData, list(f$dyvCovars))
    storage.mode(effects$parm) <- 'integer'
    storage.mode(effects$group) <- 'integer'
    storage.mode(effects$period) <- 'integer'
    effects$effectPtr <- NA
    myeffects <- split(effects, effects$name)
    ans<- .Call('effects', PACKAGE="RSiena",
                pData, myeffects)
    pModel <- ans[[1]][[1]]
        for (i in 1:length(ans[[2]])) ## ans[[2]] is a list of lists of
            ## pointers to effects. Each list corresponds to one
            ## dependent variable
        {
            effectPtr <- ans[[2]][[i]]
            myeffects[[i]]$effectPtr <- effectPtr
        }
    ans <- .Call('getTargets', PACKAGE="RSiena",
                 pData, pModel, myeffects)
    ans
}
