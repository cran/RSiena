#/******************************************************************************
# * SIENA: Simulation Investigation for Empirical Network Analysis
# *
# * Web: http://www.stats.ox.ac.uk/~snidjers/siena
# *
# * File: sienaDataCreateFromSession.r
# *
# * Description: This module contains the code for creation of a
# * Siena data object from an interactive session or a session file.
# *****************************************************************************/
trim.blanks <- function(x)
{
    tmp <- gsub("^ *", "", x)
    gsub(" *$", "", tmp)
}

sessionFromFile <- function(loadfilename, tk=FALSE)
{
    ## browser()
    dots <- max(c(0, gregexpr(".", loadfilename, fixed=TRUE)[[1]]))
    if (dots > 1)
    {
        extension <- substring(loadfilename, dots + 1)
        tablename <- substring(loadfilename, 1, (dots - 1))
    }
    else
    {
        extension <- ""
        tablename <- loadfilename
    }
   ## if (extension %in% c("xls"))
    ##{
     ##   suppressPackageStartupMessages(require(RODBC))
      ##  ch <- odbcConnectExcel(loadfilename)
       ## session <- sqlFetch(ch, sqlTables(ch)$TABLE_NAME[[1]], as.is=TRUE,
        ##                    nullstring="")
       ## apply(session, 2, trim.blanks)
   ## }
   ## else
    if (extension == 'csv')
    {
        session <- read.csv(loadfilename, comment.char='',
                             colClasses='character', strip.white=TRUE,
                            na.strings='')
        }
    else if (extension =='txt')
    {
        session <- read.delim(loadfilename, comment.char='',
                               colClasses='character', strip.white=TRUE )
    }
    else if (extension =='prn')
    {
        session <- read.table(loadfilename, comment.char='',
                               colClasses='character', strip.white=TRUE)
    }
    else
    {
        if (tk)
        {
            tkmessageBox(message="Can only read csv, txt",
                         "(tab delimited) or prn (space delimiited) files",
                         icon="error")
            return()
        }
        else
        {
            stop("Can only read csv, txt (tab delimited)",
                 "or prn (space delimiited) files")
        }
    }
    session
}

readInFiles <- function(session, edited)
{
    noFiles <- nrow(session)
    if (!any(edited))
    {
        files <- vector('list', noFiles)
    }
    for (i in 1:noFiles)
    {
        if (is.na(edited[i]) || !edited[i])
        {
            if (session$Type[i] == "exogenous event")
            {
                tmp <- readLines(session$Filename[i])
                changelist <- lapply(tmp, function(x)
                                 {
                                     x <- sub("^ +", "", x)
                                     x <- unlist(strsplit(x, " +"))
                                     x
                                 })
                lens <- sapply(changelist, length)
                tmp <- matrix(NA, ncol=max(lens), nrow=length(changelist))
                for (ii in 1:nrow(tmp))
                {
                    tmp[ii, 1:lens[ii]] <- changelist[[ii]]
                }
            }
            else
            {
                if (session$Format[i] == 'matrix')
                {
                    tmp <- as.matrix(read.table(session$Filename[i], as.is=TRUE))
                }
                else
                {
                    require(network)
                    tmp <- read.paj(session$Filename[i]) ## should be a single net
                }
            }
            files[[i]] <- tmp
        }
    }
    files
}
sienaDataCreateFromSession <- function (filename=NULL, session=NULL,
                                        modelName='Siena', edited=NULL)
{
    turnoffwarn <- function()
    {
        oldwarn <- getOption('warn')
        options(warn = -1)
        oldwarn
    }
    turnonwarn <- function(oldwarn)
    {
        options(warn = oldwarn)
    }
    validateNamesession <- function()
    {
        if (nrow(namesession) > 1)
        {
            if (any (namesession$ActorSet != namesession$ActorSet[1]))
                tkmessageBox(message =
                             "Actor set must be the same for one object")
            tmp <- sapply(namefiles, dim)
            if (!is.matrix(tmp) || nrow(tmp) !=2)
                tkmessageBox(message =
                             "Dimensions must be the same for one object")
            if (any (tmp[,] != tmp[,1]))
                tkmessageBox(message =
                             "Dimensions must be the same for one object")
        }
        nodeSets <- unlist(strsplit(namesession$ActorSet[1], ' '))
        if (length(nodeSets) > 2)
            tkmessageBox(message = "Invalid actor sets")
        k <- length(ActorSets)
        for (i in seq(along = nodeSets))
        {
            mymatch <- match(nodeSets[i], ActorSets)
            if (is.na(mymatch))
            {
                k <- k + 1
                ActorSets[k] <<- nodeSets[i]
                ActorSetsSize[k] <<- dim(namefiles[[1]])[i]
            }
            else if (dim(namefiles[[1]])[i] != ActorSetsSize[mymatch])
            {
                tkmessageBox(message = paste("Conflicting sizes for actor set",
                             nodeSets[i]))
            }

        }
    }
    ActorSets <- NULL
    ActorSetsSize <- NULL
    if (!is.null(filename))
    {
        session <- sessionFromFile(filename)
        session <- session[session$Selected == "Yes", ]
        edited <- rep(FALSE, nrow(session))
    }
    else
    {
        session <- session[session$Selected == "Yes", ]
        if (is.null(edited))
        {
            edited <- rep(FALSE, nrow(session))
        }
    }
    files <- readInFiles(session, edited)
    gps <- unique(session$Group)
    for (i in seq(along=gps))
    {
        ActorSets <- NULL
        ActorSetsSize <- NULL
        gpsession <- session[session$Group == gps[i], ]
        ops <- turnoffwarn()
        observations <- max(as.numeric(gpsession$Period), na.rm=TRUE)
        turnonwarn(ops)
        gpfiles <- files[session$Group == gps[i]]
        objnames <- unique(gpsession$Name)
        for (j in seq(along=objnames))
        {
            namesession <- gpsession[gpsession$Name== objnames[j], ]
            namefiles <- gpfiles[gpsession$Name== objnames[j]]
            validateNamesession()
            switch(
                   namesession$Type[1],
                   'network' = {
                       myarray <- array(NA, dim=c(dim(namefiles[[1]]),
                                             observations))
                       miss1 <- strsplit(namesession$MissingValues, " ")
                       nonzero <-  strsplit(namesession$NonZeroCode, ' ')
                       for (x in 1:nrow(namesession))
                       {
                           ##miss <- gsub(" ", "|", miss1[x], fixed=TRUE)
                           miss <- miss1[[x]]
                           ## namefiles[[x]][grep(miss, namefiles[[x]])] <- NA
                           namefiles[[x]][namefiles[[x]] %in% miss] <- NA
                           ## nonzero <- gsub(" ", "|", nonzero1[x], fixed=TRUE)
                           ## namefiles[[x]][grep(nonzero, namefiles[[x]])] <- 1
                           namefiles[[x]][!(is.na(namefiles[[x]]))
                                          & !(namefiles[[x]] %in%
                                              nonzero[[x]])] <- 0
                           namefiles[[x]][namefiles[[x]] %in% nonzero[[x]]] <- 1
                           ##  namefiles[[x]][!(is.na(namefiles[[x]]))
                           ##   & !(namefiles[[x]] %in% c(0, 1))] <- 0
                           myarray[ , , as.numeric(namesession$Period[x])] <-
                               namefiles[[x]]
                       }
                       tmp <- sienaNet(myarray, nodeSet=namesession[1,
                                                "ActorSet"])

                       assign(objnames[j], tmp, .GlobalEnv)
                   },
                   'bipartite' = {
                       nodesets <- strsplit(namesession[1, "ActorSet"], ' ')
                       myarray <- array(NA, dim=c(dim(namefiles[[1]]),
                                            observations))
                       tmp <- sapply(1:observations, function(x) {
                           ## miss <- gsub(" ", "|",
                           ##  namesession$MissingValues[x],
                           ##              fixed=TRUE)
                           miss <- namesession$MissingValues
                           if (miss != '')
                               namefiles[[x]][namefiles[[x]] %in% miss[x]] <- NA

                           ##   namefiles[[x]][grep(miss, namefiles[[x]])] <- NA
                           myarray[ , , x] <- namefiles[[x]]
                       })
                       tmp <- sienaNet(myarray, type='bipartite',
                                       nodeSet=nodesets[[1]])
                       assign(objnames[j], tmp, .GlobalEnv)
                   },
                   'behavior' = {
                       ##miss <- gsub(" ", "|",
                       ##             namesession$MissingValues[1],
                       ##              fixed=TRUE)
                       miss <- namesession$MissingValues
                       if (!is.na(miss) && miss != '')
                           namefiles[[1]][namefiles[[1]] %in% miss] <-  NA
                       ##  namefiles[[1]][grep(miss, namefiles[[1]])] <-  NA
                       assign(objnames[j],
                              sienaNet(namefiles[[1]], type = 'behavior',
                                       nodeSet=namesession[1, "ActorSet"]),
                              .GlobalEnv)
                   },
                   'constant covariate' = {
                       ##  miss <- gsub(" ", "|",
                       ##              namesession$MissingValues[1],
                       ##              fixed=TRUE)
                       ##   namefiles[[1]][grep(miss, namefiles[[1]])] <-  NA
                       miss <- namesession$MissingValues
                       namefiles[[1]][namefiles[[1]] %in% miss] <-  NA
                       varnames <- strsplit(objnames[j], ' ')
                       tmp <- sapply(1: ncol(namefiles[[1]]), function(x){
                           assign(varnames[[1]][x],
                                  coCovar(namefiles[[1]][, x],
                                          nodeSet=namesession[1,
                                          "ActorSet"]),
                                  .GlobalEnv)})
                   },
                   'changing covariate' = {
                     ##  miss <- gsub(" ", "|",
                     ##               namesession$MissingValues[1],
                     ##               fixed=TRUE)
                     ##  namefiles[[1]][grep(miss, namefiles[[1]])] <-  NA
                       miss <- namesession$MissingValues
                       namefiles[[1]][namefiles[[1]] %in% miss] <-  NA
                       assign(objnames[j],
                              varCovar (namefiles[[1]],
                                        nodeSet=namesession[1, "ActorSet"]),
                              .GlobalEnv)
                   },
                   'constant dyadic covariate' = {
                     ##  miss <- gsub(" ", "|",
                     ##               namesession$MissingValues[1],
                     ##               fixed=TRUE)
                     ##  namefiles[[1]][grep(miss, namefiles[[1]])] <-  NA
                       miss <- namesession$MissingValues
                       namefiles[[1]][namefiles[[1]] %in% miss] <-  NA
                       if (namesession[1, "ActorSet"] == "Actors")
                       {
                           namesession[1, "ActorSet"]<- "Actors Actors"
                       }
                       nodesets <- strsplit(namesession[1, "ActorSet"], ' ')
                       assign(objnames[j],
                              coDyadCovar (namefiles[[1]],
                                           nodeSets=nodesets[[1]]),
                              .GlobalEnv)
                   },
                   'changing dyadic covariate' = {
                       if (namesession[1, "ActorSet"] == "Actors")
                       {
                           namesession[1, "ActorSet"]<- "Actors Actors"
                       }
                       nodesets <- strsplit(namesession[1, "ActorSet"], ' ')
                       myarray <- array(NA, dim=c(dim(namefiles[[1]]),
                                             observations - 1))
                       tmp <- sapply(1:nrow(namesession), function(x){
                           ##      miss <- gsub(" ", "|",
                           ##                   namesession$MissingValues[x],
                           ##                    fixed=TRUE)
                           miss <- namesession$MissingValues
                           if (miss != '')
                           {
                               ##namefiles[[x]][grep(miss,namefiles[[x]])] <- NA
                               namefiles[[x]][namefiles[[x]] %in% miss[x]] <- NA
                           }
                           myarray[ , ,as.numeric(namesession$Period[x])] <-
                               namefiles[[x]]
                       })
                       tmp <- varDyadCovar(myarray, nodeSets=nodesets[[1]])
                       assign(objnames[j], tmp, .GlobalEnv)
                   },
                   'exogenous event' = {
                       tmp <- namefiles[[1]]

                       clist <- tapply(tmp, row(tmp), function(x)
                                  {
                                      if (any(is.na(x)))
                                      {
                                          firstNA <- min(which(is.na(x)))
                                          lastNA <- max(which(is.na(x)))
                                          if (any(!is.na(x[firstNA:lastNA])))
                                              stop("gaps in values in ",
                                                   "exogenous event file")
                                          if (lastNA < length(x))
                                          {
                                              stop("Missing data at start ",
                                                   "of exogenous event file")
                                          }
                                          x<- x[1:(firstNA - 1)]
                                      }
                                      x
                                  })
                       tmp <- sienaCompositionChange(clist,
                                                     namesession$ActorSet[[1]])
                       assign(objnames[j], tmp, .GlobalEnv)
                   },
               {
                   if (is.null(filename))
                   {
                       tkmessageBox(message = paste('File of unknown type:',
                                    gpsession$type))
                   }
                   else
                   {
                       stop(paste('File of unknown type:',
                                   gpsession$type))
                   }
                   return()
               }
                   )
        }
        ## create the node sets
        tmp <- lapply(1:length(ActorSets), function(x){
            assign(ActorSets[x], sienaNodeSet(ActorSetsSize[x],
                                              nodeSetName=ActorSets[x]),
                   .GlobalEnv)
        })
        ## create the group object
        obj <- unlist(lapply(objnames, strsplit, split=" "))
        if (any(duplicated(obj)))
        {
            if (is.null(filename))
            {
                tkmessageBox(message = paste('Duplicate names',
                             obj[duplicated(obj)]))
            }
            else
                stop(paste('Duplicate names',
                             obj[duplicated(obj)]))
        }
        objlist <- mget(obj, .GlobalEnv)
        nodeSetList <- mget(ActorSets,.GlobalEnv)
        names(nodeSetList) <- NULL
        arglist <- c(objlist, nodeSets=nodeSetList)
        assign(gps[i], do.call(sienaDataCreate,
                               c(objlist, nodeSets=list(nodeSetList))),
               .GlobalEnv)
    }
    ##join the groups
    if (length(gps) > 1)
    {
        mydata <- sienaGroupCreate(mget(gps, .GlobalEnv))
    }
    else
    {
        mydata <- get(gps)
    }
    myeff <- getEffects(mydata)
    print01Report(mydata, myeff, modelName, session)
    return(list(OK=TRUE, mydata=mydata, myeff=myeff))
}
