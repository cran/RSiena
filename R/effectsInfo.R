symmetricRateEffects <-
structure(list(EffectName = c("basic rate parameter", "constant rate (period ",
"degree effect on rate", "indegree effect on rate", "reciprocity effect on rate",
"effect 1/degree on rate"), FunctionName = c("Amount of change",
"Amount of change in period", "Amount of change x outdegrees",
"Amount of change x indegrees", "Amount of change x reciprocity",
"Amount of change x (1/outdegrees)"), ShortName = c("Rate", "Rate",
"degreeRate", "indegRate", "recipRate", "degRateInv")), .Names = c("EffectName",
"FunctionName", "ShortName"), class = "data.frame", row.names = c(NA,
-6L))
nonSymmetricRateEffects <-
structure(list(EffectName = c("basic rate parameter", "constant rate (period ",
"outdegree effect on rate", "indegree effect on rate", "reciprocity effect on rate",
"effect 1/outdegree on rate"), FunctionName = c("Amount of change",
"Amount of change in period", "Amount of change x outdegrees",
"Amount of change x indegrees", "Amount of change x reciprocity",
"Amount of change x (1/outdegrees)"), ShortName = c("Rate", "Rate",
"outRate", "inRate", "recipRate", "outRateInv")), .Names = c("EffectName",
"FunctionName", "ShortName"), class = "data.frame", row.names = c(NA,
-6L))
nonSymmetricObjEffects <-
structure(list(EffectName = c("outdegree (density)", "reciprocity",
"transitive triplets", "transitive mediated triplets", "3-cycles",
"transitive ties", "betweenness", "balance", "number of actors at distance 2",
"number pairs at doubly achieved distance 2", "dense triads",
"indegree - popularity", "indegree - popularity (sqrt)", "outdegree - popularity",
"outdegree - popularity (sqrt)", "indegree - activity", "indegree - activity (sqrt)",
"outdegree - activity", "outdegree - activity (sqrt)", "1/(outdegree + #)",
"1/(outdegree+#)(outdegree+1+#)", "out-out degree^(1/2) assortativity",
"out-in degree^(1/2) assortativity", "in-out degree^(1/2) assortativity",
"in-in degree^(1/2) assortativity"), FunctionName = c("Number of ties",
"Number of reciprocated ties", "Number of transitive triplets",
"Number of transitive triplets", "3-cycles", "Number of ties with transitive closure",
"betweenness count", "Amount of balance", "Number of directed distances equal to 2",
"Number of doubly achieved distances 2", "Sum of triads with # or more ties",
"Sum of squared indegrees", "Sum of indegrees x sqrt(indegree)",
"Sum of crossproducts indegree x outdegree", "Sum of indegrees x sqrt(outdegree)",
"Sum of crossproducts indegree x outdegree", "Sum of outdegrees x sqrt(indegree)",
"Sum of squared outdegrees", "Sum of outdegrees^(1.5)", "Sum 1/(outdegrees + #)",
"Sum 1/(outdegrees + #)(outdegrees + 1 + #)", "Sum of out-out degree^(1/2) products",
"Sum of out-in degree^(1/2) products", "Sum of in-out degree^(1/2) products",
"Sum of in-in degree^(1/2) products"), Endowment. = c(TRUE, TRUE,
TRUE, TRUE, TRUE, TRUE, FALSE, TRUE, TRUE, TRUE, FALSE, TRUE,
TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE,
TRUE, TRUE, TRUE, TRUE), ShortName = c("density", "recip", "transTrip",
"transMedTrip", "cycle3", "transTies", "between", "balance",
"nbrDist2", "nbrDist2twice", "denseTriads", "inPop", "inPopSqrt",
"outPop", "outPopSqrt", "inAct", "inActSqrt", "outAct", "outActSqrt",
"outInv", "outSqInv", "outOutAss", "outInAss", "inOutAss", "inInAss"
), parm = c(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 5, 0, 0, 0, 0, 0, 0,
0, 0, 1, 1, 2, 2, 2, 2)), .Names = c("EffectName", "FunctionName",
"Endowment.", "ShortName", "parm"), row.names = c(NA, -25L), class = "data.frame")
symmetricObjEffects <-
structure(list(EffectName = c("degree (density)", "transitive triads",
"transitive ties", "betweenness", "balance", "number of actor pairs at distance 2",
"number pairs at doubly achieved distance 2", "degree of alter",
"sqrt degree of alter", "degree^(1.5)", "1/(degree + #)", "1/(degree+#)(degree+1+#)",
"degree^(1/2) assortativity"), FunctionName = c("Number of edges",
"Number of transitive triads", "Number of ties wth transitive closure",
"betweenness count", "Amount of balance", "Number of distances equal to 2",
"Number of doubly achieved distances 2", "Sum of squared degrees",
"Sum of degrees ", "Sum of degrees^(1.5)", "Sum 1/(degrees + #)",
"Sum 1/(degrees + #)(degrees + 1 + #)", "Sum of degree^(1/2) products"
), Endowment. = c(TRUE, TRUE, TRUE, FALSE, TRUE, FALSE,
FALSE, TRUE, TRUE, FALSE, FALSE,
FALSE, TRUE), ShortName = c("density", "transTriads", "transTies",
"between", "balance", "nbrDist2", "nbrDist2Twice", "inPop",
"inPopSqrt", "outActSqrt", "outInv", "outSqInv", "outOutAss"), parm = c(0,
0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 2)), .Names = c("EffectName",
"FunctionName", "Endowment.", "ShortName", "parm"), row.names = c(NA,
13L), class = "data.frame")
behaviorObjEffects <-
structure(list(EffectName = c("linear shape", "quadratic shape",
"average similarity", "total similarity", "indegree", "outdegree",
"isolate", "ave. sim. x reciprocity", "tot. sim. x reciprocity",
"ave. sim. x popularity alter", "tot. sim. x popularity alter",
"ave. sim. x rec. x pop. (alter)", "tot. sim. x rec. x pop. (alter)",
"average alter", "average rec. alters", "dense triads <maybe wrong>",
"similarity in dense triads <maybe wrong>", "reciprocated degree",
"ave. sim. x popularity ego"), Function.Name = c("cent. sum                                      ",
"sum of cent. squares                           ", "average similarity                             ",
"total similarity                               ", "indegrees                                      ",
"outdegrees                                     ", "isolate                                        ",
"ave. similarity x reciprocity                  ", "tot. similarity x reciprocity                  ",
"ave. sim. x indegrees(one-sided)               ", "tot. sim. x indegrees(one-sided)               ",
"ave. sim. x rec. x i.d.(one-sided)             ", "tot. sim. x rec. x i.d.(one-sided)             ",
"average alters                                 ", "average rec. alters                            ",
"dense triads <<<maybe wrong>>>                 ", "homogeneity of dense triads <<<maybe wrong>>>  ",
"reciprocated degrees                           ", "ave. sim. x indegrees ego                      "
), Endowment. = c(TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE,
TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE
), ShortName = c("linear", "quad", "avSim", "totSim", "indeg",
"outdeg", "isolate", "avSimRecip", "totSimRecip", "avSimPopAlt",
"totSimPopAlt", "avSimRecPop", "totSimRecPop", "avAlt", "avRecAlt",
"denseTriads", "simDenseTriads", "recipDeg", "avSimPopEgo")), .Names = c("EffectName",
"Function.Name", "Endowment.", "ShortName"), class = "data.frame", row.names = c(NA,
-19L))
behaviorRateEffects <-
structure(list(EffectName = c("rate  (period ", "outdegree effect on rate",
"indegree effect on rate", "reciprocated effect on rate"), FunctionName = c("Amount of behavioral change on",
"x outdegree", "x indegree", "x reciprocity"), ShortName = c("rate",
"outRate", "inRate", "recipRate")), .Names = c("EffectName",
"FunctionName", "ShortName"), class = "data.frame", row.names = c(NA,
-4L))
covarBehObjEffects <-
structure(c("effect from", "influence interaction? x", "x", "influ. int. possible x",
"effFrom", "inflIntX"), .Dim = 2:3)
covarBehObjInteractions <-
structure(c("av.sim. x ", "tot. sim. x ", "av. alters x ", "avSimX",
"totSimX", "avAltX"), .Dim = c(3L, 2L))
dyadObjEffects <-
structure(c("WW=>X closure of", "WX=>X closure of", "XW=>X closure of",
"WWX", "WXX", "XWX"), .Dim = c(3L, 2L))
covarNonSymmetricObjEffects <-
structure(c("alter", "squared alter", "ego", "similarity", "similarity x reciprocity",
"Sum of indegrees x", "Sum of indegrees x squared", "Sum of outdegrees x",
"Similarity on", "Similarity x reciprocity on", "altX", "altSqX",
"egoX", "simX", "simRecipX"), .Dim = c(5L, 3L))
covarSymmetricObjEffects <-
structure(c("", "squared", "similarity", "X", "sqX", "simX"), .Dim = c(3L,
2L))
