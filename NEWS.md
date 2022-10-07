# RSiena 1.3.13

## 
   
2022-10-07

## Changes in RSiena:  

### Updates:
  * Replacements in EffectFactory.cpp of single | operator by |.

# RSiena 1.3.12

## 
   
2022-10-06

## Changes in RSiena:  

### Updates:
  * Changes to comply with new version of `Matrix` package.
  * Replacements in some C++ functions of single & and | operators by && and ||.
### Corrections:  
  * `universalOffset` initialized as 0; it was earlier initialized as
    the maximum real number (`NetworkLongitudinalData.cpp`). 
  * `thetaStore` deleted (was trash in `phase2.r`).
  * Various comparisons for vectors with 0 changed to using `all`
    to avoid warnings (`initializeFRAN.r`).
### Code modifications:
  * `sigmas` and `meansigmas` added to `sienaRI` object.
  * Print of standard deviations in the `sienaRI` object for `printSigma=TRUE` 
    changed to using averages at the variance level.
  * If `returnThetas` in the call of `siena07`, also simulated estimation statistics
    during Phase 2 (deviations from targets) are returned.
### Effects:
  * Several new effects related to primary setting:
    `nonPCompress`, `primCompress`, `primary`, `primDegAct`,
    `primDegActDiff`, `primDegActDiffSqrt`, `primDegActSqrt`,
    `primDegActLog`, `primDegActInv`.
  * `gwdspFB` effect added for two-mode networks.
  * New effects `outAct_ego`, `inAct_ego`,`reciAct_ego`, `toAny`.
  * For effects `to`, `toBack`, `toRecip`, `mixedInXW`, 
    internal effect parameter 3 now specifies truncation of the number of 
    twosteps (change to `MixedTwoStepFunction`). 
### Improvements of documentation:
  * Modified help page for `sienaRI`.
  * Small modifications of help page for `sienaGOF`.

# RSiena 1.3.11

## 
   
2022-05-30

## Changes in RSiena:  

### Corrections:   
  * Correction in `effects.r` of error that led to warning
    for multivariate networks.
  * Correction of help page for `sienaGOF` (`groupName`).
  * Correction of `igraphNetworkExtraction` in the help page for
    `sienaGOF-auxiliary`.

### Improvements of functionality: 
  * Further explanation of `mixedTriadCensus` in the help page for
    `sienaGOF-auxiliary`.


# RSiena 1.3.10

## 
   
2022-04-28

## Changes in RSiena:  

### Corrections:
  * Bug corrected that occurred when several two-mode networks were included
    in the dependent variables, with an order restriction between them.
    (Correction of `HigherFilter` and `DisjointFilter`).

### Effects:
  * New effects `avInAltW`, `avWInAlt`, `totInAltW`, `totWInAlt` 
    (with help from Robert Krause).
  * Corrected implementation of `sharedTo`.
  
### Code modifications:
  * Several modifications to enable traceback of errors occurring
    in checkSenderRange called in inTies.

# RSiena 1.3.9

## 
   
2022-03-18

## Changes in RSiena:  

### Effects:
  * Corrected implementation of `simAllNear` and `simAllFar`.

### Corrections:
  * small correction of `summary.sienaGOF`.
  * small correction of `sienaTimeTest`.

# RSiena 1.3.8

## 
   
2022-03-07

## Changes in RSiena:  

### Effects:
  * Changed internal effect parameter for `simAllNear` to 2 and for
    `simAllFar` to 4.

### Improvements of functionality: 
  * In `sienaTimeTest`, added `warn=FALSE` to `varCovar()` to avoid warnings.
  * Small improvements to help pages for `sienaGroupCreate` and `sienaGOF`.

### Corrections:
  * corrected sienaRI for behavioral variables.
    This required changes in 
    `StatisticCalculator::calculateBehaviorStatistics` and
    `StatisticCalculator::calculateBehaviorGMMStatistics`.
  * dropped exclusion of bipartite for sienaRI (only continuous excluded),
    but only if there are fewer second-mode nodes than actors.
    This required changes in 
    `StatisticCalculator::calculateNetworkEvaluationStatistics`
    and in `siena07internals::getChangeContributionStatistics'.


# RSiena 1.3.7

## 
   
2022-02-18

## Changes in RSiena:  

### Effects:
  * New effect `avDeg`.

# RSiena 1.3.6

## 
   
2022-02-16

## Changes in RSiena:  

### Effects:
  * New effects  `simAllNear`,`simAllFar`, `absOutDiffIntn`, `avDegIntn`.
  * New effects `recipRateInv`, `recipRateLog` (Steffen Triebel).
  * Default internal effect parameter for `outOutActIntn`, `outOutAvIntn`, 
    and `both` changed from 2 to 1.

### Improvements of functionality:  
   * Function `includeInteraction` now also can modify the `initialValue`
     of an effect; and the order of parameters for this function was changed,
     bringing it in line with `setEffect`.
   * Small clarifications of help pages for `includeInteraction` and
     `setEffect`.
   
# RSiena 1.3.5

##

2021-12-15

## Changes in RSiena:  

### Bug corrections  
   * Corrected the check for effects in `initializeFRAN.r` which led to errors
     if interaction effects are included, because of the changes to
     `includeInteraction` in version 1.3.4.

# RSiena 1.3.4

##

2021-12-08

## Changes in RSiena:  

### Effects:
  * New effects `inRateInv`, `inRateLog` (Steffen Triebel).

### Improvements of functionality:  
   * When an effects object with interaction effects is printed,
     the names of the interacting effects are mentioned,
     and prefixes "int." and "i3." were dropped.
   * The check of whether an interaction effect is allowed now is done
     immediately when creating the interaction effect instead of waiting
     for its use in `siena07`. 
   * For function `sienaGroupCreate` some changes were made:   
     if it is applied to a list of length 1, attributes of the single group
     are not recomputed;   
     if it is applied to a list of length larger than 1, the attributes
     "range" and "range2" of behavioral variables of individual groups are
     computed as the range of the unions of the ranges of all the groups.
     The same is done for covariates.
   * `print.sienaGroup` slightly extended.
   * Creation of covariates gives a warning (optional) if all values 
     are missing, and also if all non-missing values are the same.

### Bug corrections  
   * if a `sienaGroup` object is given to `sienaBayes` and some of the covariates
     are constant in one or more of the groups, the ``simX` effect will 
     not run into an error any more; this is achieved by the first change
     mentioned above for `sienaGroupCreate`.
   * `sienaFitThetaTable` in `sienaprint.r` was corrected for `from sienaBayes`.
   * `siena.table` was corrected for `sienaBayes` objects, and this possibility
     was mentioned in the help file.

### Coding changes:
   * In `sienaprint.r`, methods and functions relating to `sienaBayes` omitted.

### Tests:   
   * Test 18, a test for ``sienaGroupCreate`, added to `parallel.R`.
     

# RSiena 1.3.3
   
## 
   
2021-10-09

## Changes in RSiena:  

### Effects:
   * internal effect parameter for diffusion rate effects ('at least p').
   * New effect `outOutAvIntn`.
   * `outOutActIntn` also made available for non-directed explanatory networks 
     and for two-mode dependent networks.

### Improvements of functionality: 
   * `toggleProbabilities` added to output of `sienaRI`.
   * Trial values of `theta` used during Phase 2 of `siena07` added 
     to the `sienaFit` object `ans` as `ans$thetas`.
   * warnings if a data object contains only missings, or only the same value.

### Bug corrections  
   * if a data set contains a constant covariate, the simX effect will
     not run into an error any more; this is obtained by defining the `range` 
     attribute of the covariate as 1 to prevent division by zero,
     and the `simMean` attribute as 0 instead of `NaN` (`sienaDataCreate`). 
     This is relevant especially for `sienaGroup` data sets, where covariates 
     might be constant for some of the groups.`

### Corrections of documentation: 
   * In `sienaAlgorithmCreate`, the default values of `diagonalize` for MoM
     is 0.2; this was corrected in the help file.
     

# RSiena 1.3.2
   
## 
   
2021-07-29

## Changes in RSiena:  

### Improvements of functionality: 
   * Effects of type `creation` or `endow` represented in siena.table
     by `creation` and `maintenance`, respectively.    
     This was erroneously omitted in 1.3.1.

# RSiena 1.3.1
   
## 
   
2021-07-27

### Effects:
   * New effects: `crprodInActIntn` (Nynke Niezink), `XXW`.

### Improvements of functionality: 
   * Effects of type `creation` or `endow` represented in siena.table
     by `creation` and `maintenance`, respectively.
   * `updateTheta` also accepts `sienaBayesFit` objects as `prevAns`.

## Small corrections:    
   * If `upOnly` or `downOnly`, the (out)degree (density) effect is also 
     excluded for symmetric networks 
     (this was reported in `print01Report`, but not carried out).
     This happens in `effects.r`.
   * Message corrected in `sienaDataCreate` if there is an attribute `higher`.
   
   
# RSiena 1.3.0
      
## 
      
2021-05-02  

   * Drop `testthat` in tests.   

# RSiena 1.2.34
   
## 
   
2021-04-30  
   
### New functions:
   * `testSame.RSiena`.
   
### Effects:  
   * New effects: `avInSim` (thanks to Steffen Triebel), `totInSim`, 
     `avInSimPopAlt`, `totInSimPopAlt`, `constant`,
     `avAttHigher`, `avAttLower`, `totAttHigher`, `totAttLower`.
   * Changed effects: endowment and creation types for `avInSim`
     (brought in line with these types for `avSim`).
    
### Improvements of functionality:  
   * `funnelPlot` adapted to lists of `sienaFit` objects
     containing missing estimates or standard errors.
   * `plot.sienaGOF`: new parameter `position`.
   * Small improvements (length of effect names) in `meta.table` and 
    `siena.table`.  
   
### Bug corrections  
   * Restore backward compatibility with respect to checks of `x$gmm`.  
   * In test functions: correct names reported for tested effects by 
     using `ans$requestedEffects` instead of `ans$effects`.
   
### Code improvements   
   * Improved coding of `SimilarityEffect`, using new parts
     of `NetworkDependentBehaviorEffect`.   
   * Changed unsigned actors to int in `Continuousvariable` and 
     `EpochSimulation`;   
     int ...EffectCounts to unsigned in `BehaviorVariable`, 
     to avoid warnings in C++ compilation.
   * Changed name of `similarity(int i, int j)` to `actor_similarity`
     in order to avoid confusion with `similarity(double v1, double v2)`.
   
### Corrections  
   * Took out of `NAMESPACE` a few imported functions from `graphics`, 
     `stats`, `utils` that were not used.  
   * Correction of footer of `CovariateDistance2EgoAltSimNetworkFunction.h`.
  
# RSiena 1.2.33

## 

2021-03-19  


   * Adjusted `configure`, `cleanup` and `Makevars` files for just C++ checks.
   * Pandoc dropped as a system requirement.
 
# RSiena 1.2.32

## 

2021-03-16

### Effects:
   * New effects: homXTransRecTrip, toU.
   * This implied creation of a new effect class dyadANetNetObjective.
   * sqrt versions for parameter 2 for the effects to, toBack, toRecip,
     from, fromMutual.
   * Effects to, toU, toBack, toRecip, MixedInXW are dyadic.
   * Reinstated effect MixedInXW, also with sqrt version for parameter 2.
   * Dropped effect to.2 (identical to `to`) 
     and MixedInWX (identical to `toBack`).

### Improvements of functionality:
   * effectsDocumentation now also includes gmm effects (at the bottom).
   * Improved fromObjectToLaTeX in meta.table and siena.table.
   * Display of deviations from targets changed to after subtraction of targets.
   * Stop if no parameters are estimated and simOnly is FALSE (initializeFRAN).

### Reduction of functionality:
   * Vignette basicRSiena.Rmd dropped (available at website).

### Documentation:
   * Extended description of GMoM in the manual.
   * Description of toBack and toRecip in manual.
   * Changed keyword for some help pages.

### Corrections / safeguards
   * Correction in phase3.2 of a bug that sometimes led to an error message 
     if simOnly.
   * oneModeNet in effects.r: some further cases where the comparison of
     types with `behavior` is replaced by 
     comparison with c(`behavior`, `continuous`).
   * Extra check in phase1.2.
   * Temporarily drop the final part of test16,
     in view of an irreproducible error.

### More neat code:
   * Dropped MixedOutStarFunction, MixedInStarFunction, MixedTwoPathFunction,
     (their functionality replaced by MixedTwoStepFunction).
   * Dropped MixedTwoStepFunction from effects 
     (its place is in effects\generic, and that's were it is).


# RSiena 1.2-31

## 

2021-02-27

   * Generalized method of Moments implemented (Viviana Amati):
     see docs\manual\Changes\_RSiena\_GMoM.tex;
     new function includeGMoMStatistics, extended functionality of siena07.
   * Require R >= 3.5.0.
   * xtable added to `Imports` (used to be in `Suggests`).
   * dyadicCov made to accept also changing dyadic covariates.
   * Used `verbose` condition in sienaGOF also for last console output.
   * new arguments plotAboveThreshold and verbose for funnelPlot.

# RSiena 1.2.30

## 

2021-02-23

   * Resolved issue with continuous dependent behavior variables 
     (Nynke Niezink).

# RSiena 1.2-29

## 

2020-12-10

   * New effects (due to Christoph Stadtfeld):
     transtrip.FR, transtrip.FE, transtrip.EE, WWX.EE, WWX.FR, WXX.FE,
     WXX.ER, XWX.ER, XWX.FE, to.2, toBack, toRecip.
   * New effect transtripX.
   * New functions meta.table and funnelPlot.
   * For effect from.w.ind, option parameter=-1 added.
   * The to effect is an ego effect.
   * New parameter `tested` in sienaGOF.
   * For siena.table, some of the effectNames changed to nice strings,
     so that LaTeX can run without errors if type=`tex`.
   * The object produced by siena08 now has IWLS estimates more easily 
     accessible, as object$muhat and object$se.muhat.
   * Error message in sienaTimeTest for sienaFit objects produced with
     lessMem=TRUE.
   * More extensive error message for error in named vectors in algorithm object
     (checkNames in initializeFRAN).
   * For sienaDataCreate: more extensive error message, and class(...) replaced
     by class(...)[1]. 
   * multiplication factor added to print.sienaAlgorithm if maxlike.
   * In sienaAlgorithmCreate: requirements for mult corrected in help page.
   * In sienaAlgorithmCreate, use the definitions for projname=NULL
     also if any environment variable _R_CHECK* is set. 


# RSiena 1.2-28

## 

2020-09-30

   * Adapted filter `disjoint` so that it operates correctly
     also when the network is symmetric.
     Consequence: constraint that two networks are disjoint
     operates correctly also when one of the networks is symmetric
     and the other is not.
   * Adapted filter `higher` so that it operates correctly
     also when the other network is symmetric.
     Consequence: constraint that one network is at least as high
     as another network operates correctly also when 
     the higher networks is symmetric and the other is not.
   * In `CheckConstraints`, used in `sienaDataCreate`, the requirement 
     was dropped that the two networks have the same symmetry property;
     and for `higher` it is required that if the lower network
     is symmetric, the higher network is also symmetric.
   * In `sienaDataConstraint`, if type is `disjoint` or `atLeastOne`,
     the constraint is also implemented for the pair (net2, net1).
   * Vignette basicRSiena added (was earlier available as a script);
     thanks to James Hollway.

# RSiena 1.2-26

2020-09-17 

   * Changed requirement for tcltk to `Suggests`, 
     and further modified / cleaned up DESCRIPTION.
   * In siena07: if `(!requireNamespace(tcltk))` set batch to TRUE.
   * In NAMESPACE drop tcltk
   * In sienaAlgorithmCreate, use the definitions for projname=NULL
     also if any environment variable _R_CHECK* is set.


