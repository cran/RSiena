/******************************************************************************
 * SIENA: Simulation Investigation for Empirical Network Analysis
 *
 * Web: http://www.stats.ox.ac.uk/~snijders/siena/
 *
 * File: siena07models.h
 *
 * Description: This file contains prototypes for the siena simulation
 * modelling functions called from R
 *****************************************************************************/

#ifndef SIENA07MODELS_H_
#define SIENA07MODELS_H_

/**
 *  Does one forward simulation for all the data by period within group
 */
SEXP model(SEXP DERIV, SEXP DATAPTR, SEXP SEEDS,
	SEXP FROMFINITEDIFF, SEXP MODELPTR, SEXP EFFECTSLIST,
	SEXP THETA, SEXP RANDOMSEED2, SEXP RETURNDEPS, SEXP NEEDSEEDS,
	SEXP USESTREAMS, SEXP ADDCHAINTOSTORE, SEXP NEEDCHANGECONTRIBUTIONS);

/* Does one forward simulation for a specified group and period.
 * Not recommended as it seems slow. Does not currently return chains.
 */
SEXP modelPeriod(SEXP DERIV, SEXP DATAPTR, SEXP SEEDS,
	SEXP FROMFINITEDIFF, SEXP MODELPTR, SEXP EFFECTSLIST,
	SEXP THETA, SEXP RANDOMSEED2, SEXP RETURNDEPS, SEXP NEEDSEEDS,
	SEXP USESTREAMS, SEXP GROUP, SEXP PERIOD);

/** Does some MH steps for a specified group and period.
 * Designed to be used for parallel processing, and currently the only
 * function available. Loop is always constructed in R. Probably would be
 * better to do it in C unless parallel processing.
 */
SEXP mlPeriod(SEXP DERIV, SEXP DATAPTR,
	SEXP MODELPTR, SEXP EFFECTSLIST,
	SEXP THETA, SEXP RETURNDEPS, SEXP GROUP, SEXP PERIOD,
	SEXP NRUNMH);

/* Recalculates the probabilities for a single chain, corresponding to a
 * specific group and period. Optionally returns the scores and derivatives
 * also.
 *
 */
SEXP getChainProbabilitiesList(SEXP CHAIN, SEXP DATAPTR, SEXP MODELPTR,
	SEXP GROUP, SEXP PERIOD, SEXP EFFECTSLIST, SEXP THETA,
	SEXP NEEDSCORES);

/* Calculates the updated probabilities for a set of chains for a single period
 * which have been stored on the model object. Should probably be extended to
 * deal with multiple periods. Uses an EpochSimulation object rather than
 * an MLSimulation object, as needs initialization the same as a forward
 * simulation.
 */
SEXP getStoredChainProbabilities(SEXP DATAPTR, SEXP MODELPTR,
	SEXP GROUP, SEXP PERIOD, SEXP EFFECTSLIST, SEXP THETA);

SEXP clearStoredChains(SEXP MODELPTR);

SEXP getChainProbabilities(SEXP CHAIN, SEXP DATAPTR, SEXP MODELPTR,
	SEXP GROUP, SEXP PERIOD, SEXP EFFECTSLIST, SEXP THETA);

SEXP MCMCcycle(SEXP DATAPTR, SEXP MODELPTR, SEXP MLPTR,
	SEXP EFFECTSLIST, SEXP PERIOD, SEXP GROUP, SEXP SCALEFACTOR,
	SEXP NRUNMH, SEXP NRUNMHBATCHES);

#endif /*SIENA07MODELS_H_*/
