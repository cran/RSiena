/******************************************************************************
 * SIENA: Simulation Investigation for Empirical Network Analysis
 *
 * Web: http://www.stats.ox.ac.uk/~snijders/siena/
 *
 * File: siena07models.cpp
 *
 * Description: This module contains routines to interface with R,
 * running simulation models. Routines in this file are
 * visible from R.
 *****************************************************************************/
/**
 * @file
 * Runs simulations.
 */

#include <stdexcept>
#include <vector>
#include <cstring>
#include <R_ext/Random.h>
#include <Rinternals.h>
#include "siena07internals.h"
#include "siena07utilities.h"
#include "data/Data.h"
#include "model/Model.h"
#include "model/State.h"
#include "model/StatisticCalculator.h"
#include "utils/Random.h"
#include "model/EpochSimulation.h"
#include "model/variables/BehaviorVariable.h"
#include "model/variables/NetworkVariable.h"
#include "model/ml/MLSimulation.h"
#include "model/ml/Chain.h"
#include "model/ml/NetworkChange.h"
#include "model/ml/BehaviorChange.h"
using namespace std;
using namespace siena;


extern "C"
{

/**
 *  Does one forward simulation for all the data by period within group
 */

SEXP model(SEXP DERIV, SEXP DATAPTR, SEXP SEEDS,
	SEXP FROMFINITEDIFF, SEXP MODELPTR, SEXP EFFECTSLIST,
	SEXP THETA, SEXP RANDOMSEED2, SEXP RETURNDEPS, SEXP NEEDSEEDS,
	SEXP USESTREAMS, SEXP ADDCHAINTOSTORE, SEXP NEEDCHANGECONTRIBUTIONS)
{
	SEXP NEWRANDOMSEED; /* for parallel testing only */
	PROTECT(NEWRANDOMSEED = duplicate(RANDOMSEED2));

	/* create a simulation and return the observed statistics and scores */

	/* get hold of the data vector */
	vector<Data *> * pGroupData = (vector<Data *> *)
		R_ExternalPtrAddr(DATAPTR);

	/* get hold of the model object */
	Model * pModel = (Model *) R_ExternalPtrAddr(MODELPTR);
	//	Rprintf("%x %x\n", pGroupData, pModel);
	int nGroups = pGroupData->size();

	/* find total number of periods to process */
	int totObservations = totalPeriods(*pGroupData);

	int fromFiniteDiff = asInteger(FROMFINITEDIFF);
	int useStreams = asInteger(USESTREAMS);

	int addChainToStore = 0;
	int needChangeContributions = 0;

	int returnDependents = asInteger(RETURNDEPS);

	if (!isNull(ADDCHAINTOSTORE))
	{
		addChainToStore = asInteger(ADDCHAINTOSTORE);
	}

	if (!isNull(NEEDCHANGECONTRIBUTIONS))
	{
		needChangeContributions = asInteger(NEEDCHANGECONTRIBUTIONS);
	}
	int deriv = asInteger(DERIV);
	int needSeeds = asInteger(NEEDSEEDS);

	/* set the deriv flag on the model */
	pModel->needScores(deriv);
	pModel->needChangeContributions((addChainToStore == 1) ||
		(needChangeContributions == 1));

	/* set the chain flag on the model */
	pModel->needChain(returnDependents == 1 || addChainToStore == 1);

	/* update the parameters */
	updateParameters(EFFECTSLIST, THETA, pGroupData, pModel);

	/* ans will be the return value */
	SEXP ans;
	PROTECT(ans = allocVector(VECSXP, 7));

	/* count up the total number of parameters */
	int dim = 0;
	for (int i = 0; i < length(EFFECTSLIST); i++)
	{
		dim += length(VECTOR_ELT(VECTOR_ELT(EFFECTSLIST, i), 0));
	}

	/* get the random seed from R into memory */
	GetRNGstate();

	/* fra will contain the simulated statistics and must be initialised
	   to 0. Use rfra to reduce function evaluations. */
	SEXP fra;
	double * rfra;
	PROTECT(fra = allocMatrix(REALSXP, dim, totObservations));
	rfra = REAL(fra);
	for (int i = 0; i < length(fra); i++)
	{
		rfra[i] = 0;
	}
	/* ntim is the total time taken in each period */
	SEXP ntim;
	double * rntim;
	PROTECT(ntim = allocVector(REALSXP, totObservations));
	rntim = REAL(ntim);
	for (int i = 0; i < length(ntim); i++)
		rntim[i] = 0.0;

	/* sims will be the returned simulated dependent variables */
	SEXP sims;
	PROTECT(sims = allocVector(VECSXP, nGroups));
	if (returnDependents)
	{
		int nVariables = (*pGroupData)[0]->rDependentVariableData().size();
		for (int group = 0; group < nGroups; group++)
		{
			SET_VECTOR_ELT(sims, group,
				allocVector(VECSXP, nVariables));
			for (int variable = 0; variable < nVariables; variable++)
			{
				SET_VECTOR_ELT(VECTOR_ELT(sims, group), variable,
					allocVector(VECSXP, (*pGroupData)[group]->
						observationCount() - 1));
			}
		}
	}
	/* chains will be the returned chains */
	SEXP chains;
	PROTECT(chains = allocVector(VECSXP, nGroups));
	if (returnDependents)
	{
		for (int group = 0; group < nGroups; group++)
		{
			SET_VECTOR_ELT(chains, group,
				allocVector(VECSXP, (*pGroupData)[group]->
					observationCount() - 1));
		}
	}

	/* seed store is a list to save the random states */
	SEXP seedstore;
	PROTECT(seedstore = allocVector(VECSXP, nGroups));
	for (int group = 0; group < nGroups; group++)
	{
		SET_VECTOR_ELT(seedstore, group,
			allocVector(VECSXP, (*pGroupData)[group]->
				observationCount() - 1));
	}

	/* rs will allow us to access or set the .Random.seed in R */
	SEXP rs;
	PROTECT(rs = install(".Random.seed"));

	/* scores will hold the return values of the scores */
	SEXP scores;
	double *rscores;
	PROTECT(scores = allocMatrix(REALSXP, dim, totObservations));
	rscores = REAL(scores);
	for (int i = 0; i < length(scores); i++)
		rscores[i] = 0.0;

	int periodFromStart = 0;

	SEXP Cgstr = R_NilValue;
	SEXP STREAMS = R_NilValue;
	SEXP ans2, ans3, ans4, R_fcall1, R_fcall2, R_fcall3, R_fcall4;
	SEXP seedvector;

	if (useStreams)
	{
		// create an R character string
		PROTECT(Cgstr = allocVector(STRSXP,1));
		SET_STRING_ELT(Cgstr, 0, mkChar("Cg"));

		// find out which stream we are using
		PROTECT(R_fcall1 = lang1(install(".lec.GetStreams")));
		PROTECT(STREAMS = eval(R_fcall1, R_GlobalEnv));
	}


	/* group loop here */
	for (int group = 0; group < nGroups; group++)
	{
		/* random states need store (not fromFiniteDiff)
		   and restore (fromFiniteDiff) for each period
		   within each  group */
		SEXP seeds = R_NilValue;
		if (fromFiniteDiff)
		{
			seeds = VECTOR_ELT(SEEDS, group);
		}

		/* find out how many periods in this Data object */
		Data * pData = (*pGroupData)[group];
		int observations = pData->observationCount();

		/* create my epochsimulation object */
		EpochSimulation * pEpochSimulation  = new
			EpochSimulation(pData, pModel);

		for (int period = 0; period < observations - 1; period++)
		{

			periodFromStart++;

			if (!isNull(RANDOMSEED2)) /* parallel testing versus Siena3 */
			{
				// overwrite R's random number seed
				defineVar(rs, RANDOMSEED2, R_GlobalEnv);
				// get it into memory
				GetRNGstate();
				// move on one
				nextDouble();
				// write it back to R
				PutRNGstate();
			}
			else /* normal run */
			{
				if (fromFiniteDiff) /* restore state */
				{
					if (useStreams) /* using lecuyer random numbers */
					{
						// overwrite the current state in R
						PROTECT(R_fcall2 = lang4(install("[[<-"),
								install(".lec.Random.seed.table"), Cgstr,
								VECTOR_ELT(seeds, period)));
						PROTECT(ans2 = eval(R_fcall2, R_GlobalEnv));
						// get the overwritten state into C table
						PROTECT(R_fcall3 =
							lang2(install(".lec.CurrentStream"),
								STREAMS));
						PROTECT(ans3 = eval(R_fcall3, R_GlobalEnv));
						UNPROTECT(4);
					}
					else /* using normal random numbers */
					{
						// overwrite R's current state
						defineVar(rs, VECTOR_ELT(seeds, period),
							R_GlobalEnv);
						// get the value from .Random.seed into memory
						GetRNGstate();
					}
				}
				else /* save state */
				{
					if (needSeeds)
					{
						if (useStreams)
						{
							PROTECT(R_fcall2 =
								lang2(install(".lec.ResetNextSubstream"),
									STREAMS));
							PROTECT(ans2 = eval(R_fcall2, R_GlobalEnv));

							PROTECT(R_fcall3 =
								lang2(install(".lec.CurrentStream"),
									STREAMS));
							PROTECT(ans3 = eval(R_fcall3, R_GlobalEnv));
							// get the relevant current state from R
							PROTECT(R_fcall4 = lang3(install("[["),
									install(".lec.Random.seed.table"),
									Cgstr));
							ans4 = eval(R_fcall4, R_GlobalEnv);
							// value is not kept unless we duplicate it
							PROTECT(seedvector = duplicate(ans4));
							// store the Cg values
							SET_VECTOR_ELT(VECTOR_ELT(seedstore, group),
								period, seedvector);
							UNPROTECT(6);
						}
						else
						{
							PutRNGstate();
							SET_VECTOR_ELT(VECTOR_ELT(seedstore, group),
								period, findVar(rs, R_GlobalEnv));
						}
					}
				}
			}
			/* run the epoch simulation for this period */
			pEpochSimulation->runEpoch(period);
			State State(pEpochSimulation);
			StatisticCalculator Calculator(pData, pModel, &State,
				period);
			vector<double> statistic(dim);
			vector<double> score(dim);
			getStatistics(EFFECTSLIST, &Calculator,
				period, group, pData, pEpochSimulation,
				&statistic, &score);
			/* fill up matrices for  return value list */
			int iii = (periodFromStart - 1) * dim;
			for (unsigned effectNo = 0; effectNo < statistic.size();
				 effectNo++)
			{
				rfra[iii + effectNo] = statistic[effectNo];

				rscores[iii + effectNo] = score[effectNo];
			}
			if (pModel->conditional())
			{
				rntim[periodFromStart - 1] = pEpochSimulation->time();
			}
			// get simulated network
			if (returnDependents)
			{
				const vector<DependentVariable *> rVariables =
					pEpochSimulation->rVariables();
				for (unsigned i = 0; i < rVariables.size(); i++)
				{
					NetworkVariable * pNetworkVariable =
						dynamic_cast<NetworkVariable *>(rVariables[i]);
					BehaviorVariable * pBehaviorVariable =
						dynamic_cast<BehaviorVariable *>(rVariables[i]);

					if (pNetworkVariable)
					{
						const Network * pNetwork =
							pNetworkVariable->pNetwork();
						SEXP thisEdge = getEdgeList(*pNetwork);
						SET_VECTOR_ELT(VECTOR_ELT(VECTOR_ELT(sims, group),
								i), period, thisEdge);
					}
					else if (pBehaviorVariable)
					{
						SEXP theseValues =
							getBehaviorValues(*pBehaviorVariable);
						SET_VECTOR_ELT(VECTOR_ELT(VECTOR_ELT(sims,
									group), i), period, theseValues);
					}
					else
					{
						throw domain_error(
							"Unexpected class of dependent variable");
					}
				}
				pModel->needChangeContributions(addChainToStore == 0 &&
					needChangeContributions == 1);
				SEXP thisChain =
					getChainDF(*(pEpochSimulation->pChain()));
/*,
								   *pEpochSimulation);*/

				SET_VECTOR_ELT(VECTOR_ELT(chains, group), period,
					thisChain);
				pModel->needChangeContributions(addChainToStore == 1 ||
					needChangeContributions == 1);

			}
			if (addChainToStore)
			{
				pModel->chainStore(*(pEpochSimulation->pChain()),
					periodFromStart - 1);
			}
		} /* end of period */
		delete pEpochSimulation;
	} /* end of group */

	/* send the .Random.seed back to R */
	PutRNGstate();
	NEWRANDOMSEED = findVar(rs, R_GlobalEnv);

	/* set up the return object */
	if (!fromFiniteDiff)
	{
		if (needSeeds)
		{
			SET_VECTOR_ELT(ans, 2, seedstore);
		}
	}
	if (deriv)
	{
		SET_VECTOR_ELT(ans, 1, scores);
	}
	if (returnDependents)
	{
		SET_VECTOR_ELT(ans, 5, sims);/* not done in phase 2 !!test this */
	}
	SET_VECTOR_ELT(ans, 0, fra);
	SET_VECTOR_ELT(ans, 3, ntim);

	if (!isNull(RANDOMSEED2))
	{
		SET_VECTOR_ELT(ans, 4, NEWRANDOMSEED);
	}
	if (useStreams)
	{
		UNPROTECT(3);
	}
	SET_VECTOR_ELT(ans, 6, chains);
	UNPROTECT(9);
	return(ans);
}

/** Does one forward simulation for a specified group and period.
 *  Not recommended as it seems slow. Does not currently return chains.
 */
SEXP modelPeriod(SEXP DERIV, SEXP DATAPTR, SEXP SEEDS,
	SEXP FROMFINITEDIFF, SEXP MODELPTR, SEXP EFFECTSLIST,
	SEXP THETA, SEXP RANDOMSEED2, SEXP RETURNDEPS, SEXP NEEDSEEDS,
	SEXP USESTREAMS, SEXP GROUP, SEXP PERIOD)
{
	/* create a simulation and return the observed statistics and scores */

	/* get hold of the data vector */
	vector<Data *> * pGroupData = (vector<Data *> *)
		R_ExternalPtrAddr(DATAPTR);
	int group = asInteger(GROUP) - 1;

	int period = asInteger(PERIOD) - 1;

	Data * pData = (*pGroupData)[group];

	/* get hold of the model object */
	Model * pModel = (Model *) R_ExternalPtrAddr(MODELPTR);

	int fromFiniteDiff = asInteger(FROMFINITEDIFF);
	int useStreams = asInteger(USESTREAMS); /* always true */
	if (!useStreams)
	{
		error("function modelPeriod called with useStreams FALSE");
	}

	int returnDependents = asInteger(RETURNDEPS);

	int deriv = asInteger(DERIV);
	int needSeeds = asInteger(NEEDSEEDS);

	/* set the deriv flag on the model */
	pModel->needScores(deriv);

	/* set the chain flag on the model */
	pModel->needChain(returnDependents);

	/* update the parameters */
	updateParameters(EFFECTSLIST, THETA, pGroupData, pModel);

	/* count up the total number of parameters */
	int dim = 0;
	for (int i = 0; i < length(EFFECTSLIST); i++)
	{
		dim += length(VECTOR_ELT(VECTOR_ELT(EFFECTSLIST, i), 0));
	}

	/* fra will contain the simulated statistics and must be initialised
	   to 0. Use rfra to reduce function evaluations. */
	SEXP fra;
	double * rfra;
	PROTECT(fra = allocVector(REALSXP, dim));
	rfra = REAL(fra);
	for (int i = 0; i < length(fra); i++)
	{
		rfra[i] = 0;
	}

	/* ntim is the total time taken in this period */
	SEXP ntim;
	double * rntim;
	PROTECT(ntim = allocVector(REALSXP, 1));
	rntim = REAL(ntim);
	rntim[0] = 0.0;

	/* ans will be the return value */
	SEXP ans;
	PROTECT(ans = allocVector(VECSXP, 7));

	/* sims will be the returned simulated dependent variables */
	SEXP sims;
	int nVariables = (*pGroupData)[0]->rDependentVariableData().size();
	PROTECT(sims = allocVector(VECSXP, nVariables));

	/* seed store is a list to save the random state */
	SEXP seedstore;
	PROTECT(seedstore = allocVector(VECSXP, 1));

	/* scores will hold the return values of the scores */
	SEXP scores;
	double *rscores;
	PROTECT(scores = allocVector(REALSXP, dim));
	rscores = REAL(scores);
	for (int i = 0; i < length(scores); i++)
		rscores[i] = 0.0;

	/* random states need store (not fromFiniteDiff)
	   and restore (fromFiniteDiff) for each period
	   within each  group */

	/* create my epochsimulation object */
	EpochSimulation * pEpochSimulation  = new
		EpochSimulation(pData, pModel);

	SEXP Cgstr, ans2, ans3, ans4, STREAMS, R_fcall1, R_fcall2,
		R_fcall3, R_fcall4;

	// create an R character string
	PROTECT(Cgstr = allocVector(STRSXP,1));
	SET_STRING_ELT(Cgstr, 0, mkChar("Cg"));

	// find out which stream we are using
	PROTECT(R_fcall1 = lang1(install(".lec.GetStreams")));
	PROTECT(STREAMS = eval(R_fcall1, R_GlobalEnv));

	if (!isNull(RANDOMSEED2))
	{
		error("non null randomseed2");
	}
	else
	{
		if (fromFiniteDiff) /* restore state */
		{
			// overwrite the current state in R
			PROTECT(R_fcall2 = lang4(install("[[<-"),
					install(".lec.Random.seed.table"), Cgstr,
					VECTOR_ELT(SEEDS, 0)));
			PROTECT(ans2 = eval(R_fcall2, R_GlobalEnv));
			// get the overwritten state into C table
			PROTECT(R_fcall3 = lang2(install(".lec.CurrentStream"),
					STREAMS));
			PROTECT(ans3 = eval(R_fcall3, R_GlobalEnv));
			UNPROTECT(4);
		}
		else /* save state */
		{
			if (needSeeds)
			{
				// move on to next substream in R copy of our stream
				PROTECT(R_fcall2 =
					lang2(install(".lec.ResetNextSubstream"),
						STREAMS));
				PROTECT(ans2 = eval(R_fcall2, R_GlobalEnv));
				// now make this stream the current one so C table
				// contains these values
				PROTECT(R_fcall3 = lang2(install(".lec.CurrentStream"),
						STREAMS));
				PROTECT(ans3 = eval(R_fcall3, R_GlobalEnv));
				// get the relevant current state from R
				PROTECT(R_fcall4 = lang3(install("[["),
						install(".lec.Random.seed.table"), Cgstr));
				PROTECT(ans4 = eval(R_fcall4, R_GlobalEnv));
				// store the Cg values
				SET_VECTOR_ELT(seedstore, 0, ans4);
				UNPROTECT(6);
			}
		}
	}

	/* run the epoch simulation for this period */
	pEpochSimulation->runEpoch(period);

	State State(pEpochSimulation);
	StatisticCalculator Calculator(pData, pModel, &State,
		period);
	vector<double> statistic(dim);
	vector<double> score(dim);
	getStatistics(EFFECTSLIST, &Calculator,
		period, group, pData, pEpochSimulation,
		&statistic, &score);
	/* fill up vector for  return value list */
	for (unsigned effectNo = 0; effectNo < statistic.size();
		 effectNo++)
	{
		rfra[effectNo] = statistic[effectNo];

		rscores[effectNo] = score[effectNo];
	}
	if (pModel->conditional())
	{
		rntim[0] = pEpochSimulation->time();
	}
	// get simulated network
	if (returnDependents)
	{
		const vector<DependentVariable *> rVariables =
			pEpochSimulation->rVariables();
		for (unsigned i = 0; i < rVariables.size(); i++)
		{
			NetworkVariable * pNetworkVariable =
				dynamic_cast<NetworkVariable *>(rVariables[i]);
			BehaviorVariable * pBehaviorVariable =
				dynamic_cast<BehaviorVariable *>(rVariables[i]);

			if (pNetworkVariable)
			{
				const Network * pNetwork =
					pNetworkVariable->pNetwork();
				SEXP thisEdge = getEdgeList(*pNetwork);
				SET_VECTOR_ELT(sims, i, thisEdge);
			}
			else if (pBehaviorVariable)
			{
				SEXP theseValues =
					getBehaviorValues(*pBehaviorVariable);
				SET_VECTOR_ELT(sims, i, theseValues);
			}
			else
			{
				throw domain_error("Unexpected class of dependent variable");
			}
		}
	}
	delete pEpochSimulation;

	/* set up the return object */
	if (!fromFiniteDiff)
	{
		if (needSeeds)
		{
			SET_VECTOR_ELT(ans, 2, seedstore);
		}
	}
	if (deriv)
	{
		SET_VECTOR_ELT(ans, 1, scores);
	}
	if (returnDependents)
	{
		SET_VECTOR_ELT(ans, 5, sims);/* not done in phase 2 !!!!test this*/
	}
	SET_VECTOR_ELT(ans, 0, fra);
	SET_VECTOR_ELT(ans, 3, ntim);

	UNPROTECT(9);
	return(ans);
}

/** Does some MH steps for a specified group and period.
 * Designed to be used for parallel processing, and currently the only
 * function available. Loop is always constructed in R. Probably would be
 * better to do it in C unless parallel processing.
 */
SEXP mlPeriod(SEXP DERIV, SEXP DATAPTR, SEXP MODELPTR, SEXP EFFECTSLIST,
	SEXP THETA, SEXP RETURNDEPS, SEXP GROUP, SEXP PERIOD,
	SEXP NRUNMH, SEXP ADDCHAINTOSTORE, SEXP NEEDCHANGECONTRIBUTIONS)
{
	/* do some MH steps and return the scores and derivs of the chain
	   at the end */

	/* get hold of the data vector */
	vector<Data *> * pGroupData = (vector<Data *> *)
		R_ExternalPtrAddr(DATAPTR);

	int group = asInteger(GROUP) - 1;

	int period = asInteger(PERIOD) - 1;

	int groupPeriod = periodFromStart(*pGroupData, group, period);

	Data * pData = (*pGroupData)[group];

	/* get hold of the model object */
	Model * pModel = (Model *) R_ExternalPtrAddr(MODELPTR);

	/* create the ML simulation object */
	MLSimulation * pMLSimulation = new MLSimulation(pData, pModel);

	pMLSimulation->simpleRates(pModel->simpleRates());

	// next calls are ambiguous unless I use a const pModel
	const Model * pConstModel = pModel;

	pMLSimulation->
		missingNetworkProbability(pConstModel->
			missingNetworkProbability(groupPeriod));
	pMLSimulation->
		missingBehaviorProbability(pConstModel->
			missingBehaviorProbability(groupPeriod));

	// copy chain  back for this period
	Chain * pChain = pModel->rChainStore(groupPeriod).back();
	pMLSimulation->pChain(pChain->copyChain());

	// ready to recreate after the simulation
	pModel->deleteLastChainStore(groupPeriod);

	int addChainToStore = 0;
	int needChangeContributions = 0;

	int returnDependents = 1; //asInteger(RETURNDEPS);

	if (!isNull(ADDCHAINTOSTORE))
	{
		addChainToStore = asInteger(ADDCHAINTOSTORE);
	}

	if (!isNull(NEEDCHANGECONTRIBUTIONS))
	{
		needChangeContributions = asInteger(NEEDCHANGECONTRIBUTIONS);
	}

	int deriv = asInteger(DERIV);

	pModel->needChangeContributions((addChainToStore == 1) ||
		(needChangeContributions == 1));

	/* set the chain flag on the model */
	pModel->needChain(returnDependents == 1 || addChainToStore == 1);

	/* set the deriv flag on the model */
	pModel->needDerivatives(deriv);

	/* update the parameters */
	updateParameters(EFFECTSLIST, THETA, pGroupData, pModel);

	/* count up the total number of parameters */
	int dim = 0;
	for (int i = 0; i < length(EFFECTSLIST); i++)
	{
		dim += length(VECTOR_ELT(VECTOR_ELT(EFFECTSLIST, i), 0));
	}

	/* fra will contain the scores and must be initialised
	   to 0. Use rfra to reduce function evaluations. */
	SEXP fra;
	double * rfra;
	PROTECT(fra = allocVector(REALSXP, dim));
	rfra = REAL(fra);
	for (int i = 0; i < length(fra); i++)
	{
		rfra[i] = 0;
	}

	/* ans will be the return value */
	SEXP ans;
	PROTECT(ans = allocVector(VECSXP, 9));

	/* sims will be the returned chain */
	SEXP sims;
	PROTECT(sims = allocVector(VECSXP, 1));

	/* rs will allow us to access or set the .Random.seed in R */
	SEXP rs;
	PROTECT(rs = install(".Random.seed"));

	/* dff will hold the return values of the derivatives */
	SEXP dff;
	double *rdff;
	PROTECT(dff = allocVector(REALSXP, dim * dim));
	rdff = REAL(dff);
	for (int i = 0; i < length(dff); i++)
	{
		rdff[i] = 0.0;
	}

	// get the value from .Random.seed into memory
	GetRNGstate();
	pModel->needScores(false);
	pModel->needDerivatives(false);

	int nrunMH = asInteger(NRUNMH);
	pModel->numberMLSteps(nrunMH);

	/* run the epoch simulation for this period */
	pMLSimulation->runEpoch(period);

	/* run through current state of chain and calculate
	   scores and derivatives
	*/

	PutRNGstate();
	pModel->needScores(true);
	pModel->needDerivatives(deriv);

	pMLSimulation->updateProbabilities(pMLSimulation->pChain(),
		pMLSimulation->pChain()->pFirst()->pNext(),
		pMLSimulation->pChain()->pLast()->pPrevious());

	/* collect the scores and derivatives */
	State State(pMLSimulation);
	vector<double> derivs(dim * dim);
	vector<double> score(dim);

	getScores(EFFECTSLIST, 	period, group, pData, pMLSimulation,
		&derivs, &score);

	/* get hold of the statistics for accept and reject */
	SEXP accepts;
	PROTECT(accepts = allocVector(INTSXP, 7));
	SEXP rejects;
	PROTECT(rejects = allocVector(INTSXP, 7));
	int * iaccepts = INTEGER(accepts);
	int * irejects = INTEGER(rejects);

	for (int i = 0; i < 7; i++)
	{
		iaccepts[i] = pMLSimulation->acceptances(i);
		irejects[i] = pMLSimulation->rejections(i);
	}

	/* fill up vectors for  return value list */
	for (unsigned effectNo = 0; effectNo < score.size();
		 effectNo++)
	{
		rfra[effectNo] = score[effectNo];
	}

	for (unsigned ii = 0; ii < derivs.size(); ii++)
	{
		rdff[ii] = derivs[ii];
	}
	// get chain
	//	if (returnDependents)
	//	{
	SEXP theseValues =
		getChainList(*(pMLSimulation->pChain()), *pMLSimulation);

	SET_VECTOR_ELT(sims, 0, theseValues);
	//}

	/* set up the return object */
	if (deriv)
	{
		SET_VECTOR_ELT(ans, 6, dff);
	}
	//if (returnDependents)
	//{
	SET_VECTOR_ELT(ans, 5, sims);/* not done in phase 2 !!!!test this*/
	//}
	SET_VECTOR_ELT(ans, 0, fra);
	SET_VECTOR_ELT(ans, 7, accepts);
	SET_VECTOR_ELT(ans, 8, rejects);

	/* store chain on Model */
	pModel->chainStore(*(pMLSimulation->pChain()), groupPeriod);

	delete pMLSimulation;

	UNPROTECT(7);
	return(ans);
}

/** Recalculates the probabilities for a single chain, corresponding to a
 * specific group and period. Optionally returns the scores and derivatives
 * also.
 *
 */
SEXP getChainProbabilitiesList(SEXP CHAIN, SEXP DATAPTR, SEXP MODELPTR,
	SEXP GROUP, SEXP PERIOD, SEXP EFFECTSLIST, SEXP THETA,
	SEXP NEEDSCORES)
{
	/* need to make sure the parameters have been updated first */
	/* get hold of the data vector */
	vector<Data *> * pGroupData = (vector<Data *> *)
		R_ExternalPtrAddr(DATAPTR);
	int group = asInteger(GROUP) - 1;
	int period = asInteger(PERIOD) - 1;
	Data * pData = (*pGroupData)[group];

	/* get hold of the model object */
	Model * pModel = (Model *) R_ExternalPtrAddr(MODELPTR);

	/* See if we need the scores too */
	int needScores = asInteger(NEEDSCORES);
	pModel->needScores(needScores);

	/* create an ml simulation object */
	MLSimulation * pMLSimulation = new MLSimulation(pData, pModel);

	/* chain is a list: aspect, varname, ego, alter, diff */

	/* create a chain */
	Chain * pChain = new Chain(pData);

	/* set period */
	pChain->period(period);

	for (int i = 0; i < length(CHAIN); i++)
	{
		SEXP THISCHAIN;
		THISCHAIN = VECTOR_ELT(CHAIN, i);
		if (strcmp(CHAR(STRING_ELT(VECTOR_ELT(THISCHAIN, 0), 0)),
				"Network") == 0)
		{
			NetworkChange * pNetworkChange = new NetworkChange
				(pData->pNetworkData(CHAR(STRING_ELT(VECTOR_ELT(THISCHAIN,
								2), 0))),
					asInteger(VECTOR_ELT(THISCHAIN, 3)),
					asInteger(VECTOR_ELT(THISCHAIN, 4)));
			pChain->insertBefore(pNetworkChange, pChain->pLast());
		}
		else
		{
			BehaviorChange * pBehaviorChange = new BehaviorChange
				(pData->pBehaviorData(CHAR(STRING_ELT(VECTOR_ELT(THISCHAIN,
								2), 0))),
					asInteger(VECTOR_ELT(THISCHAIN, 3)),
					asInteger(VECTOR_ELT(THISCHAIN, 5)));
			pChain->insertBefore(pBehaviorChange, pChain->pLast());
		}
	}

	/* update the parameters */
	updateParameters(EFFECTSLIST, THETA, pGroupData, pModel);
	pMLSimulation->updateParameters();

	/* calculate the probabilities: this uses runEpoch so we need to
	   set the number of steps to zero first */
	pModel->numberMLSteps(0);
	pMLSimulation->pChainProbabilities(pChain, period);

	/* count up the total number of parameters */
	int dim = 0;
	for (int i = 0; i < length(EFFECTSLIST); i++)
	{
		dim += length(VECTOR_ELT(VECTOR_ELT(EFFECTSLIST, i), 0));
	}
	/* scores will hold the return values of the scores */
	SEXP scores, dff;
	double *rscores, *rdff;
	PROTECT(scores = allocVector(REALSXP, dim));
	rscores = REAL(scores);
	PROTECT(dff = allocMatrix(REALSXP, dim, dim));
	rdff = REAL(dff);

	/* collect the scores and derivatives */
	State State(pMLSimulation);
	vector<double> derivs(dim * dim);
	vector<double> score(dim);
	getScores(EFFECTSLIST, 	period, group, pData, pMLSimulation,
		&derivs, &score);

	for (unsigned effectNo = 0; effectNo < score.size();
		 effectNo++)
	{
		rscores[effectNo] = score[effectNo];
	}

	/* get the chain with probs */
	SEXP ans;
	PROTECT(ans = getChainList(*(pMLSimulation->pChain()), *pMLSimulation));

	delete pMLSimulation;

	SEXP returnval;
	returnval = allocVector(VECSXP, 2);
	SET_VECTOR_ELT(returnval, 0, ans);
	SET_VECTOR_ELT(returnval, 1, scores);
	UNPROTECT(3);
	return  returnval;
}

/** Calculates the updated probabilities for a set of chains for a single period
 * which have been stored on the model object. Should probably be extended to
 * deal with multiple periods. Uses an EpochSimulation object rather than
 * an MLSimulation object, as needs initialization the same as a forward
 * simulation.
 */
SEXP getStoredChainProbabilities(SEXP DATAPTR, SEXP MODELPTR,
	SEXP GROUP, SEXP PERIOD, SEXP EFFECTSLIST, SEXP THETA)
{
	/* need to make sure the parameters have been updated first */
	/* get hold of the data vector */
	vector<Data *> * pGroupData = (vector<Data *> *)
		R_ExternalPtrAddr(DATAPTR);
	int group = asInteger(GROUP) - 1;
	int period = asInteger(PERIOD) - 1;
	Data * pData = (*pGroupData)[group];

	/* get hold of the model object */
	Model * pModel = (Model *) R_ExternalPtrAddr(MODELPTR);

	/* create a simulation  object */
	EpochSimulation * pEpochSimulation = new
		EpochSimulation(pData, pModel);

	/* initialize for this period */
	pEpochSimulation->initialize(period);

	/* update the parameters */
	updateParameters(EFFECTSLIST, THETA, pGroupData, pModel);
	pEpochSimulation->updateParameters();

	int periodNo = periodFromStart(*pGroupData, group, period);
	unsigned numberChains = (pModel->rChainStore(periodNo)).size();
	SEXP logprob;
	PROTECT(logprob = allocVector(REALSXP, numberChains));
	double * rlogprob = REAL(logprob);

	/* get the probabilities: */

	for (unsigned i = 0; i < numberChains; i++)
 	{
		vector<Chain *> rChainStore = pModel->rChainStore(periodNo);
 		Chain * pChain = rChainStore[i];
 		rlogprob[i] = pEpochSimulation->
 			calculateChainProbabilities(pChain);
 	}
	UNPROTECT(1);
	delete pEpochSimulation;
	return logprob;
}

/** Clears the chains that have been stored on a model. Leave final one for ML.
 */
SEXP clearStoredChains(SEXP MODELPTR, SEXP MAXLIKE)
{
	int maxlike = asInteger(MAXLIKE);
	/* get hold of the model object */
	Model * pModel = (Model *) R_ExternalPtrAddr(MODELPTR);
	if (maxlike > 0)
	{
		pModel->partClearChainStore();
	}
	else
	{
		pModel->clearChainStore();
	}
	return R_NilValue;
}

/** Reads a chain from a data frame into C, updates theta and recalculates
 * the probabilities.
 */
SEXP getChainProbabilities(SEXP CHAIN, SEXP DATAPTR, SEXP MODELPTR,
	SEXP GROUP, SEXP PERIOD, SEXP EFFECTSLIST, SEXP THETA)
{
	/* need to make sure the parameters have been updated first */
	/* get hold of the data vector */
	vector<Data *> * pGroupData = (vector<Data *> *)
		R_ExternalPtrAddr(DATAPTR);
	int group = asInteger(GROUP) - 1;
	int period = asInteger(PERIOD) - 1;
	Data * pData = (*pGroupData)[group];

	/* get hold of the model object */
	Model * pModel = (Model *) R_ExternalPtrAddr(MODELPTR);

	/* create an ml simulation object */
	MLSimulation * pMLSimulation = new MLSimulation(pData, pModel);

	/* chain is a data.frame: aspect, varname, ego, alter, diff */
	SEXP ASPECT = VECTOR_ELT(CHAIN, 0);
	SEXP VARNAME = VECTOR_ELT(CHAIN, 2);
	SEXP EGO = VECTOR_ELT(CHAIN, 3);
	SEXP ALTER = VECTOR_ELT(CHAIN, 4);
	SEXP DIFF = VECTOR_ELT(CHAIN, 5);
	int * ego = INTEGER(EGO);
	int * alter = INTEGER(ALTER);
	int * difference = INTEGER(DIFF);

	/* create a chain */
	Chain * pChain = new Chain(pData);

	/* set period */
	pChain->period(period);

	for (int i = 0; i < length(ASPECT); i++)
	{
		if (strcmp(CHAR(STRING_ELT(ASPECT, i)), "Network") == 0)
		{
			NetworkChange * pNetworkChange = new NetworkChange
				(pData->pNetworkData(CHAR(STRING_ELT(VARNAME, i))),
					ego[i], alter[i]);
			pChain->insertBefore(pNetworkChange, pChain->pLast());
		}
		else
		{
			BehaviorChange * pBehaviorChange = new BehaviorChange
				(pData->pBehaviorData(CHAR(STRING_ELT(VARNAME, i))),
					ego[i], difference[i]);
			pChain->insertBefore(pBehaviorChange, pChain->pLast());
		}
	}

	/* update the parameters */
	updateParameters(EFFECTSLIST, THETA, pGroupData, pModel);

	/* calculate the probabilities: this uses runEpoch so we need to
	   set the number of steps to zero first */
	pModel->numberMLSteps(0);
	pMLSimulation->pChainProbabilities(pChain, period);

	/* get the chain with probs */
	SEXP ans;
	ans = getChainDF(*(pMLSimulation->pChain()));
	delete pMLSimulation;
	return  ans;
}

/** Samples from parameter distributions, then runs MH steps.
 */
SEXP MCMCcycle(SEXP DATAPTR, SEXP MODELPTR,
	SEXP EFFECTSLIST, SEXP PERIOD, SEXP GROUP, SEXP SCALEFACTOR,
	SEXP NRUNMH, SEXP NRUNMHBATCHES)
{
	/* get hold of the data vector */
	vector<Data *> * pGroupData = (vector<Data *> *)
		R_ExternalPtrAddr(DATAPTR);

	int period = asInteger(PERIOD) - 1;
	int group = asInteger(GROUP) - 1;;
	int groupPeriod = periodFromStart(*pGroupData, group, period);
	Data * pData = (*pGroupData)[group];

	/* get hold of the model object */
	Model * pModel = (Model *) R_ExternalPtrAddr(MODELPTR);


	int nrunMH = asInteger(NRUNMH);
	int nrunMHBatches = asInteger(NRUNMHBATCHES);
	double scaleFactor = asReal(SCALEFACTOR);

	pModel->numberMHBatches(nrunMHBatches) ;

	pModel->numberMLSteps(nrunMH);

	pModel->BayesianScaleFactor(scaleFactor);

	/* set up storage for statistics for MH accept and reject */
	SEXP accepts;
	PROTECT(accepts = allocMatrix(INTSXP, nrunMHBatches, 7));
	SEXP rejects;
	PROTECT(rejects = allocMatrix(INTSXP, nrunMHBatches, 7));
	int * iaccepts = INTEGER(accepts);
	int * irejects = INTEGER(rejects);

	GetRNGstate();

	/* create the ML simulation object */
	MLSimulation * pMLSimulation = new MLSimulation(pData, pModel);

	pMLSimulation->simpleRates(pModel->simpleRates());

	// next calls are ambiguous unless I use a const pModel
	const Model * pConstModel = pModel;

	pMLSimulation->
		missingNetworkProbability(pConstModel->
			missingNetworkProbability(groupPeriod));
	pMLSimulation->
		missingBehaviorProbability(pConstModel->
			missingBehaviorProbability(groupPeriod));

	// copy chain  back for this period
	Chain * pChain = pModel->rChainStore(groupPeriod).back();
	pMLSimulation->pChain(pChain->copyChain());

	// ready to recreate after the simulation
	pModel->deleteLastChainStore(groupPeriod);

	pMLSimulation->initializeMCMCcycle();
	for (int batch = 0; batch < nrunMHBatches; batch++)
	{
		pMLSimulation->runEpoch(period);
		// store the MH statistics for this batch
		for (int i = 0; i < 7; i++)
		{
			iaccepts[ i * nrunMHBatches + batch] =
				pMLSimulation->acceptances(i);
			irejects[ i * nrunMHBatches + batch] =
				pMLSimulation->rejections(i);
		}

		pMLSimulation->MHPstep();
	}

	PutRNGstate();
	/* store chain on Model */
	pModel->chainStore(*(pMLSimulation->pChain()), groupPeriod);


	SEXP ans;  // main return list
	PROTECT(ans = allocVector(VECSXP, 5));

	SEXP BayesAccepts; // vector of acceptances
	PROTECT(BayesAccepts = allocVector(INTSXP, nrunMHBatches));
	int * iBayes = INTEGER(BayesAccepts);
	for (int i = 0; i < nrunMHBatches; i++)
	{
		iBayes[i] =  pMLSimulation->BayesAcceptances(i);
	}

	/* count up the total number of parameters */
	int dim = 0;
	for (int i = 0; i < length(EFFECTSLIST); i++)
	{
		dim += length(VECTOR_ELT(VECTOR_ELT(EFFECTSLIST, i), 0));
	}

	SEXP BayesCandidates; // matrix of candidate parameter values
	PROTECT(BayesCandidates = allocMatrix(REALSXP, nrunMHBatches, dim));
	double * rcandidates = REAL(BayesCandidates);

	SEXP BayesShapes;
	PROTECT(BayesShapes = allocMatrix(INTSXP, nrunMHBatches, dim));
	int * ibayesshapes = INTEGER(BayesShapes);

	vector<int> shapes(nrunMHBatches * dim);
	vector<double> candidates(nrunMHBatches * dim);

	getCandidatesAndShapes(EFFECTSLIST, period, group, pData, pMLSimulation,
		&candidates, &shapes, nrunMHBatches);
	/* fill up vectors for  return value list */
	for (unsigned effectNo = 0; effectNo < candidates.size();
		 effectNo++)
	{
		rcandidates[effectNo] = candidates[effectNo];
		ibayesshapes[effectNo] = shapes[effectNo];
	}
	delete pMLSimulation;

	SET_VECTOR_ELT(ans, 0, BayesAccepts);
	SET_VECTOR_ELT(ans, 1, BayesCandidates);
	SET_VECTOR_ELT(ans, 2, BayesShapes);
	SET_VECTOR_ELT(ans, 3, accepts);
	SET_VECTOR_ELT(ans, 4, rejects);
	UNPROTECT(6);
	return  ans;
}


}