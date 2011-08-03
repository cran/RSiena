/******************************************************************************
 * SIENA: Simulation Investigation for Empirical Network Analysis
 *
 * Web: http://www.stats.ox.ac.uk/~snijders/siena/
 *
 * File: siena07setup.cpp
 *
 * Description: This module contains routines to interface with R,
 * setting up the Data object with data from R. Routines in this file are
 * visible from R.
 *****************************************************************************/
/**
 * @file
 * Sets up the Data object with data from R
 */

#include <vector>
#include <R_ext/Random.h>
#include <Rinternals.h>
#include "siena07internals.h"
#include "siena07utilities.h"
#include "data/Data.h"
#include "data/NetworkLongitudinalData.h"
#include "model/Model.h"
#include "model/ml/Chain.h"
#include "model/State.h"
#include "model/StatisticCalculator.h"
#include "data/ActorSet.h"
#include "model/ml/MLSimulation.h"

using namespace std;
using namespace siena;


extern "C"
{
/**
 *  Creates an array of pointers to Data objects, one for each group
 *  and returns the address of the array to R. Also creates the actor sets
 *  for each group.
 */
SEXP setupData(SEXP OBSERVATIONSLIST, SEXP ACTORSLIST)
{
	/* make error messages go back to R nicely */
	set_terminate(Rterminate);

	int nGroups = length(OBSERVATIONSLIST);

	vector<Data *> *pGroupData = new vector <Data *>;

	for (int group = 0; group < nGroups; group++)
	{
		int observations = INTEGER(VECTOR_ELT(OBSERVATIONSLIST, group))[0];

		pGroupData->push_back(new Data(observations));
		int nNodeSets = length(VECTOR_ELT(ACTORSLIST, group));

		for (int nodeSet = 0; nodeSet < nNodeSets; nodeSet++)
		{
			SEXP nsn;
			PROTECT(nsn = install("nodeSetName"));
			SEXP nodeSetName = getAttrib(VECTOR_ELT(VECTOR_ELT(ACTORSLIST,
						group),
					nodeSet), nsn);
			(*pGroupData)[group]->
				createActorSet(CHAR(STRING_ELT(nodeSetName, 0)),
					length(VECTOR_ELT(VECTOR_ELT(ACTORSLIST,
								group), nodeSet)));
			UNPROTECT(1);
		}
	}
	SEXP RpData;
	RpData = R_MakeExternalPtr((void *) pGroupData, R_NilValue,
		R_NilValue);
	return RpData;
}

/**
 *  Creates all the groups of one mode networks in the data
 *
 */
SEXP OneMode(SEXP RpData, SEXP ONEMODELIST)
{
	/* retrieve the address of our data */
	vector<Data *> * pGroupData = (vector<Data *> *)
		R_ExternalPtrAddr(RpData);

	int nGroups = pGroupData->size();
	/* one mode networks are passed in as list of edgelists with attributes
	   giving the size of the network */
	if (nGroups != length(ONEMODELIST) )
	{
		error ("wrong number of groups");
	}
	for (int group = 0; group < nGroups; group++)
	{
		setupOneModeGroup(VECTOR_ELT(ONEMODELIST, group),
			(*pGroupData)[group]);
	}
	return R_NilValue;
}
/**
 *  Creates all the groups of bipartite networks in the data
 *
 */
SEXP Bipartite(SEXP RpData, SEXP BIPARTITELIST)
{
/* retrieve the address of our data */
	vector<Data *> * pGroupData = (vector<Data *> *)
		R_ExternalPtrAddr(RpData);
	int nGroups = pGroupData->size();
/* bipartite networks are passed in as list of edgelists with attributes
   giving the size of the network */
	if (nGroups != length(BIPARTITELIST) )
	{
		error ("wrong number of groups");
	}
	for (int group = 0; group < nGroups; group++)
	{
		setupBipartiteGroup(VECTOR_ELT(BIPARTITELIST, group),
			(*pGroupData)[group]);
	}
	return R_NilValue;
}
/**
 *  Creates all the groups of behavior networks in the data
 */
SEXP Behavior(SEXP RpData, SEXP BEHLIST)
{
	vector<Data *> * pGroupData = (vector<Data *> *)
		R_ExternalPtrAddr(RpData);
	int nGroups = pGroupData->size();
/* behavior networks are passed in a list of lists of two matrices,
   one of values, one of missing values (boolean) */
	if (nGroups != length(BEHLIST) )
	{
		error ("wrong number of groups");
	}
	for (int group = 0; group < nGroups; group++)
	{
		setupBehaviorGroup(VECTOR_ELT(BEHLIST, group),
			(*pGroupData)[group]);
	}
	return R_NilValue;
}
/**
 *  Creates all the groups of constant covariates in the data
 */

SEXP ConstantCovariates(SEXP RpData, SEXP COCOVARLIST)
{
	vector<Data *> * pGroupData = (vector<Data *> *)
		R_ExternalPtrAddr(RpData);
	int nGroups = pGroupData->size();
/* constant covariates are passed in as vectors with embedded missing values */
/* ignore missings for now */
	if (nGroups != length(COCOVARLIST) )
	{
		error ("wrong number of groups");
	}
	for (int group = 0; group < nGroups; group++)
	{
		setupConstantCovariateGroup(VECTOR_ELT(COCOVARLIST,
				group),
			(*pGroupData)[group]);
	}
	return R_NilValue;
}
/**
 *  Creates all the groups of constant covariates in the data
 */

SEXP ChangingCovariates(SEXP RpData, SEXP VARCOVARLIST)
{
	vector<Data *> * pGroupData = (vector<Data *> *)
		R_ExternalPtrAddr(RpData);
	int nGroups = pGroupData->size();
/* changing covariates are passed in as matrices with embedded missing values */
/* ignore missings for now */
	if (nGroups != length(VARCOVARLIST) )
	{
		error ("wrong number of groups");
	}
	for (int group = 0; group < nGroups; group++)
	{
		setupChangingCovariateGroup(VECTOR_ELT(VARCOVARLIST,
				group),
			(*pGroupData)[group]);
	}
	return R_NilValue;
}
/**
 *  Creates all the groups of constant dyadic covariates in the data
 */

SEXP DyadicCovariates(SEXP RpData, SEXP DYADVARLIST)
{
	vector<Data *> * pGroupData = (vector<Data *> *)
		R_ExternalPtrAddr(RpData);
	int nGroups = pGroupData->size();
/* dyadic covariates are passed in as edgelists with embedded missing values */
/* ignore missings for now */
	if (nGroups != length(DYADVARLIST) )
	{
		error ("wrong number of groups");
	}
	for (int group = 0; group < nGroups; group++)
	{
		setupDyadicCovariateGroup(VECTOR_ELT(DYADVARLIST,
				group),
			(*pGroupData)[group]);
	}
	return R_NilValue;
}
/**
 *  Creates all the groups of changing dyadic covariates in the data
 */

SEXP ChangingDyadicCovariates(SEXP RpData, SEXP VARDYADLIST)
{
	vector<Data *> * pGroupData = (vector<Data *> *)
		R_ExternalPtrAddr(RpData);
	int nGroups = pGroupData->size();
/* dyadic covariates are passed in as lists of edgelists
   with embedded missing values */
/* ignore missings for now */
	if (nGroups != length(VARDYADLIST) )
	{
		error ("wrong number of groups");
	}
	for (int group = 0; group < nGroups; group++)
	{
		setupChangingDyadicCovariateGroup(VECTOR_ELT(VARDYADLIST,
				group),
			(*pGroupData)[group]);
	}
	return R_NilValue;
}
/**
 *  Creates all the composition change events in the data
 */
SEXP ExogEvent(SEXP RpData, SEXP EXOGEVENTLIST)
{
	vector<Data *> * pGroupData = (vector<Data *> *)
		R_ExternalPtrAddr(RpData);
	int nGroups = pGroupData->size();

	/* data for each actor set in each group consists of
	   two items: first a list of events: event type, period, actor, time.
	   Secondly a matrix of booleans indicating whether active at the start
	   of the period. Final period exists in the latter but probably is not
	   necessary. */

	if (nGroups != length(EXOGEVENTLIST) )
	{
		error ("wrong number of groups");
	}
	for (int group = 0; group < nGroups; group++)
	{
		setupExogenousEventGroup(VECTOR_ELT(EXOGEVENTLIST, group),
			(*pGroupData)[group]);
	}
	return R_NilValue;

}

/**
 *  Sets the pairwise constraints for the data
 */
SEXP Constraints(SEXP RpData, SEXP FROMHIGHERLIST, SEXP TOHIGHERLIST,
	SEXP FROMDISJOINTLIST, SEXP TODISJOINTLIST,
	SEXP FROMATLEASTONELIST, SEXP TOATLEASTONELIST)
{
	vector<Data *> * pGroupData = (vector<Data *> *)
		R_ExternalPtrAddr(RpData);
	int nGroups = pGroupData->size();

	for (int group = 0; group < nGroups; group++)
	{
		Data * pData = (*pGroupData)[group];

		/* higher */
		for (int i = 0; i < length(FROMHIGHERLIST); i++)
		{
			pData->
				addNetworkConstraint(CHAR(STRING_ELT(FROMHIGHERLIST, i)),
					CHAR(STRING_ELT(TOHIGHERLIST, i)), HIGHER);
		}
		/* disjoint */
		for (int i = 0; i < length(FROMDISJOINTLIST); i++)
		{
			pData->
				addNetworkConstraint(CHAR(STRING_ELT(FROMDISJOINTLIST, i)),
					CHAR(STRING_ELT(TODISJOINTLIST, i)), DISJOINT);
		}
		/* at least one */
		for (int i = 0; i < length(FROMATLEASTONELIST); i++)
		{
			pData->
				addNetworkConstraint(CHAR(STRING_ELT(FROMATLEASTONELIST,
							i)),
					CHAR(STRING_ELT(TOATLEASTONELIST, i)), AT_LEAST_ONE);
		}
	}
	return R_NilValue;
}


/**
 *  creates the requested basic effects
 */

SEXP effects(SEXP RpData, SEXP EFFECTSLIST)
{
	vector<Data *> * pGroupData = (vector<Data *> *)
		R_ExternalPtrAddr(RpData);

	Model * pModel = new Model();

	/* find the number of columns of the data frame (all will be the same
	   as they are split in R just before the call) */
	//	int n = length(VECTOR_ELT(EFFECTSLIST, 0));
	// get the column names from the names attribute
	SEXP cols;
	PROTECT(cols = install("names"));
	SEXP Names = getAttrib(VECTOR_ELT(EFFECTSLIST, 0), cols);

	int netTypeCol; /* net type */
	int nameCol; /* network name */
	int effectCol;  /* short name of effect */
	int parmCol;
	int int1Col;
	int int2Col;
	int initValCol;
	int typeCol;
	int groupCol;
	int periodCol;
	int pointerCol;
	int rateTypeCol;
	int intptr1Col;
	int intptr2Col;
	int intptr3Col;

	getColNos(Names, &netTypeCol, &nameCol, &effectCol,
		&parmCol, &int1Col, &int2Col, &initValCol,
		&typeCol, &groupCol, &periodCol, &pointerCol,
		&rateTypeCol, &intptr1Col, &intptr2Col, &intptr3Col);

	/* create a structure for the return values */
	SEXP pointers;
	PROTECT(pointers = allocVector(VECSXP, length(EFFECTSLIST)));

	/* loop over the different dependent variables */
	for (int i = 0; i < length(EFFECTSLIST); i++)
	{
		const char * networkName =  CHAR(STRING_ELT(
				VECTOR_ELT(VECTOR_ELT(
						EFFECTSLIST, i),
					nameCol), 0));

		SEXP ptrs =
			createEffects(VECTOR_ELT(EFFECTSLIST, i), pModel, pGroupData,
				networkName, effectCol, parmCol, int1Col,
				int2Col, initValCol, typeCol, groupCol,
				periodCol, pointerCol, rateTypeCol, netTypeCol);

		SET_VECTOR_ELT(pointers, i, ptrs);

	}
	SEXP RpModel;
	PROTECT (RpModel = allocVector(VECSXP, 1));
	SET_VECTOR_ELT(RpModel, 0, R_MakeExternalPtr((void *) pModel,
			R_NilValue,
			R_NilValue));


	/* ans will be the return value */
	SEXP ans;
	PROTECT(ans = allocVector(VECSXP, 2));

	SET_VECTOR_ELT(ans, 1, pointers);
	SET_VECTOR_ELT(ans, 0, RpModel);

	UNPROTECT(4);

	return ans;
}
/**
 *  creates the requested interaction effects
 */

SEXP interactionEffects(SEXP RpData, SEXP RpModel, SEXP EFFECTSLIST)
{
	vector<Data *> * pGroupData = (vector<Data *> *)
		R_ExternalPtrAddr(RpData);
	Model * pModel = (Model *) R_ExternalPtrAddr(RpModel);

	// get the column names from the names attribute

	SEXP cols;
	PROTECT(cols = install("names"));
	SEXP Names = getAttrib(VECTOR_ELT(EFFECTSLIST, 0), cols);

	int netTypeCol; /* net type */
	int nameCol; /* network name */
	int effectCol;  /* short name of effect */
	int parmCol;
	int int1Col;
	int int2Col;
	int initValCol;
	int typeCol;
	int groupCol;
	int periodCol;
	int pointerCol;
	int rateTypeCol;
	int intptr1Col;
	int intptr2Col;
	int intptr3Col;

	getColNos(Names, &netTypeCol, &nameCol, &effectCol,
		&parmCol, &int1Col, &int2Col, &initValCol,
		&typeCol, &groupCol, &periodCol, &pointerCol,
		&rateTypeCol, &intptr1Col, &intptr2Col, &intptr3Col);

	/* create a structure for the return values */
	SEXP pointers;
	PROTECT(pointers = allocVector(VECSXP, length(EFFECTSLIST)));

	/* loop over the different dependent variables */
	for (int i = 0; i < length(EFFECTSLIST); i++)
	{
		if (length(VECTOR_ELT(VECTOR_ELT(EFFECTSLIST, i), 0)) > 0)
		{
			const char * networkName =
				CHAR(STRING_ELT(VECTOR_ELT(VECTOR_ELT(EFFECTSLIST, i),
							nameCol), 0));

			SEXP ptrs =
				createInteractionEffects(VECTOR_ELT(EFFECTSLIST, i),
					pModel, pGroupData,
					networkName, effectCol, parmCol, int1Col,
					int2Col, initValCol, typeCol, groupCol,
					periodCol, pointerCol, rateTypeCol, netTypeCol,
					intptr1Col, intptr2Col, intptr3Col);

			SET_VECTOR_ELT(pointers, i, ptrs);
		}
		else
		{
			SET_VECTOR_ELT(pointers, i,
				R_MakeExternalPtr((void *) 0,
					R_NilValue, R_NilValue));
		}
	}
	/* ans will be the return value */
	SEXP ans;
	PROTECT(ans = allocVector(VECSXP, 1));

	SET_VECTOR_ELT(ans, 0, pointers);

	UNPROTECT(3);

	return ans;
}
/**
 *  removes the objects created for the data.
 */

SEXP deleteData(SEXP RpData)
{

	vector<Data *> * pGroupData = (vector<Data *> *)
		R_ExternalPtrAddr(RpData);
	vector<Data *>::iterator it = pGroupData->begin();
	while(it != pGroupData->end())
	{
		delete *it;
		pGroupData->erase(it); /* not sure I need to do this */
	}
	//	Rprintf("%d delete\n", pGroupData->size());
	delete pGroupData;
	return R_NilValue;
}
/**
 *  removes the model object.
 */

SEXP deleteModel(SEXP RpModel)
{
	Model * pModel = (Model *) R_ExternalPtrAddr(RpModel);
	delete pModel;
	//	Rprintf("deleteModel\n");
	return R_NilValue;
}

/**
 *  sets up the model options of MAXDEGREE, CONDITIONAL
 */
SEXP setupModelOptions(SEXP DATAPTR, SEXP MODELPTR, SEXP MAXDEGREE,
	SEXP CONDVAR, SEXP CONDTARGETS, SEXP PROFILEDATA, SEXP PARALLELRUN,
	SEXP MODELTYPE)
{
	/* get hold of the data vector */
	vector<Data *> * pGroupData = (vector<Data *> *)
		R_ExternalPtrAddr(DATAPTR);
	int nGroups = pGroupData->size();

	/* get hold of the model object */
	Model * pModel = (Model *) R_ExternalPtrAddr(MODELPTR);

	int totObservations = totalPeriods(*pGroupData);

	pModel->numberOfPeriods(totObservations);

	if (!isNull(CONDVAR))
	{
		int *change = INTEGER(CONDTARGETS);
		pModel->conditional(true);
		pModel->conditionalDependentVariable(CHAR(STRING_ELT(CONDVAR,0)));

		int i = 0;

		for (int group = 0; group < nGroups; group++)
		{
			Data * pData = (*pGroupData)[group];

			for (int period = 0;
				 period < pData->observationCount() - 1;
				 period++)
			{
				pModel->targetChange(pData, period, change[i]);
				i++;
			}
		}
	}
	/* get names vector for max degree */
	if (!isNull(MAXDEGREE))
	{
		SEXP Names = getAttrib(MAXDEGREE, R_NamesSymbol);

		for (int group = 0; group < nGroups; group++)
		{
			for (int i = 0; i < length(Names); i++)
			{
				Data * pData = (*pGroupData)[group];
				NetworkLongitudinalData * pNetworkData =
					pData->pNetworkData(CHAR(STRING_ELT(Names, i)));
				pNetworkData->maxDegree(INTEGER(MAXDEGREE)[i]);
			}
		}
	}
	/* set the parallel run flag on the model */
	if (!isNull(PARALLELRUN))
	{
		pModel->parallelRun(true);
	}
	if (!isNull(MODELTYPE))
	{
		pModel->modelType(asInteger(MODELTYPE));
	}
	// print out Data for profiling
	if (asInteger(PROFILEDATA))
	{
		printOutData((*pGroupData)[0]);
	}
	return R_NilValue;

}
/**
 *  Gets target values relative to the input data
 */

SEXP getTargets(SEXP DATAPTR, SEXP MODELPTR, SEXP EFFECTSLIST)
{
	/* get hold of the data vector */
	vector<Data *> * pGroupData = (vector<Data *> *)
		R_ExternalPtrAddr(DATAPTR);

	/* get hold of the model object */
	Model * pModel = (Model *) R_ExternalPtrAddr(MODELPTR);

	int nGroups = pGroupData->size();

	int totObservations = totalPeriods(*pGroupData);

	/* find the number of effects over all dependent variables:
	   sum of lengths of first columns:
	   for dimension of return vector */
	int nEffects = 0;
	for (int i = 0; i < length(EFFECTSLIST); i++)
	{
		nEffects += length(VECTOR_ELT(VECTOR_ELT(EFFECTSLIST, i), 0));
	}

	/* fra will contain the simulated statistics and must be initialised
	   to 0. Use rfra to reduce function evaluations. */
	SEXP fra;
	double * rfra;
	PROTECT(fra = allocMatrix(REALSXP, nEffects, totObservations));
	rfra = REAL(fra);
	for (int i = 0; i < length(fra); i++)
	{
		rfra[i] = 0;
	}
	/* find the targets: for each data object separately:
	   add them up on return to R (easier to check!) */
	int periodFromStart = 0;

	for (int group = 0; group < nGroups; group++)
	{
		Data * pData = (*pGroupData)[group];

		for (int period = 0; period < pData->observationCount() - 1;
			 period++)
		{
			periodFromStart++;
			//	EpochSimulation  Simulation(pData, pModel);
			//Simulation.initialize(period + 1);
			State State (pData, period + 1);
			//State State (&Simulation);
			StatisticCalculator Calculator (pData, pModel, &State,
				period);
			vector<double> statistic(nEffects);
			vector<double> score(nEffects); /* not used */

			getStatistics(EFFECTSLIST, &Calculator, period,
				group, pData, (EpochSimulation *) 0,
				&statistic, &score);
			//getStatistics(EFFECTSLIST, &Calculator, period,
			//			group, pData, &Simulation,
			//&statistic, &score);
			//	Rprintf("%f %f \n",statistic[1], statistic[2]);
			/* fill up matrices for  return value list */
			int iii = (periodFromStart - 1) * nEffects;
			for (unsigned effectNo = 0; effectNo < statistic.size();
				 effectNo++)
			{
				rfra[iii + effectNo] = statistic[effectNo];
			}
		}
	}
	UNPROTECT(1);
	return fra;
}
/** Sets up a minimal chain and does pre burnin and burnin.
 * Processes a complete set of data objects, crewating a chain for each
 * period and storing them on the model object.
 */
SEXP mlMakeChains(SEXP DATAPTR, SEXP MODELPTR, SEXP SIMPLERATES,
	SEXP PROBS, SEXP PRMIN, SEXP PRMIB, SEXP MINIMUMPERM,
	SEXP MAXIMUMPERM, SEXP INITIALPERM)
{
	/* get hold of the data vector */
	vector<Data *> * pGroupData = (vector<Data *> *)
		R_ExternalPtrAddr(DATAPTR);
	int nGroups = pGroupData->size();

	/* find total number of periods to process */
	int totObservations = totalPeriods(*pGroupData);

	/* get hold of the model object */
	Model * pModel = (Model *) R_ExternalPtrAddr(MODELPTR);
	// create chain storage
	pModel->setupChainStore(totObservations);

	/* copy permutation lengths to the model */

	pModel->maximumPermutationLength(REAL(MAXIMUMPERM)[0]);
	pModel->minimumPermutationLength(REAL(MINIMUMPERM)[0]);
	pModel->initialPermutationLength(REAL(INITIALPERM)[0]);
	/* set probability flags */
	pModel->insertDiagonalProbability(REAL(PROBS)[0]);
	pModel->cancelDiagonalProbability(REAL(PROBS)[1]);
	pModel->permuteProbability(REAL(PROBS)[2]);
	pModel->insertPermuteProbability(REAL(PROBS)[3]);
	pModel->deletePermuteProbability(REAL(PROBS)[4]);
	pModel->insertRandomMissingProbability(REAL(PROBS)[5]);
	PrintValue(PROBS);
	pModel->deleteRandomMissingProbability(REAL(PROBS)[6]);

	double * prmin = REAL(PRMIN);
	double * prmib = REAL(PRMIB);

	SEXP minimalChains;
	PROTECT(minimalChains = allocVector(VECSXP, totObservations));
	SEXP currentChains;
	PROTECT(currentChains = allocVector(VECSXP, totObservations));

	int periodFromStart = 0;

	for (int group = 0; group < nGroups; group++)
	{
		Data * pData = (*pGroupData)[group];
		int observations = pData->observationCount() - 1;

		/* create the ML simulation object */
		MLSimulation * pMLSimulation = new MLSimulation(pData, pModel);

		pModel->simpleRates(asInteger(SIMPLERATES));
		pMLSimulation->simpleRates(pModel->simpleRates());

		for (int period = 0; period < observations; period ++)
		{
			// store for later on model
			pModel->missingNetworkProbability(prmin[periodFromStart]);
			pModel->missingBehaviorProbability(prmib[periodFromStart]);

			// put ones for this period on simulation object
			pMLSimulation->
				missingNetworkProbability(prmin[periodFromStart]);
			pMLSimulation->
				missingBehaviorProbability(prmib[periodFromStart]);

			/* initialize the chain: this also initializes the data */
			// does not initialise with previous period missing values yet
			pMLSimulation->connect(period);

			SEXP ch;
			PROTECT(ch = getChainDF(*(pMLSimulation->pChain())));
			SET_VECTOR_ELT(minimalChains, periodFromStart, ch);
			UNPROTECT(1);

			/* get the chain up to a reasonable length */

			GetRNGstate();

			pMLSimulation->preburnin();

			/* do some more steps */

			pMLSimulation->setUpProbabilityArray();

			int numSteps = 500;
			for (int i = 0; i < numSteps; i++)
			{
				pMLSimulation->MLStep();
			}

			/* store chain on Model after creating difference vectors */
			Chain * pChain = pMLSimulation->pChain();
			pChain->createInitialStateDifferences();
			pMLSimulation->createEndStateDifferences();
			pModel->chainStore(*pChain, periodFromStart);

			/* return chain as a list. */
 			SEXP ch1;
 			PROTECT(ch1 = getChainList(*pChain,	*pMLSimulation));
 			SET_VECTOR_ELT(currentChains, periodFromStart, ch1);
 			UNPROTECT(1);

			periodFromStart++;
		}
		delete pMLSimulation;
	}

	SEXP ans;
	PROTECT(ans = allocVector(VECSXP, 2));
	SET_VECTOR_ELT(ans, 0, minimalChains);
	SET_VECTOR_ELT(ans, 1, currentChains);
	UNPROTECT(3);
	PutRNGstate();
	return ans;
}

/** Sets up a minimal chain and does pre burnin and burnin.
 * Processes a complete set of data objects, crewating a chain for each
 * period and storing them on the model object.
 */
SEXP mlInitializeSubProcesses(SEXP DATAPTR, SEXP MODELPTR, SEXP SIMPLERATES,
	SEXP PROBS, SEXP PRMIN, SEXP PRMIB, SEXP MINIMUMPERM,
	SEXP MAXIMUMPERM, SEXP INITIALPERM, SEXP CHAINS, SEXP MISSINGCHAINS)
{
	/* get hold of the data vector */
	vector<Data *> * pGroupData = (vector<Data *> *)
		R_ExternalPtrAddr(DATAPTR);

	int nGroups = pGroupData->size();

	/* find total number of periods to process */
	int totObservations = totalPeriods(*pGroupData);

	/* get hold of the model object */
	Model * pModel = (Model *) R_ExternalPtrAddr(MODELPTR);

	// create chain storage
	pModel->setupChainStore(totObservations);

	/* copy permutation lengths to the model */

	pModel->maximumPermutationLength(REAL(MAXIMUMPERM)[0]);
	pModel->minimumPermutationLength(REAL(MINIMUMPERM)[0]);
	pModel->initialPermutationLength(REAL(INITIALPERM)[0]);
	/* set probability flags */
	pModel->insertDiagonalProbability(REAL(PROBS)[0]);
	pModel->cancelDiagonalProbability(REAL(PROBS)[1]);
	pModel->permuteProbability(REAL(PROBS)[2]);
	pModel->insertPermuteProbability(REAL(PROBS)[3]);
	pModel->deletePermuteProbability(REAL(PROBS)[4]);
	pModel->insertRandomMissingProbability(REAL(PROBS)[5]);
	pModel->deleteRandomMissingProbability(REAL(PROBS)[6]);

	double * prmin = REAL(PRMIN);
	double * prmib = REAL(PRMIB);

	int periodFromStart = 0;

	for (int group = 0; group < nGroups; group++)
	{
		Data * pData = (*pGroupData)[group];
		int observations = pData->observationCount() - 1;

		pModel->simpleRates(asInteger(SIMPLERATES));

		for (int period = 0; period < observations; period ++)
		{
			// store for later on model
			pModel->missingNetworkProbability(prmin[periodFromStart]);
			pModel->missingBehaviorProbability(prmib[periodFromStart]);

			/* copy the chain for this period onto the model */
			Chain * pChain = makeChainFromList(pData,
				VECTOR_ELT(CHAINS, periodFromStart), period);
			pModel->chainStore(*pChain, periodFromStart);

			periodFromStart++;
		}
	}

	return R_NilValue;
}
}

