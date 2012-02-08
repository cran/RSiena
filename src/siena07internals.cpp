/******************************************************************************
 * SIENA: Simulation Investigation for Empirical Network Analysis
 *
 * Web: http://www.stats.ox.ac.uk/~snijders/siena/
 *
 * File: siena07internals.cpp
 *
 * Description: This module contains routines used in to set up data in C.
 * Internal functions to C++: not visible from R.
 *****************************************************************************/
/**
 * @file
 * Internal routines used to set up the Data object with data from R.
 * Not visible from R.
 */
#include <vector>
#include <cstring>
#include <Rinternals.h>
#include "siena07internals.h"
#include "data/Data.h"
#include "data/LongitudinalData.h"
#include "data/NetworkLongitudinalData.h"
#include "data/OneModeNetworkLongitudinalData.h"
#include "data/BehaviorLongitudinalData.h"
#include "data/ChangingDyadicCovariate.h"
#include "data/ConstantDyadicCovariate.h"
#include "data/ChangingCovariate.h"
#include "data/ConstantCovariate.h"
#include "data/ExogenousEvent.h"
#include "model/Model.h"
#include "model/EffectInfo.h"
#include "model/State.h"
#include "model/StatisticCalculator.h"
#include "data/ActorSet.h"
#include "model/EpochSimulation.h"
#include "model/variables/DependentVariable.h"
#include "model/variables/BehaviorVariable.h"
#include "model/variables/NetworkVariable.h"
#include "model/ml/MLSimulation.h"

using namespace std;
using namespace siena;


/**
 * Matches column names with indices. The object Names is the names of the
 * effects data frame.
 */
void getColNos(SEXP Names, int * netTypeCol, int * nameCol, int * effectCol,
			   int * parmCol, int * int1Col, int * int2Col, int * initValCol,
			   int * typeCol, int * groupCol, int * periodCol, int * pointerCol,
	int * rateTypeCol, int * intptr1Col, int * intptr2Col, int * intptr3Col)
{
	*netTypeCol = -1; /* net type */
	*nameCol = -1; /* network name */
	*effectCol = -1;  /* short name of effect */
	*parmCol = -1;
	*int1Col = -1;
	*int2Col = -1;
	*initValCol = -1;
	*typeCol = -1;
	*groupCol = -1;
	*periodCol = -1;
	*pointerCol = -1;
	*rateTypeCol = -1;
	*intptr1Col = -1;
	*intptr2Col = -1;
	*intptr3Col = -1;

	int n = length(Names);
	for (int j = 0; j < n; j++)
	{
		if (strcmp(CHAR(STRING_ELT(Names, j)), "netType") == 0)
		{
			*netTypeCol = j;
		}
		if (strcmp(CHAR(STRING_ELT(Names, j)), "name") == 0)
		{
			*nameCol = j;
		}
		if (strcmp(CHAR(STRING_ELT(Names, j)), "shortName") == 0)
		{
			*effectCol = j;
		}
		if (strcmp(CHAR(STRING_ELT(Names, j)), "parm") == 0)
		{
			*parmCol = j;
		}

		if (strcmp(CHAR(STRING_ELT(Names, j)), "interaction1") == 0)
		{
			*int1Col = j;
		}
		if (strcmp(CHAR(STRING_ELT(Names, j)), "interaction2") == 0)
		{
			*int2Col = j;
		}
		if (strcmp(CHAR(STRING_ELT(Names, j)), "initialValue") == 0)
		{
			*initValCol = j;
		}
		if (strcmp(CHAR(STRING_ELT(Names, j)), "type") == 0)
		{
			*typeCol = j;
		}
		if (strcmp(CHAR(STRING_ELT(Names, j)), "group") == 0)
		{
			*groupCol = j;
		}
		if (strcmp(CHAR(STRING_ELT(Names, j)), "period") == 0)
		{
			*periodCol = j;
		}
		if (strcmp(CHAR(STRING_ELT(Names, j)), "effectPtr") == 0)
		{
			*pointerCol = j;
		}
		if (strcmp(CHAR(STRING_ELT(Names, j)), "rateType") == 0)
		{
			*rateTypeCol = j;
		}
		if (strcmp(CHAR(STRING_ELT(Names, j)), "effect1") == 0)
		{
			*intptr1Col = j;
		}
		if (strcmp(CHAR(STRING_ELT(Names, j)), "effect2") == 0)
		{
			*intptr2Col = j;
		}
		if (strcmp(CHAR(STRING_ELT(Names, j)), "effect3") == 0)
		{
			*intptr3Col = j;
		}
	}
	if (*netTypeCol < 0)
	{
		error("cannot find nettype");
	}
	if (*nameCol < 0)
	{
		error("cannot find network name");
	}
	if (*effectCol < 0)
	{
		error("cannot find effectName");
	}
	if (*parmCol < 0)
	{
		error("cannot find internal parameter");
	}
	if (*int1Col < 0)
	{
		error("cannot find interaction1");
        }

	if (*int2Col < 0)
	{
		error("cannot find interaction2");
	}
	if (*initValCol < 0)
	{
		error("cannot find initial value");
	}
	if (*groupCol < 0)
	{
		error("cannot find group");
	}
	if (*periodCol < 0)
	{
		error("cannot find period");
	}
	if (*pointerCol < 0)
	{
		error("cannot find effect pointer");
	}
	if (*rateTypeCol < 0)
	{
		error("cannot find rate type");
	}
	if (*intptr1Col < 0)
	{
		error("cannot find effect1");
	}
	if (*intptr2Col < 0)
	{
		error("cannot find effect2");
	}
	if (*intptr3Col < 0)
	{
		error("cannot find effect3");
	}
//Rprintf("%d parmcol\n", *parmCol);
}


/**
 *  updates the parameter values for each of the effects.
 */
void updateParameters(SEXP EFFECTSLIST, SEXP THETA, vector<Data *> *
						  pGroupData, Model * pModel)
{
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

	int thetasub = -1;
	/* find each effect and update its weight */
	for (int net = 0; net < length(EFFECTSLIST); net++)
	{
		const char * networkName =
			CHAR(STRING_ELT(VECTOR_ELT(VECTOR_ELT(EFFECTSLIST, net),
									   nameCol), 0));
		SEXP EFFECTS = VECTOR_ELT(EFFECTSLIST, net);
		for (int eff = 0; eff < length(VECTOR_ELT(EFFECTS,0)); eff++)
		{
			thetasub = thetasub + 1;
			const char * effectName =
				CHAR(STRING_ELT(VECTOR_ELT(EFFECTS, effectCol),  eff));
			double currentValue = REAL(THETA)[thetasub];
			const char * effectType =
				CHAR(STRING_ELT(VECTOR_ELT(EFFECTS, typeCol), eff));
			const char * netType =
				CHAR(STRING_ELT(VECTOR_ELT(EFFECTS, netTypeCol), eff));

			if (strcmp(effectType, "rate") == 0 &&
				strcmp(effectName, "Rate") == 0)
			{
				int group =
					INTEGER(VECTOR_ELT(EFFECTS, groupCol))[eff]  - 1;
				int period =
					INTEGER(VECTOR_ELT(EFFECTS, periodCol))[eff] - 1;
				/* find the network */
				Data * pData = (*pGroupData)[group];

				if (strcmp(netType, "behavior") == 0)
				{
					BehaviorLongitudinalData * pNetwork =
						pData->pBehaviorData(networkName);
					pModel->basicRateParameter(pNetwork,
						period,
						currentValue);
				}
				else
				{
					NetworkLongitudinalData * pNetwork =
						pData->pNetworkData(networkName);
					pModel->basicRateParameter(pNetwork,
						period,
						currentValue);
				}
			}
			else
			{
				EffectInfo * pEffectInfo =
					(EffectInfo *) R_ExternalPtrAddr(
						VECTOR_ELT(VECTOR_ELT(EFFECTS, pointerCol), eff));
				pEffectInfo->parameter(currentValue);
			}
		}
	}

	UNPROTECT(1);
	return;
}


/**
 * Create one observation for a one mode Network: ties, missing, structural
 *
 */
void setupOneModeNetwork(SEXP ONEMODE,
	OneModeNetworkLongitudinalData * pNetworkData,
	int observation)
{
	/* one mode networks are passed in as list of edgelists with attributes
	 giving the size of the network - not checked yet*/

	// Tie values
	//Rprintf("%x\n", pNetworkData);
	SEXP ONEMODEVALS = VECTOR_ELT(ONEMODE, 0);
	int *start = INTEGER(ONEMODEVALS);
	int listlen = ncols(ONEMODEVALS);
	int pos = 0;

	for (int row = 0; row < listlen; row++)
	{
		int i;
		int j;
		int val;
		i = start[pos++];
		j = start[pos++];
		val = start[pos++];
		pNetworkData->tieValue(i - 1, j - 1, observation, val);
	}

	// Missingness

	ONEMODEVALS = VECTOR_ELT(ONEMODE, 1);
	start = INTEGER(ONEMODEVALS);
	listlen = ncols(ONEMODEVALS);
	pos = 0;

	for (int row = 0; row < listlen; row++)
	{
		int i;
		int j;
		int val;
		i = start[pos++];
		j = start[pos++];
		val = start[pos++];
		pNetworkData->missing(i - 1, j - 1, observation, val);
	}

	// Structural ties

	ONEMODEVALS = VECTOR_ELT(ONEMODE, 2);
	start = INTEGER(ONEMODEVALS);
	listlen = ncols(ONEMODEVALS);
	pos = 0;

	for (int row = 0; row < listlen; row++)
	{
		int i;
		int j;
		int val;
		i = start[pos++];
		j = start[pos++];
		val = start[pos++];
		pNetworkData->structural(i - 1, j - 1, observation, val);
	}
}


/**
 * Create all observations for a one mode Network
 *
 */
void setupOneModeObservations(SEXP ONEMODES,
			      OneModeNetworkLongitudinalData *
                              pOneModeNetworkLongitudinalData)

{
    int observations = length(ONEMODES);
    if (observations != pOneModeNetworkLongitudinalData->observationCount())
    {
		error ("wrong number of observations in OneMode");
    }
    SEXP uo;
    PROTECT(uo = install("uponly"));
    SEXP uponly = getAttrib(ONEMODES, uo);
    SEXP dow;
    PROTECT(dow = install("downonly"));
    SEXP downonly = getAttrib(ONEMODES, dow);
    for (int period = 0; period < (observations - 1); period++)
    {
        pOneModeNetworkLongitudinalData->upOnly(period,
                                         LOGICAL(uponly)[period]);
        pOneModeNetworkLongitudinalData->downOnly(period,
                                         LOGICAL(downonly)[period]);
    }
    for (int period = 0; period < observations; period++)
    {
    	setupOneModeNetwork(VECTOR_ELT(ONEMODES, period),
			pOneModeNetworkLongitudinalData,
			period);
	}
    UNPROTECT(2);
}
/**
 * Create one group of one mode Networks
 *
 */
void setupOneModeGroup(SEXP ONEMODEGROUP, Data * pData)
{
    int nOneMode = length(ONEMODEGROUP);

    for (int oneMode = 0; oneMode < nOneMode; oneMode++)
    {
        SEXP as;
        PROTECT(as = install("nodeSet"));
        SEXP actorSet = getAttrib(VECTOR_ELT(ONEMODEGROUP, oneMode), as);
        SEXP symm;
        PROTECT(symm = install("symmetric"));
        SEXP symmetric = getAttrib(VECTOR_ELT(ONEMODEGROUP, oneMode), symm);
		SEXP balm;
        PROTECT(balm = install("balmean"));
        SEXP balmean = getAttrib(VECTOR_ELT(ONEMODEGROUP, oneMode), balm);
		SEXP strm;
        PROTECT(strm = install("structmean"));
        SEXP structmean = getAttrib(VECTOR_ELT(ONEMODEGROUP, oneMode), strm);
		SEXP avin;
        PROTECT(avin = install("averageInDegree"));
        SEXP averageInDegree = getAttrib(VECTOR_ELT(ONEMODEGROUP, oneMode),
			avin);
		SEXP avout;
        PROTECT(avout = install("averageOutDegree"));
        SEXP averageOutDegree = getAttrib(VECTOR_ELT(ONEMODEGROUP, oneMode),
			avout);
		SEXP nm;
        PROTECT(nm = install("name"));
        SEXP name = getAttrib(VECTOR_ELT(ONEMODEGROUP, oneMode), nm);
        const ActorSet * myActorSet = pData->pActorSet(CHAR(STRING_ELT(
					actorSet, 0)));
		OneModeNetworkLongitudinalData *  pOneModeNetworkLongitudinalData =
			pData->createOneModeNetworkData(CHAR(STRING_ELT(name, 0)),
                                     myActorSet);
        pOneModeNetworkLongitudinalData->symmetric(*(LOGICAL(symmetric)));
        pOneModeNetworkLongitudinalData->balanceMean(*(REAL(balmean)));
        pOneModeNetworkLongitudinalData->structuralMean(*(REAL(structmean)));
        pOneModeNetworkLongitudinalData->
			averageInDegree(*(REAL(averageInDegree)));
        pOneModeNetworkLongitudinalData->
			averageOutDegree(*(REAL(averageOutDegree)));
		setupOneModeObservations(VECTOR_ELT(ONEMODEGROUP, oneMode),
			pOneModeNetworkLongitudinalData);
		//	Rprintf("%f %f\n", pOneModeNetworkLongitudinalData->
		//	averageInDegree(),  pOneModeNetworkLongitudinalData->
		//	averageOutDegree());

		// Once all network data has been stored, calculate some
		// statistical properties of that data.
		pOneModeNetworkLongitudinalData->calculateProperties();
		//Rprintf("%f %f\n", pOneModeNetworkLongitudinalData->
		//	averageInDegree(), pOneModeNetworkLongitudinalData->
		//	averageOutDegree());
        UNPROTECT(7);
    }
}

/**
 * Create one observation for a bipartite Network: ties, missing, structural
 *
 */
void setupBipartiteNetwork(SEXP BIPARTITE,
	NetworkLongitudinalData * pNetworkData,
	int observation)
{
	/* bipartite networks are passed in as list of edgelists with attributes
	 giving the size of the network - not checked yet*/

	// Tie values

	SEXP BIPARTITEVALS = VECTOR_ELT(BIPARTITE, 0);
	int *start = INTEGER(BIPARTITEVALS);
	int listlen = ncols(BIPARTITEVALS);
	int pos = 0;

	for (int row = 0; row < listlen; row++)
	{
		int i;
		int j;
		int val;
		i = start[pos++];
		j = start[pos++];
		val = start[pos++];
		pNetworkData->tieValue(i - 1, j - 1, observation, val);
	}

	// Missingness

	BIPARTITEVALS = VECTOR_ELT(BIPARTITE, 1);
	start = INTEGER(BIPARTITEVALS);
	listlen = ncols(BIPARTITEVALS);
	pos = 0;

	for (int row = 0; row < listlen; row++)
	{
		int i;
		int j;
		int val;
		i = start[pos++];
		j = start[pos++];
		val = start[pos++];
		pNetworkData->missing(i - 1, j - 1, observation, val);
	}

	// Structural ties

	BIPARTITEVALS = VECTOR_ELT(BIPARTITE, 2);
	start = INTEGER(BIPARTITEVALS);
	listlen = ncols(BIPARTITEVALS);
	pos = 0;

	for (int row = 0; row < listlen; row++)
	{
		int i;
		int j;
		int val;
		i = start[pos++];
		j = start[pos++];
		val = start[pos++];
		pNetworkData->structural(i - 1, j - 1, observation, val);
	}
}


/**
 * Create all observations for a bipartite Network
 *
 */
void setupBipartiteObservations(SEXP BIPARTITES,
	NetworkLongitudinalData *
	pNetworkLongitudinalData)

{
    int observations = length(BIPARTITES);
    if (observations != pNetworkLongitudinalData->observationCount())
    {
		error ("wrong number of observations in bipartite");
    }
    SEXP uo;
    PROTECT(uo = install("uponly"));
    SEXP uponly = getAttrib(BIPARTITES, uo);
    SEXP dow;
    PROTECT(dow = install("downonly"));
    SEXP downonly = getAttrib(BIPARTITES, dow);
    for (int period = 0; period < (observations - 1); period++)
    {
        pNetworkLongitudinalData->upOnly(period,
			LOGICAL(uponly)[period]);
        pNetworkLongitudinalData->downOnly(period,
			LOGICAL(downonly)[period]);
    }
    for (int period = 0; period < observations; period++)
    {
    	setupBipartiteNetwork(VECTOR_ELT(BIPARTITES, period),
			pNetworkLongitudinalData,
			period);
    }
    UNPROTECT(2);
}
/**
 * Create one group of bipartite Networks
 *
 */
void setupBipartiteGroup(SEXP BIPARTITEGROUP, Data * pData)
{
    int nBipartite = length(BIPARTITEGROUP);

    for (int bipartite = 0; bipartite < nBipartite; bipartite++)
    {
        SEXP as;
        PROTECT(as = install("nodeSet"));
        SEXP actorSet = getAttrib(VECTOR_ELT(BIPARTITEGROUP, bipartite), as);
        SEXP nm;
        PROTECT(nm = install("name"));
        SEXP name = getAttrib(VECTOR_ELT(BIPARTITEGROUP, bipartite), nm);
		SEXP avout;
        PROTECT(avout = install("averageOutDegree"));
        SEXP averageOutDegree = getAttrib(VECTOR_ELT(BIPARTITEGROUP,
				bipartite), avout);
        const ActorSet * pSenders = pData->pActorSet(CHAR(STRING_ELT(
					actorSet, 0)));
        const ActorSet * pReceivers = pData->pActorSet(CHAR(STRING_ELT(
					actorSet, 1)));
		NetworkLongitudinalData *  pNetworkLongitudinalData =
			pData->createNetworkData(CHAR(STRING_ELT(name, 0)),
				pSenders, pReceivers);
        pNetworkLongitudinalData->averageOutDegree(*(REAL(averageOutDegree)));
		setupBipartiteObservations(VECTOR_ELT(BIPARTITEGROUP, bipartite),
			pNetworkLongitudinalData);

		// Once all network data has been stored, calculate some
		// statistical properties of that data.

		pNetworkLongitudinalData->calculateProperties();
        UNPROTECT(3);
    }
}

/**
 * Create all observations for a behavior Network
 *
 */
void setupBehavior(SEXP BEHAVIOR, BehaviorLongitudinalData * pBehaviorData)

{
    int observations = ncols(VECTOR_ELT(BEHAVIOR, 0));

    if (observations != pBehaviorData->observationCount())
    {
		error ("wrong number of observations in Behavior");
    }
    int nActors = nrows(VECTOR_ELT(BEHAVIOR, 0));

    if (nActors != pBehaviorData->n())
    {
        error ("wrong number of actors");
    }
    int * start = INTEGER(VECTOR_ELT(BEHAVIOR, 0));
	int * missing = LOGICAL(VECTOR_ELT(BEHAVIOR, 1));

    for (int period = 0; period < observations; period++)
    {
        for (int actor = 0; actor < nActors; actor++)
        {
            pBehaviorData->value(period, actor, *start++);
			pBehaviorData->missing(period, actor, *missing++);
       }
    }
    SEXP uo;
    PROTECT(uo = install("uponly"));
    SEXP uponly = getAttrib(VECTOR_ELT(BEHAVIOR, 0), uo);
    SEXP dow;
    PROTECT(dow = install("downonly"));
    SEXP downonly = getAttrib(VECTOR_ELT(BEHAVIOR,0), dow);
    for (int period = 0; period < (observations - 1); period++)
    {
        pBehaviorData->upOnly(period, LOGICAL(uponly)[period]);
        pBehaviorData->downOnly(period, LOGICAL(downonly)[period]);
    }
    SEXP sim;
    PROTECT(sim = install("simMean"));
    SEXP simMean = getAttrib(VECTOR_ELT(BEHAVIOR,0), sim);
	pBehaviorData->similarityMean(REAL(simMean)[0]);
	SEXP sims;
	PROTECT(sims = install("simMeans"));
	SEXP simMeans = getAttrib(VECTOR_ELT(BEHAVIOR, 0), sims);
	SEXP simNames;
	PROTECT(simNames = getAttrib(simMeans, R_NamesSymbol));
	int numberNetworks = length(simMeans);
	for (int net = 0; net < numberNetworks; net++)
	{
		pBehaviorData->similarityMeans(REAL(simMeans)[net],
			CHAR(STRING_ELT(simNames, net)));
	}

    // Now that the values are set, calculate some important statistics
	pBehaviorData->calculateProperties();
	UNPROTECT(5);
}
/**
 * Create one group of Behavior Networks
 *
 */
void setupBehaviorGroup(SEXP BEHGROUP, Data *pData)
{
    int nBehavior = length(BEHGROUP);

    for (int behavior= 0; behavior < nBehavior; behavior++)
    {
		SEXP as;
		PROTECT(as = install("nodeSet"));
        SEXP actorSet = getAttrib(VECTOR_ELT(VECTOR_ELT(BEHGROUP, behavior), 0),
								  as);

        SEXP nm;
        PROTECT(nm = install("name"));
        SEXP name = getAttrib(VECTOR_ELT(VECTOR_ELT(BEHGROUP, behavior), 0),
							  nm);

        const ActorSet * myActorSet = pData->pActorSet(CHAR(STRING_ELT(
                                                                actorSet, 0)));
	BehaviorLongitudinalData * pBehaviorData =
	    pData->createBehaviorData(CHAR(STRING_ELT(name, 0)), myActorSet);
//	Rprintf("%x\n", pBehaviorData);
	setupBehavior(VECTOR_ELT(BEHGROUP, behavior), pBehaviorData);
        UNPROTECT(2);
    }
}
/**
 * Create a constant covariate
 *
 */
void setupConstantCovariate(SEXP COCOVAR, ConstantCovariate *
			    pConstantCovariate)

{
    int nActors = length(COCOVAR);
	// Rprintf("%x\n", pConstantCovariate);
    double * start = REAL(COCOVAR);
    for (int actor = 0; actor < nActors; actor++)
    {
		double value = *start++;
		if (ISNAN(value))
		{
			pConstantCovariate->value(actor, 0);
			pConstantCovariate->missing(actor, 1);
		}
		else
		{
			pConstantCovariate->value(actor, value);
			pConstantCovariate->missing(actor, 0);

		}
    }
}
/**
 * Create one group of constant covariates
 *
 */
void setupConstantCovariateGroup(SEXP COCOVARGROUP, Data *pData)
{
    int nConstantCovariate = length(COCOVARGROUP);
//    Rprintf("nConstantCovariate %d\n", nConstantCovariate);
    for (int constantCovariate = 0; constantCovariate < nConstantCovariate;
	 constantCovariate++)
    {
        SEXP as;
        PROTECT(as = install("nodeSet"));
        SEXP actorSet = getAttrib(VECTOR_ELT(COCOVARGROUP, constantCovariate),
                                  as);
        SEXP nm;
        PROTECT(nm = install("name"));
        SEXP name = getAttrib(VECTOR_ELT(COCOVARGROUP, constantCovariate), nm);
        const ActorSet * myActorSet = pData->pActorSet(CHAR(STRING_ELT(
                                                                actorSet, 0)));
        int nActors = length(VECTOR_ELT(COCOVARGROUP, constantCovariate));
//    Rprintf("nactors %d\n", nActors);

        if (nActors != myActorSet->n())
        {
            error ("wrong number of actors");
        }
		ConstantCovariate * pConstantCovariate =
			pData->createConstantCovariate(CHAR(STRING_ELT(name, 0)),
				myActorSet);
        setupConstantCovariate(VECTOR_ELT(COCOVARGROUP,	constantCovariate),
			pConstantCovariate);
		SEXP sim;
		PROTECT(sim = install("simMean"));
		SEXP simMean = getAttrib(VECTOR_ELT(COCOVARGROUP, constantCovariate),
			sim);
		pConstantCovariate->similarityMean(REAL(simMean)[0]);
		SEXP sims;
		PROTECT(sims = install("simMeans"));
		SEXP simMeans = getAttrib(VECTOR_ELT(COCOVARGROUP, constantCovariate),
			sims);
		SEXP simNames;
		PROTECT(simNames = getAttrib(simMeans, R_NamesSymbol));
		int numberNetworks = length(simMeans);
		for (int net = 0; net < numberNetworks; net++)
		{
			pConstantCovariate->similarityMeans(REAL(simMean)[net],
				CHAR(STRING_ELT(simNames, net)));
		}
		SEXP range;
		PROTECT(range = install("range"));
		SEXP Range = getAttrib(VECTOR_ELT(COCOVARGROUP, constantCovariate),
			range);
		pConstantCovariate->range(REAL(Range)[0]);
        UNPROTECT(6);
    }
}
/**
 * Create all observations for a changing covariate
 *
 */
void setupChangingCovariate(SEXP VARCOVAR,
			    ChangingCovariate * pChangingCovariate)

{
    int observations = ncols(VARCOVAR);
    int nActors = nrows(VARCOVAR);
	double * start = REAL(VARCOVAR);
	for (int period = 0; period < observations; period++)
    {
		for (int actor = 0; actor < nActors; actor++)
		{
			double value = *start++;
			if (ISNAN(value))
			{
				pChangingCovariate->value(actor, period,
					0);
				pChangingCovariate->missing(actor, period,
					1);
			}
			else
			{
				pChangingCovariate->value(actor, period,
					value);
				pChangingCovariate->missing(actor, period,
					0);
			}
		}
    }
}
/**
 * Create one group of changing covariates
 *
 */
void setupChangingCovariateGroup(SEXP VARCOVARGROUP, Data *pData)
{
    if (length(VARCOVARGROUP) == 0)
        return;
    int observations = ncols(VECTOR_ELT(VARCOVARGROUP,0));
    if (observations != pData->observationCount() - 1)
    {
		error ("wrong number of observations in Changing Covariate");
    }
    int nChangingCovariate = length(VARCOVARGROUP);
	for (int changingCovariate = 0;
		 changingCovariate < nChangingCovariate;
		 changingCovariate++)
	{
		SEXP as;
		PROTECT(as = install("nodeSet"));
		SEXP actorSet = getAttrib(VECTOR_ELT(VARCOVARGROUP, changingCovariate),
			as);
		SEXP nm;
		PROTECT(nm = install("name"));
		SEXP name = getAttrib(VECTOR_ELT(VARCOVARGROUP, changingCovariate),
			nm);
		const ActorSet * myActorSet = pData->pActorSet(CHAR(STRING_ELT(
					actorSet, 0)));
		int nActors = nrows(VECTOR_ELT(VARCOVARGROUP,changingCovariate));

		if (nActors != myActorSet->n())
		{
			error ("wrong number of actors");
		}
		ChangingCovariate * pChangingCovariate =
			pData->createChangingCovariate(CHAR(STRING_ELT(name, 0)),
				myActorSet);
		setupChangingCovariate(VECTOR_ELT(VARCOVARGROUP, changingCovariate),
			pChangingCovariate);
		SEXP sim;
		PROTECT(sim = install("simMean"));
		SEXP simMean = getAttrib(VECTOR_ELT(VARCOVARGROUP, changingCovariate),
			sim);
		pChangingCovariate->similarityMean(REAL(simMean)[0]);
		SEXP sims;
		PROTECT(sims = install("simMeans"));
		SEXP simMeans = getAttrib(VECTOR_ELT(VARCOVARGROUP, changingCovariate),
			sims);
		SEXP simNames;
		PROTECT(simNames = getAttrib(simMeans, R_NamesSymbol));
		int numberNetworks = length(simMeans);
		for (int net = 0; net < numberNetworks; net++)
		{
			pChangingCovariate->similarityMeans(REAL(simMean)[net],
				CHAR(STRING_ELT(simNames, net)));
		}
		SEXP range;
		PROTECT(range = install("range"));
		SEXP Range = getAttrib(VECTOR_ELT(VARCOVARGROUP, changingCovariate),
			range);
		pChangingCovariate->range(REAL(Range)[0]);
		UNPROTECT(6);
	}
}
/**
 * Create a constant dyadic covariate
 *
 */
void setupDyadicCovariate(SEXP DYADVAR,
                          ConstantDyadicCovariate * pConstantDyadicCovariate)
{
    double *start = REAL(VECTOR_ELT(DYADVAR, 0));
    double *missingstart = REAL(VECTOR_ELT(DYADVAR, 1));
    int listlen = ncols(VECTOR_ELT(DYADVAR, 0));
//	Rprintf("listlen =  %d\n", listlen);
    int pos = 0;
    for (int row = 0; row < listlen; row++)
    {
		int i;
		int j;
		double val;
		i = start[pos++];
		j = start[pos++];
		val = start[pos++];
		pConstantDyadicCovariate->value(i-1, j-1, val);
	}
    listlen = ncols(VECTOR_ELT(DYADVAR, 1));
    pos = 0;
    for (int row = 0; row < listlen; row++)
    {
		int i;
		int j;
		double val;
		i = missingstart[pos++];
		j = missingstart[pos++];
		val = missingstart[pos++];
		pConstantDyadicCovariate->missing(i-1, j-1, val);
	}
}

/**
 * Create one group of constant dyadic covariates
 *
 */
void setupDyadicCovariateGroup(SEXP DYADVARGROUP, Data *pData)
{
    int nDyadicCovariate = length(DYADVARGROUP);
//    Rprintf("nDyadicCovariate %d\n", nDyadicCovariate);
    for (int dyadicCovariate = 0; dyadicCovariate < nDyadicCovariate;
	 dyadicCovariate++)
    {
        SEXP as;
        PROTECT(as = install("nodeSet"));
        SEXP actorSet = getAttrib(VECTOR_ELT(DYADVARGROUP, dyadicCovariate),
                                  as);
        SEXP nm;
        PROTECT(nm = install("name"));
        SEXP name = getAttrib(VECTOR_ELT(DYADVARGROUP, dyadicCovariate),
                              nm);
        const ActorSet * myActorSet1 = pData->pActorSet(CHAR(STRING_ELT(
                                                                actorSet, 0)));
        const ActorSet * myActorSet2 = pData->pActorSet(CHAR(STRING_ELT(
                                                                 actorSet, 1)));
		ConstantDyadicCovariate * pConstantDyadicCovariate =
			pData->createConstantDyadicCovariate(CHAR(STRING_ELT(name, 0)),
				myActorSet1, myActorSet2);
		setupDyadicCovariate(VECTOR_ELT(DYADVARGROUP, dyadicCovariate),
			pConstantDyadicCovariate);
		SEXP mean;
		PROTECT(mean = install("mean"));
		SEXP Mean = getAttrib(VECTOR_ELT(DYADVARGROUP, dyadicCovariate),
			mean);
		pConstantDyadicCovariate->mean(REAL(Mean)[0]);
        UNPROTECT(3);
    }
}

/**
 * Unpack one set of values for a changing dyadic covariate
 *
 */
void unpackChangingDyadicPeriod(SEXP VARDYADVALS, ChangingDyadicCovariate *
                                pChangingDyadicCovariate, int period)
{
    double *start = REAL(VECTOR_ELT(VARDYADVALS, 0));
	int listlen = ncols(VECTOR_ELT(VARDYADVALS, 0));
//	Rprintf("listlen =  %d\n", listlen);
    int pos = 0;
    for (int row = 0; row < listlen; row++)
    {
		int i;
		int j;
		double val;
		i = start[pos++];
		j = start[pos++];
		val = start[pos++];
		pChangingDyadicCovariate->value(i - 1, j - 1, period, val);
    }
    double *missingstart = REAL(VECTOR_ELT(VARDYADVALS, 1));
    listlen = ncols(VECTOR_ELT(VARDYADVALS, 1));
//	Rprintf("listlen =  %d\n", listlen);
	pos = 0;
    for (int row = 0; row < listlen; row++)
    {
		int i;
		int j;
		double val;
		i = missingstart[pos++];
		j = missingstart[pos++];
		val = missingstart[pos++];
		pChangingDyadicCovariate->missing(i - 1, j - 1, period, val);
    }
}
/**
 * Create all observations for a changing dyadic covariate
 *
 */
void setupChangingDyadicObservations(SEXP VARDYAD,
			      ChangingDyadicCovariate *
                              pChangingDyadicCovariate)

{
    int observations = length(VARDYAD);
    //   if (observations != pworkLongitudinalData->observationCount())
    // {
    //	error ("wrong number of observations in OneMode");
    //  }
    for (int period = 0; period < (observations - 1); period++)
    {
		unpackChangingDyadicPeriod(VECTOR_ELT(VARDYAD, period),
			pChangingDyadicCovariate, period);
    }
}
/**
 * Create one group of changing dyadic covariates
 *
 */
void setupChangingDyadicCovariateGroup(SEXP VARDYADGROUP, Data * pData)
{
    int nChangingDyadic = length(VARDYADGROUP);

    for (int changingDyadic = 0; changingDyadic < nChangingDyadic;
         changingDyadic++)
    {
        SEXP as;
        PROTECT(as = install("nodeSet"));
        SEXP actorSet = getAttrib(VECTOR_ELT(VARDYADGROUP, changingDyadic), as);
        SEXP nm;
        PROTECT(nm = install("name"));
        SEXP name = getAttrib(VECTOR_ELT(VARDYADGROUP, changingDyadic), nm);
        const ActorSet * myActorSet1 = pData->pActorSet(CHAR(STRING_ELT(
                                                                actorSet, 0)));
        const ActorSet * myActorSet2 = pData->pActorSet(CHAR(STRING_ELT(
                                                                actorSet, 1)));
		ChangingDyadicCovariate *  pChangingDyadicCovariate =
			pData->createChangingDyadicCovariate(CHAR(STRING_ELT(name, 0)),
				myActorSet1, myActorSet2);
		setupChangingDyadicObservations(VECTOR_ELT(VARDYADGROUP,
				changingDyadic),
			pChangingDyadicCovariate);
 		SEXP mean;
		PROTECT(mean = install("mean"));
		SEXP Mean = getAttrib(VECTOR_ELT(VARDYADGROUP, changingDyadic),
			mean);
		pChangingDyadicCovariate->mean(REAL(Mean)[0]);
		UNPROTECT(3);
    }
}
/**
 * Create the exogenous composition change events for one actor set within
 * one group.
 */
void setupExogenousEventSet(SEXP EXOGEVENTSET, Data *pData)
{
    /* pass in the data for one actor set as two items: first
	   a list of columns: event type, period, actor, time.
	   Secondly a matrix of booleans, indicating whether active at start
	   of period.*/

	/* first find the actor set */
	SEXP as;
	PROTECT(as = install("nodeSet"));
	SEXP actorSet = getAttrib(EXOGEVENTSET, as);

	/* now process the events */
	SEXP EVENTS = VECTOR_ELT(EXOGEVENTSET, 0);
    int nEvents = length(VECTOR_ELT(EVENTS, 0));
    //Rprintf("number of rows of data frame %d\n",nEvents);
	//Rprintf("%d\n", length(EVENTS));
    int * type = INTEGER(VECTOR_ELT(EVENTS, 0));
	//Rprintf("type %d\n",*type);
	int * period = INTEGER(VECTOR_ELT(EVENTS, 1));
	//Rprintf("period %d\n",*period);
	int * actor = INTEGER(VECTOR_ELT(EVENTS, 2));
	//Rprintf("actor %d\n",*actor);
    double * time = REAL(VECTOR_ELT(EVENTS, 3));
	//Rprintf("time %5.4f\n",*time);
	const ActorSet * myActorSet = pData->pActorSet(CHAR(STRING_ELT(actorSet,
															   0)));
    for (int event = 0; event < nEvents; event++)
    {
		if (*type == 1)
		{
				pData->addJoiningEvent(*period-1, myActorSet, *actor-1, *time);
		}
		else
		{
			pData->addLeavingEvent(*period-1, myActorSet, *actor-1, *time);
		}
		type++;
 		period++;
 		actor++;
 		time++;
  }
/* retrieve some to check*/
//     const EventSet * myeventset= pData->pEventSet(0);
// 	EventSet::iterator myit = myeventset->begin();
// 	Rprintf("period 1 first event? %d %3.2f\n",(*myit)->actor(),
// 			(*myit)->time());

	/* initialise the active flags */
	SEXP ACTIVES= VECTOR_ELT(EXOGEVENTSET, 1);

	/* this is a matrix with column for each observation and row per actor*/

	int nActors = myActorSet->n();
	int *active = LOGICAL(ACTIVES);
	for (int period = 0; period < pData->observationCount(); period++)
	{
			for (int actor = 0; actor < nActors; actor++)
			{
				pData->active(myActorSet, actor, period, *active);
				active++;
			}

	}
	UNPROTECT(1);

}
/**
 * Create one group of exogenous composition change events
 *
 */
void setupExogenousEventGroup(SEXP EXOGEVENTGROUP, Data *pData)
{
    /* pass in the data for each actor set as two items: first
	   a list of columns: event type, period, actor, time.
	   Secondly a matrix of booleans, indicating whether active at start
	   of period.*/

    int nActorSets = length(EXOGEVENTGROUP);
//	Rprintf("number of actor sets %d\n", nActorSets);

	for (int actorSet = 0; actorSet < nActorSets; actorSet++)
	{
		setupExogenousEventSet(VECTOR_ELT(EXOGEVENTGROUP, actorSet), pData);
	}

}
/**
 *  Creates all the basic effects for one network
 */
SEXP createEffects(SEXP EFFECTS, Model *pModel, vector<Data *> * pGroupData,
					   const char *networkName, int effectCol,
					   int parmCol, int int1Col, int int2Col,
					   int initValCol, int typeCol, int groupCol,
					   int periodCol, int rateTypeCol,
					   int netTypeCol)
    {
        // find out how many effects there are
        int nEffects = length(VECTOR_ELT(EFFECTS, 0));

        // create the effects

		/* set up a vector to return the pointers in */
		SEXP effectPtrs;
		PROTECT(effectPtrs = allocVector(VECSXP, nEffects));

		for (int i = 0; i < nEffects; i++)
		{
			EffectInfo * pEffectInfo = 0;

			const char * effectName =
				CHAR(STRING_ELT(VECTOR_ELT(EFFECTS, effectCol), i));
			int parm1 = INTEGER(VECTOR_ELT(EFFECTS, parmCol))[i];
			double parm = parm1;
			const char * interaction1 =
				CHAR(STRING_ELT(VECTOR_ELT(EFFECTS, int1Col), i));
			const char * interaction2 =
				CHAR(STRING_ELT(VECTOR_ELT(EFFECTS, int2Col), i));
            double initialValue = REAL(VECTOR_ELT(EFFECTS, initValCol))[i];
			const char * effectType =
				CHAR(STRING_ELT(VECTOR_ELT(EFFECTS, typeCol), i));
			const char * rateType =
				CHAR(STRING_ELT(VECTOR_ELT(EFFECTS, rateTypeCol), i));
			const char * netType =
				CHAR(STRING_ELT(VECTOR_ELT(EFFECTS, netTypeCol), i));

			if (strcmp(effectType, "rate") == 0 &&
				strcmp(effectName, "Rate") == 0)
			{
				/* find the network */
				int group = INTEGER(VECTOR_ELT(EFFECTS, groupCol))[i] - 1;
				int period = INTEGER(VECTOR_ELT(EFFECTS, periodCol))[i] - 1;

				Data * pData = (*pGroupData)[group];

				if (strcmp(netType, "behavior") != 0)
				{
					NetworkLongitudinalData * pNetwork =
						pData->pNetworkData(networkName);
					pModel->basicRateParameter(pNetwork, period, initialValue);
				}
				else //if (strcmp(netType, "Behavior") == 0)
				{
					BehaviorLongitudinalData * pNetwork =
						pData->pBehaviorData(networkName);
					pModel->basicRateParameter(pNetwork, period, initialValue);
				}
			}
			else
			{
				pEffectInfo = pModel->addEffect(networkName,
					effectName,
					effectType,
					initialValue,
					parm,
					interaction1,
					interaction2,
					rateType);
			}

			SET_VECTOR_ELT(effectPtrs, i,
				R_MakeExternalPtr((void *) pEffectInfo,
					R_NilValue, R_NilValue));
		}

		UNPROTECT(1);
		return effectPtrs;
	}
/**
 *  Creates all the interaction effects for one network
 */
SEXP createInteractionEffects(SEXP EFFECTS, Model *pModel,
	const char *networkName, int effectCol, int initValCol,
	int typeCol, int intptr1Col, int intptr2Col, int intptr3Col)
{
	// find out how many effects there are
	int nEffects = length(VECTOR_ELT(EFFECTS, 0));

	// create the effects

	/* set up a vector to return the pointers in */
	SEXP effectPtrs;
	PROTECT(effectPtrs = allocVector(VECSXP, nEffects));

	for (int i = 0; i < nEffects; i++)
	{
		EffectInfo * pEffectInfo = 0;

		const char * effectName =
			CHAR(STRING_ELT(VECTOR_ELT(EFFECTS, effectCol), i));
		double initialValue = REAL(VECTOR_ELT(EFFECTS, initValCol))[i];
		const char * effectType =
			CHAR(STRING_ELT(VECTOR_ELT(EFFECTS, typeCol), i));
		EffectInfo * pEffect1 = (EffectInfo *) R_ExternalPtrAddr(
			VECTOR_ELT(VECTOR_ELT(EFFECTS, intptr1Col), i));
		EffectInfo * pEffect2 = (EffectInfo *) R_ExternalPtrAddr(
			VECTOR_ELT(VECTOR_ELT(EFFECTS, intptr2Col), i));
		EffectInfo * pEffect3 = 0;
		if (!isNull(VECTOR_ELT(VECTOR_ELT(EFFECTS, intptr3Col), i)))
		{
			pEffect3 = (EffectInfo *) R_ExternalPtrAddr(
				VECTOR_ELT(VECTOR_ELT(EFFECTS, intptr3Col), i));
		}

		pEffectInfo = pModel->addInteractionEffect(networkName,
			effectName,
			effectType,
			initialValue,
			pEffect1,
			pEffect2,
			pEffect3);

		SET_VECTOR_ELT(effectPtrs, i,
			R_MakeExternalPtr((void *) pEffectInfo,
				R_NilValue, R_NilValue));
	}

	UNPROTECT(1);
	return effectPtrs;
}

/**
 *  Retrieves the values of the statistics and scores for each of the effects,
 *  for one period. The call will relate to one group only, although all effects
 *  are the same apart from the basic rates. Not used in maximum likelihood.
 */
void getStatistics(SEXP EFFECTSLIST,
	const StatisticCalculator * pCalculator,
	int period, int group, const Data *pData,
	const EpochSimulation * pEpochSimulation,
	vector<double> * rfra, vector<double> *rscore)
{

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


	double statistic = 0;
	double score = 0;
	int istore = 0;

	for (int ii = 0; ii < length(EFFECTSLIST); ii++)
	{
		const char * networkName =
			CHAR(STRING_ELT(VECTOR_ELT(VECTOR_ELT(EFFECTSLIST, ii),
						nameCol), 0));
		SEXP EFFECTS = VECTOR_ELT(EFFECTSLIST, ii);

		for (int i = 0; i < length(VECTOR_ELT(EFFECTS,0)); i++)
		{
			const char * effectName =
				CHAR(STRING_ELT(VECTOR_ELT(EFFECTS, effectCol),  i));
			//	int parm = INTEGER(VECTOR_ELT(EFFECTS, parmCol))[i];
			const char * interaction1 =
				CHAR(STRING_ELT(VECTOR_ELT(EFFECTS, int1Col),i));
			//	const char * interaction2 =
			//CHAR(STRING_ELT(VECTOR_ELT(EFFECTS, int2Col), i));
			//double initialValue =
			//REAL(VECTOR_ELT(EFFECTS, initValCol))[i];
			const char * effectType =
				CHAR(STRING_ELT(VECTOR_ELT(EFFECTS, typeCol), i));
			const char * netType =
				CHAR(STRING_ELT(VECTOR_ELT(EFFECTS, netTypeCol), i));
			const char * rateType =
				CHAR(STRING_ELT(VECTOR_ELT(EFFECTS, rateTypeCol), i));
			//	Rprintf("%s %s \n", effectType, netType);
			if (strcmp(effectType, "rate") == 0)
			{
				if (strcmp(effectName, "Rate") == 0)
				{
					int groupno =
						INTEGER(VECTOR_ELT(EFFECTS, groupCol))[i];
					int periodno =
						INTEGER(VECTOR_ELT(EFFECTS, periodCol))[i];
					if ((periodno - 1) == period && (groupno - 1) == group)
					{

						if (strcmp(netType, "behavior") == 0)
						{
							LongitudinalData * pNetworkData =
								pData->pBehaviorData(networkName);
							statistic = pCalculator->distance(pNetworkData,
								period);
							//	Rprintf("%f behavior dist\n", statistic);
							if (pEpochSimulation)
							{
								const DependentVariable * pVariable =
									pEpochSimulation->pVariable(networkName);
								score = pVariable->basicRateScore();
							}
							else
							{
								score = 0;
							}
						}
						else
						{
							LongitudinalData * pNetworkData =
								pData->pNetworkData(networkName);
							statistic = pCalculator->distance(pNetworkData,
								period);
							if (pEpochSimulation)
							{
								/* find the dependent variable for the score */
								const DependentVariable * pVariable =
									pEpochSimulation->pVariable(networkName);
								score = pVariable->basicRateScore();
							}
							else
							{
								score = 0;
							}
						}
					}
					else
					{
						statistic = 0;
						score = 0;
					}
				}
				else if (strcmp(rateType, "structural") == 0)
				{
					EffectInfo * pEffectInfo = (EffectInfo *)
						R_ExternalPtrAddr(
							VECTOR_ELT(VECTOR_ELT(EFFECTS,
									pointerCol), i));
					statistic = pCalculator->statistic(pEffectInfo);
					if (pEpochSimulation)
					{
						const DependentVariable * pVariable =
							pEpochSimulation->pVariable(networkName);
						const NetworkVariable * pNetworkVariable;
						if (strcmp(interaction1, "") == 0)
						{
							pNetworkVariable =
								(const NetworkVariable *)
								pEpochSimulation->pVariable(networkName);
						}
						else
						{
							pNetworkVariable =
								(const NetworkVariable *)
								pEpochSimulation->pVariable(interaction1);
						}
						if (strcmp(effectName, "outRate") == 0)
						{
							score =
								pVariable->outDegreeScore(pNetworkVariable);
						}
						else if (strcmp(effectName, "inRate") == 0)
						{
							score =
								pVariable->inDegreeScore(pNetworkVariable);
						}
						else if (strcmp(effectName, "recipRate") == 0)
						{
							score =
								pVariable->reciprocalDegreeScore(pNetworkVariable);
						}
						else if (strcmp(effectName, "outRateInv") == 0)
						{
							score =
								pVariable->inverseOutDegreeScore(pNetworkVariable);
						}
						else
						{

							error("Unexpected rate effect %s\n",
								effectName);
						}
					}
					else
					{
						score = 0;
					}
				}
				else
				{
					EffectInfo * pEffectInfo = (EffectInfo *)
						R_ExternalPtrAddr(
							VECTOR_ELT(VECTOR_ELT(EFFECTS,
									pointerCol), i));
					statistic = pCalculator->statistic(pEffectInfo);

					if (pEpochSimulation)
					{
						ConstantCovariate * pConstantCovariate =
							pData->pConstantCovariate(interaction1);
						ChangingCovariate * pChangingCovariate =
							pData->pChangingCovariate(interaction1);
						BehaviorVariable * pBehavior =
							(BehaviorVariable *)
							pEpochSimulation->pVariable(interaction1);
						//find the network

						const DependentVariable * pVariable =
							pEpochSimulation->pVariable(networkName);

						if (pConstantCovariate)
						{
							score = pVariable->constantCovariateScore(
								pConstantCovariate);
						}
						else if (pChangingCovariate)
						{
							score = pVariable->changingCovariateScore(
								pChangingCovariate);
						}
						else if (pBehavior)
						{
							score = pVariable->behaviorVariableScore(
								pBehavior);
						}
						else
						{
							error("No individual covariate named %s.",
								interaction1);
						}
					}
					else
					{
						score = 0;
					}
				}

			}
			else
			{
				if (strcmp(effectType, "eval") == 0)
				{
					EffectInfo * pEffectInfo = (EffectInfo *)
						R_ExternalPtrAddr(
							VECTOR_ELT(VECTOR_ELT(EFFECTS,
									pointerCol), i));
					statistic =	pCalculator->statistic(pEffectInfo);
					if (pEpochSimulation)
					{
						score = pEpochSimulation->score(pEffectInfo);
					}
					else
					{
						score = 0;
					}
				}
				else if (strcmp(effectType, "endow") == 0)
				{
					EffectInfo * pEffectInfo = (EffectInfo *)
						R_ExternalPtrAddr(
							VECTOR_ELT(VECTOR_ELT(EFFECTS,
									pointerCol), i));
					statistic = pCalculator->statistic(pEffectInfo);
					if (strcmp(netType, "behavior") != 0)
					{
						statistic = -1 * statistic;
					}
					if (pEpochSimulation)
					{
						score = pEpochSimulation->score(pEffectInfo);
					}
					else
					{
						score = 0;
					}
				}
				else if (strcmp(effectType, "creation") == 0)
				{
					EffectInfo * pEffectInfo = (EffectInfo *)
						R_ExternalPtrAddr(
							VECTOR_ELT(VECTOR_ELT(EFFECTS,
									pointerCol), i));
					statistic = pCalculator->statistic(pEffectInfo);
					if (pEpochSimulation)
					{
						score = pEpochSimulation->score(pEffectInfo);
					}
					else
					{
						score = 0;
					}
				}
				else
				{
					error("invalid effect type %s\n", effectType);
				}
			}
			(*rfra)[istore] = statistic;
			if (pEpochSimulation)
			{
				//		Rprintf("%f %d \n", score, istore);
				(*rscore)[istore] = score;
			}
			istore++; /* keep forgetting to move the ++ */
		}
	}
	UNPROTECT(1);
}

/**
 *  retrieves the values of the scores and derivatives for each of the effects,
 *  for one period. The call will relate to one group only, although all effects
 *  are the same apart from the basic rates. Only used in maximum likelihood.
 */
void getScores(SEXP EFFECTSLIST, int period, int group,
	const MLSimulation * pMLSimulation,
	vector<double> * rderiv, vector<double> *rscore)
{

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

	int storescore = 0;
	int storederiv = 0;

	for (int ii = 0; ii < length(EFFECTSLIST); ii++)
	{
		const char * networkName =
			CHAR(STRING_ELT(VECTOR_ELT(VECTOR_ELT(EFFECTSLIST, ii),
						nameCol), 0));
		SEXP EFFECTS = VECTOR_ELT(EFFECTSLIST, ii);

		for (int i = 0; i < length(VECTOR_ELT(EFFECTS,0)); i++)
		{
			const char * effectName =
				CHAR(STRING_ELT(VECTOR_ELT(EFFECTS, effectCol),  i));
			const char * effectType =
				CHAR(STRING_ELT(VECTOR_ELT(EFFECTS, typeCol), i));
			if (strcmp(effectType, "rate") == 0)
			{
				if (strcmp(effectName, "Rate") == 0)
				{
					int groupno =
						INTEGER(VECTOR_ELT(EFFECTS, groupCol))[i];
					int periodno =
						INTEGER(VECTOR_ELT(EFFECTS, periodCol))[i];
					if ((periodno - 1) == period && (groupno - 1) == group)
					{
						const DependentVariable * pVariable =
							pMLSimulation->pVariable(networkName);
						(*rscore)[storescore++] = pVariable->basicRateScore();
						(*rderiv)[storederiv++] =
							pVariable->basicRateDerivative();
					}
					else
					{
						(*rscore)[storescore++] = 0;
						(*rderiv)[storederiv++] = 0;
					}
				}
				else
				{
					error("Non constant rate effects are not yet %s",
						"implemented for maximum likelihood.");
				}
			}
			else
			{
				EffectInfo * pEffectInfo = (EffectInfo *)
					R_ExternalPtrAddr(VECTOR_ELT(VECTOR_ELT(EFFECTS,
								pointerCol), i));
				(*rscore)[storescore++] = pMLSimulation->score(pEffectInfo);

				// get the map of derivatives
				map<const EffectInfo *, double > deriv =
					pMLSimulation->derivative(pEffectInfo);

				for (int j = 0; j < length(VECTOR_ELT(EFFECTS,0)); j++)
				{
					const char * effectType =
						CHAR(STRING_ELT(VECTOR_ELT(EFFECTS, typeCol), j));

					if (!strcmp(effectType, "rate") == 0)
					{
						//	Rprintf("%s %s \n", effectType, netType);
						EffectInfo * pEffectInfo2 = (EffectInfo *)
							R_ExternalPtrAddr(
								VECTOR_ELT(VECTOR_ELT(EFFECTS,
										pointerCol), j));

						(*rderiv)[storederiv++] =
							pMLSimulation->derivative(pEffectInfo,
								pEffectInfo2);
					}
				}
			}
		}
	}
	UNPROTECT(1);
}


