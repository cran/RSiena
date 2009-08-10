/******************************************************************************
 * SIENA: Simulation Investigation for Empirical Network Analysis
 *
 * Web: http://www.stats.ox.ac.uk/~snijders/siena/
 *
 * File: siena07.cpp
 *
 * Description: This module contains routines to interface with R,
 * setting up the Data object with data from R.
 *****************************************************************************/
/**
 * @file
 * Sets up the Data object with data from R
 */
#include <iostream>
#include <fstream>
#include <valarray>
#include <vector>
#include <set>
#include <string>
#include <R.h>
#include <Rinternals.h>
#include <Rmath.h>
#include "data/Data.h"
#include "data/OneModeNetwork.h"
#include "data/TieIterator.h"
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
#include "model/Function.h"
#include "model/State.h"
#include "model/StatisticCalculator.h"
#include "utils/Random.h"
#include "data/ActorSet.h"
#include "model/EpochSimulation.h"
#include "model/variables/DependentVariable.h"
#include "model/variables/NetworkVariable.h"
using namespace std;
using namespace siena;

/**
 * print out the data
 *
 */
void printOutData(Data *pData)
{

	ofstream myfile ("data.txt");

	if (myfile.is_open())
	{
		myfile << pData->observationCount() << endl;
		const std::vector<LongitudinalData * > rVariables =
			pData->rDependentVariableData();
		int nActors = rVariables[0]->n();

		myfile << rVariables[0]->n();
		for (unsigned i = 0; i < pData->rDependentVariableData().size(); i++)
		{
			NetworkLongitudinalData * pNetworkData =
				dynamic_cast<NetworkLongitudinalData *>(rVariables[i]);
			BehaviorLongitudinalData * pBehaviorData =
				dynamic_cast<BehaviorLongitudinalData *>(rVariables[i]);

			myfile << rVariables[i]->name() << endl;

			for(int period = 0; period < pData->observationCount(); period++)
			{
				if (pNetworkData)
				{
					if (period == 0) myfile << "oneMode" << endl;
					myfile << pNetworkData->pNetwork(period)->tieCount() << endl;
					for (TieIterator
							 iter=pNetworkData-> pNetwork(period)->ties();
						 iter.valid();
						 iter.next())
					{
						myfile << iter.ego() << " "
							   << iter.alter() << " "
						 << iter.value() << endl;
					}
					myfile << pNetworkData->
						pMissingTieNetwork(period)->tieCount() << endl;
					for (TieIterator
							 iter=pNetworkData->
							 pMissingTieNetwork(period)->ties();
						 iter.valid();
						 iter.next())
					{
						myfile << iter.ego() << " "
							   << iter.alter() << " "
						 << iter.value() << endl;
					}
					myfile << pNetworkData->
						pStructuralTieNetwork(period)->tieCount() << endl;
					for (TieIterator
							 iter=pNetworkData->
							 pStructuralTieNetwork(period)->ties();
						 iter.valid();
						 iter.next())
					{
						myfile << iter.ego() << " "
							   << iter.alter() << " "
						 << iter.value() << endl;
					}
					// other attributes uponly downonly
					myfile << pNetworkData->upOnly(period) <<  " " <<
						pNetworkData->downOnly(period) << endl;
				}

				else if (pBehaviorData)
				{
					if (period ==0 )
					{
						myfile << "behavior" << endl;
						myfile <<  pBehaviorData->n() << endl;
					}
					for (int ii = 0; ii < pBehaviorData->n(); ii++)
					{
						myfile << ii << " " <<
							pBehaviorData->value(period, ii) << endl;
					}
					for (int ii = 0; ii < pBehaviorData->n(); ii++)
					{
						myfile << ii << " " <<
							pBehaviorData->missing(period, ii) << endl;
					}
					// other attributes similarityMean uponly downonly
					myfile << pBehaviorData->similarityMean() << endl;
					myfile << pBehaviorData->upOnly(period) << " " <<
						pBehaviorData->downOnly(period) << endl;
				}
				else
				{
					throw ("Unexpected class of dependent variable");
				}

			}
		}
		//Rprintf("gothere\n");
		const std::vector<ConstantCovariate * > rConstantCovariates =
			pData->rConstantCovariates();
		for (unsigned i = 0; i < pData->rConstantCovariates().size(); i++)
		{
			ConstantCovariate * pCovariate = rConstantCovariates[i];

			myfile << pCovariate->name() << endl;
			myfile << "constantcovariate" << endl;
			myfile <<  nActors << endl;

			for (int ii = 0; ii < nActors; ii++)
			{
				Rprintf("gothere  %d %f\n", ii, pCovariate->value(ii));
				myfile << ii << " " <<	pCovariate->value(ii) << endl;
			}
			for (int ii = 0; ii < nActors; ii++)
			{
				myfile << ii << " " <<	pCovariate->missing(ii) << endl;
			}

			myfile << pCovariate->similarityMean() << endl;
			myfile << pCovariate->range() << endl;
		}
		const std::vector<ChangingCovariate * > rChangingCovariates =
			pData->rChangingCovariates();

		for (unsigned i = 0; i < pData->rChangingCovariates().size(); i++)
		{
			ChangingCovariate * pCovariate = rChangingCovariates[i];

			myfile << pCovariate->name() << endl;
			myfile << "changingcovariate" << endl;
			myfile <<  nActors << endl;

			for (int period = 0; period < pData->observationCount() - 1;
				 period++)
			{
				for (int ii = 0; ii < nActors; ii++)
				{
					myfile << ii << " " <<
						pCovariate->value(ii, period) << endl;
				}
				for (int ii = 0; ii < nActors; ii++)
				{
					myfile << ii << " " <<
						pCovariate->missing(ii, period) << endl;
				}
			}
			myfile << pCovariate->similarityMean() << endl;
			myfile << pCovariate->range() << endl;
		}

		const std::vector<ConstantDyadicCovariate * >
			rConstantDyadicCovariates =	pData->rConstantDyadicCovariates();

		for (unsigned i = 0; i < pData->rConstantDyadicCovariates().size(); i++)
		{
			ConstantDyadicCovariate * pCovariate = rConstantDyadicCovariates[i];

			myfile << pCovariate->name() << endl;
			myfile << "constantdyadiccovariate" << endl;
			myfile <<  nActors << endl;

			for (int j = 0; j < nActors; j++)
			{
				for (int k = 0; k < nActors; k++)
				{
					myfile << j << " " << k << " " << pCovariate->value(j, k)
						   << endl;
					myfile << j << " " << k << " " <<
						pCovariate->missing(j, k) << endl;
				}
			}
			myfile << pCovariate->mean() << endl;
		}
		const std::vector<ChangingDyadicCovariate * >
			rChangingDyadicCovariates =	pData->rChangingDyadicCovariates();

		for (unsigned i = 0; i < pData->rChangingDyadicCovariates().size(); i++)
		{
			ChangingDyadicCovariate * pCovariate = rChangingDyadicCovariates[i];

			myfile << pCovariate->name() << endl;
			myfile << "changingdyadiccovariate" << endl;
			myfile <<  nActors << endl;

			for (int period = 0; period < pData->observationCount() - 1;
				 period++)
			{
				for (int j = 0; j < nActors; j++)
				{
					for (int k = 0; k < nActors; k++)
					{
						myfile << j << " " << k << " " <<
							pCovariate->value(j, k, period)
							   << endl;
						myfile << j << " " << k << " " <<
							pCovariate->missing(j, k, period) << endl;
					}
				}
			}
			myfile << pCovariate->mean() << endl;
		}
	}
}
/**
 * Traps errors so R can stop the function rather than being stoppped itself.
 *
 */
void Rterminate()
{
    try
    {
		throw;
    }
    catch(exception& e)
    {
		error(e.what());
    }
}

int mysample(const int n, const valarray<double>& p)
{
    double r= runif(0,1);
    int lo= -1, hi= n-1, mid;
    while((hi-lo)>1)
    {
	mid=hi- (hi-lo)/2;
	if (mid<0)
	{
	    return(0);
	}
	else
	{
	    if (p[mid+1] < r*p[n])
	    {
		lo= mid;
	    }
	    else
	    {
		hi=mid;
	    }
	}
    }
    return(hi);
}
int sample(int n)
{
    double r= runif(0, 1);
    int sel=(int) n*r ;
    return(sel);
}
SEXP getAdjacency(const Network& net)
{
    SEXP ans;
    int n=net.n();
    PROTECT( ans=allocMatrix(INTSXP,n,n));
    int *ians = INTEGER(ans);
    /* initialise the memory: possibly only neccesary in case of error! */
    for (int i = 0; i<n*n;i++)
	ians[i]=0;
    for (TieIterator iter=net.ties(); iter.valid(); iter.next())
    {
	ians[iter.ego()+(iter.alter()*n)] = iter.value();
    }

    UNPROTECT(1);
    return(ans) ;
}
SEXP getEdgeList(const Network& net)
{
    SEXP ans;
	int nties = net.tieCount();
    PROTECT(ans = allocMatrix(INTSXP, nties, 2));
    int *ians = INTEGER(ans);
    /* initialise the memory: possibly only neccesary in case of error! */
    for (int i = 0; i < nties * 2; i++)
		ians[i] = 0;
	int irow = 0;
    for (TieIterator iter=net.ties(); iter.valid(); iter.next())
    {
		ians[irow ] = iter.ego() + 1;
		ians[ nties + irow] = iter.alter() + 1;
		irow ++;
    }

    UNPROTECT(1);
    return(ans) ;
}


/**
 * Calculates change between two networks. Ignores any values missing in
 *  either network. Does not use structural values as yet.
 */
int calcChange(const Network* net0, const Network* net1, const Network* miss0,
	       const Network* miss1)
{
    int change = 0;
    TieIterator iter0 = net0->ties();
    TieIterator iter1 = net1->ties();
    while(iter0.valid() || iter1.valid())
	{
	    if (iter0.valid() && (!iter1.valid() ||
				  (iter0.ego() < iter1.ego()) ||
				  (iter0.ego() == iter1.ego() &&
				   iter0.alter() < iter1.alter())))
	    {
		if (!miss0->tieValue(iter0.ego(), iter0.alter()) &&
		    !miss1->tieValue(iter0.ego(), iter0.alter()))
		    change++;
		iter0.next();
	    }
	    else
	    {
		if (iter1.valid() && (!iter0.valid() ||
				      (iter1.ego() < iter0.ego()) ||
				      (iter0.ego() == iter1.ego() &&
				       iter1.alter() < iter0.alter())))
		{
		    if (!miss0->tieValue(iter1.ego(),
					 iter1.alter()) &&
			!miss1->tieValue(iter1.ego(),
					 iter1.alter()))
			change++;
		    iter1.next();
		}
		else
		{
		    if (iter0.value() != iter1.value())
		    {
			if (!miss0->tieValue(iter0.ego(),
					     iter0.alter()) &&
			    !miss1->tieValue(iter0.ego(),
					     iter0.alter()))
			    change++;
		    }
		    iter0.next();
		    iter1.next();
		}
	    }
	}
    return change;
}
/**
 * Matches column names with indices.
 */
void getColNos(SEXP Names, int * netTypeCol, int * nameCol, int * effectCol,
			   int * parmCol, int * int1Col, int * int2Col, int * initValCol,
			   int * typeCol, int * groupCol, int * periodCol, int * pointerCol,
			   int * rateTypeCol)
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
//	Rprintf("%d parmcol\n", *parmCol);
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

	getColNos(Names, &netTypeCol, &nameCol, &effectCol,
			  &parmCol, &int1Col, &int2Col, &initValCol,
			  &typeCol, &groupCol, &periodCol, &pointerCol,
			  &rateTypeCol);
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
			//	int parm = INTEGER(VECTOR_ELT(EFFECTS, parmCol))[eff];
			//const char * interaction1 =
			//	CHAR(STRING_ELT(VECTOR_ELT(EFFECTS, int1Col), eff));
			//const char * interaction2 =
			//	CHAR(STRING_ELT(VECTOR_ELT(EFFECTS, int2Col), eff));
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
 * Unpack one set of ties for a onemode network
 */
void unpackOneModeNetwork(SEXP ONEMODEVALS, OneModeNetwork * pNetwork)
{
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
	pNetwork->setTieValue(i-1, j-1, val);
    }
//	Rprintf("end 1,14 %d\n",pNetwork->tieValue(0, 13));
}
/**
 * Create one observation for a one mode Network: ties, missing, structural
 *
 */
void setupOneModeNetwork(SEXP ONEMODE, OneModeNetwork * pNetwork,
			 OneModeNetwork * pMissingTieNetwork,
			 OneModeNetwork * pStructuralTieNetwork)

{
/* one mode networks are passed in as list of edgelists with attributes
   giving the size of the network - not checked yet*/
    unpackOneModeNetwork(VECTOR_ELT(ONEMODE, 0), pNetwork);
    unpackOneModeNetwork(VECTOR_ELT(ONEMODE, 1), pMissingTieNetwork);
    unpackOneModeNetwork(VECTOR_ELT(ONEMODE, 2), pStructuralTieNetwork);
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
	OneModeNetwork * pNetwork = (OneModeNetwork *)
            pOneModeNetworkLongitudinalData->pNetwork(period);
	OneModeNetwork * pMissingTieNetwork = (OneModeNetwork *)
	    pOneModeNetworkLongitudinalData->pMissingTieNetwork(period);
	OneModeNetwork * pStructuralTieNetwork = (OneModeNetwork *)
	    pOneModeNetworkLongitudinalData->pStructuralTieNetwork(period);
	setupOneModeNetwork(VECTOR_ELT(ONEMODES, period), pNetwork,
			    pMissingTieNetwork, pStructuralTieNetwork);
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
		setupOneModeObservations(VECTOR_ELT(ONEMODEGROUP, oneMode),
			pOneModeNetworkLongitudinalData);
        UNPROTECT(4);
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

    // Now that the values are set, calculate some important statistics
	pBehaviorData->calculateStatistics();
	UNPROTECT(3);
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
  Rprintf("%x\n", pConstantCovariate);
    double * start = REAL(COCOVAR);
    for (int actor = 0; actor < nActors; actor++)
    {
		double value = *start++;
		if (ISNA(value))
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
		SEXP range;
		PROTECT(range = install("range"));
		SEXP Range = getAttrib(VECTOR_ELT(COCOVARGROUP, constantCovariate),
			range);
		pConstantCovariate->range(REAL(Range)[0]);
        UNPROTECT(4);
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
			if (ISNA(value))
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
		SEXP range;
		PROTECT(range = install("range"));
		SEXP Range = getAttrib(VECTOR_ELT(VARCOVARGROUP, changingCovariate),
			range);
		pChangingCovariate->range(REAL(Range)[0]);
		UNPROTECT(4);
	}
}
/**
 * Create a constant dyadic covariate
 *
 */
void setupDyadicCovariate(SEXP DYADVAR,
                          ConstantDyadicCovariate * pConstantDyadicCovariate)
{
    double *start = REAL(DYADVAR);
    int listlen = ncols(DYADVAR);
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
 	pConstantDyadicCovariate->missing(i-1, j-1, 0);
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
    double *start = REAL(VARDYADVALS);
    int listlen = ncols(VARDYADVALS);
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
extern "C"
{

/**
 *  Creates an array of pointers to Data objects, one for each group
 *  and returns the address of the array to R. Also creates the actor sets
 *  for each group
 */
    SEXP setupData(SEXP OBSERVATIONSLIST, SEXP ACTORSLIST)
    {
/* make error messages go back to R nicely */
		set_terminate(Rterminate);

		int nGroups = length(OBSERVATIONSLIST);
//	Rprintf("%d\n", nGroups);

		vector<Data *> *pGroupData = new vector <Data *>;

		for (int group = 0; group < nGroups; group++)
		{
			int observations = INTEGER(VECTOR_ELT(OBSERVATIONSLIST, group))[0];
            //   Rprintf("%d\n", observations);

            pGroupData->push_back(new Data(observations));
			int nNodeSets = length(VECTOR_ELT(ACTORSLIST, group));

			for (int nodeSet = 0; nodeSet < nNodeSets; nodeSet++)
			{
                SEXP nsn;
                PROTECT(nsn = install("nodeSetName"));
                SEXP nodeSetName = getAttrib(VECTOR_ELT(VECTOR_ELT(ACTORSLIST,
                                                                   group),
                                                        nodeSet), nsn);
                //   Rprintf("%s\n", CHAR(STRING_ELT(nodeSetName, 0)));
				// const ActorSet *pActors =
                    (*pGroupData)[group]->
                    createActorSet(CHAR(STRING_ELT(nodeSetName, 0)),
                                   length(VECTOR_ELT(VECTOR_ELT(ACTORSLIST,
                                                                group),
                                                     nodeSet)));
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
 *  Creates all the effects for one network
 */
	SEXP createEffects(SEXP EFFECTS, Model *pModel, vector<Data *> * pGroupData,
					   const char *networkName, int effectCol,
					   int parmCol, int int1Col, int int2Col,
					   int initValCol, int typeCol, int groupCol,
					   int periodCol, int pointerCol, int rateTypeCol,
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

				if (strcmp(netType, "oneMode") == 0)
				{
					NetworkLongitudinalData * pNetwork =
						pData->pNetworkData(networkName);
					pModel->basicRateParameter(pNetwork, period, initialValue);
				}
				else if (strcmp(netType, "Behavior") == 0)
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
 *  creates the requested effects
 */

    SEXP effects(SEXP RpData, SEXP EFFECTSLIST)
    {
		vector<Data *> * pGroupData = (vector<Data *> *)
			R_ExternalPtrAddr(RpData);

		/* need to know the total number of observations in the run
		   in order to set up the rates correctly*/

        int nGroups = pGroupData->size();
        /* this should be changed when the data structures are
           changed to allow for the 'virtual' dependent variables?? */

        int totObservations = 0;
        for (int group = 0; group < nGroups; group++)
            totObservations += (*pGroupData)[group]->observationCount() - 1;

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

		getColNos(Names, &netTypeCol, &nameCol, &effectCol,
				  &parmCol, &int1Col, &int2Col, &initValCol,
				  &typeCol, &groupCol, &periodCol, &pointerCol,
				  &rateTypeCol);

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
 *  removes the objects created for the data. TODO more things to delete?
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
 *  removes the model object. TODO more things to delete?
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
		SEXP CONDVAR, SEXP CONDTARGETS, SEXP PROFILEDATA)
    {
        /* get hold of the data vector */
		vector<Data *> * pGroupData = (vector<Data *> *)
			R_ExternalPtrAddr(DATAPTR);
		int nGroups = pGroupData->size();
		int totObservations = 0;
		for (int group = 0; group < nGroups; group++)
			totObservations += (*pGroupData)[group]->observationCount() - 1;

        /* get hold of the model object */
        Model * pModel = (Model *) R_ExternalPtrAddr(MODELPTR);

        if (!isNull(CONDVAR))
        {
			int *change = INTEGER(CONDTARGETS);
            pModel->conditional(true);
			pModel->conditionalDependentVariable(CHAR(STRING_ELT(CONDVAR,0)));
			for (int i = 0; i < totObservations; i++)
			{
				pModel->addTargetChange(change[i]);
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
					pNetworkData->maxDegree(INTEGER(VECTOR_ELT(MAXDEGREE, i))[0]);
				}
			}
		}
		// print out Data for profiling
		if (asInteger(PROFILEDATA))
		{
			printOutData((*pGroupData)[0]);
		}
		return R_NilValue;

    }
/**
 *  retrieves the values of the statistics for each of the effects,
 *  for one period. The call will relate to one group only, although all effects
 *  are the same apart from the basic rates.
 */
	void getStatistics(SEXP EFFECTSLIST, StatisticCalculator * pCalculator,
					   int period, int group, Data *pData, EpochSimulation *
					   pEpochSimulation,
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

		getColNos(Names, &netTypeCol, &nameCol, &effectCol,
				  &parmCol, &int1Col, &int2Col, &initValCol,
				  &typeCol, &groupCol, &periodCol, &pointerCol,
				  &rateTypeCol);


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
							Rprintf("here\n");
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

        int totObservations = 0;
        for (int group = 0; group < nGroups; group++)
            totObservations += (*pGroupData)[group]->observationCount() - 1;

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
					EpochSimulation  Simulation(pData, pModel);
					Simulation.initialize(period + 1);
					//State State (pData, period + 1);
				State State (&Simulation);
				StatisticCalculator Calculator (pData, pModel, &State,
					period);
 				vector<double> statistic(nEffects);
 				vector<double> score(nEffects); /* not used */

// 				getStatistics(EFFECTSLIST, &Calculator, period,
// 					group, pData, (EpochSimulation *) 0,
// 					&statistic, &score);
  				getStatistics(EFFECTSLIST, &Calculator, period,
					group, pData, &Simulation,
					&statistic, &score);
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
/************************************************************************
*******************************************************************************
*************************************************************************/

/**
 *  Does one simulation
 */

    SEXP model(SEXP DERIV, SEXP DATAPTR, SEXP SEEDS,
			   SEXP FROMFINITEDIFF, SEXP MODELPTR, SEXP EFFECTSLIST,
		SEXP THETA, SEXP RANDOMSEED2, SEXP RETURNDEPS)
    {
		SEXP NEWRANDOMSEED; /* for parallel testing only */
		PROTECT(NEWRANDOMSEED = duplicate(RANDOMSEED2));

		/* create a simulation and return the observed statistics and scores */

        /* get hold of the data vector */
		vector<Data *> * pGroupData = (vector<Data *> *)
			R_ExternalPtrAddr(DATAPTR);

       /* get hold of the model object */
        Model * pModel = (Model *) R_ExternalPtrAddr(MODELPTR);

        int nGroups = pGroupData->size();
        /* the group loop should be removed when the data structures are
           changed to allow for the 'virtual' dependent variables */


		int totObservations = 0;
        for (int group = 0; group < nGroups; group++)
            totObservations += (*pGroupData)[group]->observationCount() - 1;

		int fromFiniteDiff = asInteger(FROMFINITEDIFF);

		int returnDependents = asInteger(RETURNDEPS);

		int deriv = asInteger(DERIV);

		/* set the deriv flag on the model */
		pModel->needScores(deriv);

		/* update the parameters */
		updateParameters(EFFECTSLIST, THETA, pGroupData, pModel);

		int dim = 0;
		for (int i = 0; i < length(EFFECTSLIST); i++)
		{
			dim += length(VECTOR_ELT(VECTOR_ELT(EFFECTSLIST, i), 0));
		}

		/* get the random seed from R into memory */
		GetRNGstate();

        /* fra will contain the simulated statistics and must be initialised
           to 0. Use ifra to reduce function evaluations. */
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

        /* ans will be the return value */
        SEXP ans;
        PROTECT(ans = allocVector(VECSXP, 6));

		/* nets will be the returned simulated networks */
		SEXP nets;
        PROTECT(nets = allocVector(VECSXP, 2));

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
        /* group loop here */
        for (int group = 0; group < nGroups; group++)
        {
			/* random states need store (not fromFiniteDiff)
			   and restore (fromFiniteDiff) for each period
			   within each  group */
            SEXP seeds = 0;
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
				if (!isNull(RANDOMSEED2))
				{
					defineVar(rs, NEWRANDOMSEED, R_GlobalEnv);
					GetRNGstate();
					//double dummy =
					nextDouble();
					PutRNGstate();
					NEWRANDOMSEED = findVar(rs, R_GlobalEnv);
//  					Rprintf("%d %d %d %d\n",INTEGER(NEWRANDOMSEED)[0],
// 							INTEGER(NEWRANDOMSEED)[1],
// 							INTEGER(NEWRANDOMSEED)[2],
//  							INTEGER(NEWRANDOMSEED)[3]);
				}
				else
				{
					if (fromFiniteDiff) /* restore state */
					{
						defineVar(rs, VECTOR_ELT(seeds, period), R_GlobalEnv);
						GetRNGstate();
					}
					else /* save state */
					{
						PutRNGstate();
						SET_VECTOR_ELT(VECTOR_ELT(seedstore, group),
									   period, findVar(rs, R_GlobalEnv));
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
// 				if (deriv)
// 				{
//                     rscores[period + iii] = score;
// 					rscores[totObservations + iii] = score[1];
// 					rscores[totObservations + 1 + iii] = score[2];
// 					rscores[totObservations + 2 + iii] = score[3];
//                 }
				if (pModel->conditional())
				{
                    rntim[periodFromStart - 1] = pEpochSimulation->time();
				}
				// get simulated network
				if (returnDependents)
				{
					const DependentVariable * thisnet =
						pEpochSimulation->rVariables()[0];
					const NetworkVariable * thisnetv =
						(const NetworkVariable * ) thisnet;
					const Network * thisn = thisnetv->pNetwork();
					SEXP thisedge = getEdgeList(*thisn);
					SET_VECTOR_ELT(nets, period, thisedge);
				}
			} /* end of period */
			delete pEpochSimulation;
   } /* end of group */
       /* send the .Random.seed back to R */
        PutRNGstate();

        /* set up the return object */

        if (!fromFiniteDiff)
        {
            SET_VECTOR_ELT(ans, 2, seedstore);
        }
		if (deriv)
        {
            SET_VECTOR_ELT(ans, 1, scores);
 		SET_VECTOR_ELT(ans, 5, nets);
       }
      SET_VECTOR_ELT(ans, 0, fra);
		SET_VECTOR_ELT(ans, 3, ntim);
// 			Rprintf("r2 %d %d %d %d\n",INTEGER(NEWRANDOMSEED)[0],
// 					INTEGER(NEWRANDOMSEED)[1],
// 					INTEGER(NEWRANDOMSEED)[2],
// 					INTEGER(NEWRANDOMSEED)[3]);
		if (!isNull(RANDOMSEED2))
		{
			SET_VECTOR_ELT(ans, 4, NEWRANDOMSEED);
		}
        UNPROTECT(8);
        return(ans);
    }


}


