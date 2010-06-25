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
#include <stdexcept>
#include <iostream>
#include <fstream>
#include <valarray>
#include <vector>
#include <set>
#include <cstring>
#include <R_ext/Print.h>
#include <R_ext/Random.h>
#include <Rinternals.h>
#include <Rmath.h>
#include "data/Data.h"
#include "network/OneModeNetwork.h"
#include "network/TieIterator.h"
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
#include "utils/Utils.h"
#include "data/ActorSet.h"
#include "model/EpochSimulation.h"
#include "model/variables/DependentVariable.h"
#include "model/variables/BehaviorVariable.h"
#include "model/variables/NetworkVariable.h"
#include "model/ml/MLSimulation.h"
#include "model/ml/Chain.h"
#include "model/ml/MiniStep.h"
#include "model/ml/NetworkChange.h"
#include "model/ml/BehaviorChange.h"
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
    return hi;
}
int sample(int n)
{
    double r= runif(0, 1);
    int sel=(int) n*r ;
    return(sel);
}
SEXP getBehaviorValues(const BehaviorVariable & behavior)
{
    SEXP ans;
    int n = behavior.n();
    PROTECT(ans = allocVector(INTSXP, n));
    int *ians = INTEGER(ans);
	const int *pValues = behavior.values();
    for (int i = 0; i < n; i++)
	{
		ians[i] = pValues[i];
	}
	UNPROTECT(1);
    return ans  ;
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
    return ans;
}
SEXP getEdgeList(const Network& net)
{
    SEXP ans;
	int nties = net.tieCount();
    PROTECT(ans = allocMatrix(INTSXP, nties, 3));
    int *ians = INTEGER(ans);
    /* initialise the memory: possibly only neccesary in case of error! */
    for (int i = 0; i < nties * 3; i++)
		ians[i] = 0;
	int irow = 0;
    for (TieIterator iter=net.ties(); iter.valid(); iter.next())
    {
		ians[irow ] = iter.ego() + 1;
		ians[nties + irow] = iter.alter() + 1;
		ians[2 * nties + irow] = iter.value();
		irow ++;
    }

    UNPROTECT(1);
    return ans;
}
/**
 * utilities to access chains and ministeps
 *
 */
namespace siena
{
SEXP getMiniStepDF(const MiniStep& miniStep)
{
	SEXP MINISTEP, classname, dimnames, colnames;

	PROTECT(colnames = allocVector(STRSXP, 9));
	SET_STRING_ELT(colnames, 0, mkChar("Aspect"));
	SET_STRING_ELT(colnames, 1, mkChar("Var"));
	SET_STRING_ELT(colnames, 2, mkChar("VarName"));
	SET_STRING_ELT(colnames, 3, mkChar("Ego"));
	SET_STRING_ELT(colnames, 4, mkChar("Alter"));
	SET_STRING_ELT(colnames, 5, mkChar("Diff"));
	SET_STRING_ELT(colnames, 6, mkChar("ReciRate"));
	SET_STRING_ELT(colnames, 7, mkChar("LogOptionSetProb"));
	SET_STRING_ELT(colnames, 8, mkChar("LogChoiceProb"));

	PROTECT(MINISTEP = allocVector(VECSXP, 9));

	if (miniStep.networkMiniStep())
	{
		const NetworkChange& networkChange =
			dynamic_cast<const NetworkChange &>(miniStep);
		SET_VECTOR_ELT(MINISTEP, 0, mkString("Network"));
		SET_VECTOR_ELT(MINISTEP, 4, ScalarInteger(networkChange.alter()));
		SET_VECTOR_ELT(MINISTEP, 5, ScalarInteger(0));
	}
	else
	{
		const BehaviorChange& behaviorChange =
			dynamic_cast<const BehaviorChange &>(miniStep);
		SET_VECTOR_ELT(MINISTEP, 0, mkString("Behavior"));
		SET_VECTOR_ELT(MINISTEP, 4, ScalarInteger(0));
		SET_VECTOR_ELT(MINISTEP, 5,
			ScalarInteger(behaviorChange.difference()));
	}
	SET_VECTOR_ELT(MINISTEP, 1, ScalarInteger(miniStep.variableId()));
	SET_VECTOR_ELT(MINISTEP, 2, mkString(miniStep.variableName().c_str()));
	SET_VECTOR_ELT(MINISTEP, 3, ScalarInteger(miniStep.ego()));
	SET_VECTOR_ELT(MINISTEP, 6, ScalarReal(miniStep.reciprocalRate()));
	SET_VECTOR_ELT(MINISTEP, 7,
		ScalarReal(miniStep.logOptionSetProbability()));
	SET_VECTOR_ELT(MINISTEP, 8,
		ScalarReal(miniStep.logChoiceProbability()));

	namesgets(MINISTEP, colnames);

	PROTECT(dimnames = allocVector(INTSXP, 2));
	int * idimnames = INTEGER(dimnames);
	idimnames[0] = NA_INTEGER;
	idimnames[1] = -1;
	setAttrib(MINISTEP, R_RowNamesSymbol, dimnames);

	PROTECT(classname = allocVector(STRSXP, 1));
	SET_STRING_ELT(classname, 0, mkChar("data.frame"));
	classgets(MINISTEP, classname);

	UNPROTECT(4);
	return MINISTEP;
}
SEXP getChainDF(const Chain& chain)
{
	SEXP ans, col0, col1, col2, col3, col4, col5, col6, col7, col8,
		colnames, dimnames, classname;
	PROTECT(colnames = allocVector(STRSXP, 9));
	SET_STRING_ELT(colnames, 0, mkChar("Aspect"));
	SET_STRING_ELT(colnames, 1, mkChar("Var"));
	SET_STRING_ELT(colnames, 2, mkChar("VarName"));
	SET_STRING_ELT(colnames, 3, mkChar("Ego"));
	SET_STRING_ELT(colnames, 4, mkChar("Alter"));
	SET_STRING_ELT(colnames, 5, mkChar("Diff"));
	SET_STRING_ELT(colnames, 6, mkChar("ReciRate"));
	SET_STRING_ELT(colnames, 7, mkChar("LogOptionSetProb"));
	SET_STRING_ELT(colnames, 8, mkChar("LogChoiceProb"));

	PROTECT(ans = allocVector(VECSXP, 9));
	PROTECT(col0 = allocVector(STRSXP, chain.ministepCount() - 1));

	PROTECT(col1 = allocVector(INTSXP, chain.ministepCount() - 1));
	int * icol1 = INTEGER(col1);
	PROTECT(col2 = allocVector(STRSXP, chain.ministepCount() - 1));
	PROTECT(col3 = allocVector(INTSXP, chain.ministepCount() - 1));
	int * icol3 = INTEGER(col3);
	PROTECT(col4 = allocVector(INTSXP, chain.ministepCount() - 1));
	int * icol4 = INTEGER(col4);
	PROTECT(col5 = allocVector(INTSXP, chain.ministepCount() - 1));
	int * icol5 = INTEGER(col5);
	PROTECT(col6 = allocVector(REALSXP, chain.ministepCount() - 1));
	double * rcol6 = REAL(col6);
	PROTECT(col7 = allocVector(REALSXP, chain.ministepCount() - 1));
	double * rcol7 = REAL(col7);
	PROTECT(col8 = allocVector(REALSXP, chain.ministepCount() - 1));
	double * rcol8 = REAL(col8);

	MiniStep *pMiniStep = chain.pFirst()->pNext();
	for (int i = 0; i < chain.ministepCount() - 1; i++)
	{
		SEXP ministep;
		PROTECT(ministep = getMiniStepDF(*pMiniStep));
		//put them in the data frame
		//	PrintValue(VECTOR_ELT(ministep, 0));
		//	Rprintf("%s her at la\n",CHAR(STRING_ELT(VECTOR_ELT(ministep, 0), 0)));
		SET_STRING_ELT(col0, i, STRING_ELT(VECTOR_ELT(ministep, 0), 0));
		icol1[i] =  INTEGER(VECTOR_ELT(ministep, 1))[0];
		SET_STRING_ELT(col2, i, STRING_ELT(VECTOR_ELT(ministep, 2), 0));
		icol3[i] =  INTEGER(VECTOR_ELT(ministep, 3))[0];
		icol4[i] =  INTEGER(VECTOR_ELT(ministep, 4))[0];
		icol5[i] =  INTEGER(VECTOR_ELT(ministep, 5))[0];
		rcol6[i] =  REAL(VECTOR_ELT(ministep, 6))[0];
		rcol7[i] =  REAL(VECTOR_ELT(ministep, 7))[0];
		rcol8[i] =  REAL(VECTOR_ELT(ministep, 8))[0];
		//		PrintValue(ministep);
		pMiniStep = pMiniStep->pNext();
		UNPROTECT(1);
	}
	SET_VECTOR_ELT(ans, 0, col0);
	SET_VECTOR_ELT(ans, 1, col1);
	SET_VECTOR_ELT(ans, 2, col2);
	SET_VECTOR_ELT(ans, 3, col3);
	SET_VECTOR_ELT(ans, 4, col4);
	SET_VECTOR_ELT(ans, 5, col5);
	SET_VECTOR_ELT(ans, 6, col6);
	SET_VECTOR_ELT(ans, 7, col7);
	SET_VECTOR_ELT(ans, 8, col8);

	namesgets(ans, colnames);

	PROTECT(dimnames = allocVector(INTSXP, 2));
	int * idimnames = INTEGER(dimnames);
	idimnames[0] = NA_INTEGER;
	idimnames[1] = -chain.ministepCount() + 1;
	setAttrib(ans, R_RowNamesSymbol, dimnames);

	PROTECT(classname = allocVector(STRSXP, 1));
	SET_STRING_ELT(classname, 0, mkChar("data.frame"));
	classgets(ans, classname);
	//PrintValue(ans);
	UNPROTECT(13);
	return ans;
}


SEXP getMiniStepList(const MiniStep& miniStep, int period,
	const EpochSimulation& epochSimulation)
{
	SEXP MINISTEP, EVAL, ENDOW;
	int nEvaluationEffects, nEndowmentEffects, nRows;
	double * reval, * rendow;
	PROTECT(MINISTEP = allocVector(VECSXP, 12));
	SET_VECTOR_ELT(MINISTEP, 3, ScalarInteger(miniStep.ego()));
	nEvaluationEffects =
		epochSimulation.pModel()->
		rEvaluationEffects(miniStep.variableName()).size();
	nEndowmentEffects =
		epochSimulation.pModel()->
		rEndowmentEffects(miniStep.variableName()).size();
	if (miniStep.networkMiniStep())
	{
		const NetworkChange& networkChange =
			dynamic_cast<const NetworkChange &>(miniStep);
		SET_VECTOR_ELT(MINISTEP, 0, mkString("Network"));
		SET_VECTOR_ELT(MINISTEP, 4, ScalarInteger(networkChange.alter()));
		SET_VECTOR_ELT(MINISTEP, 5, ScalarInteger(0));
		if (epochSimulation.pModel()->needChangeContributions())
		{
			nRows = epochSimulation.pVariable(miniStep.variableName())->m();
			PROTECT(EVAL = allocMatrix(REALSXP, nRows, nEvaluationEffects));
			PROTECT(ENDOW = allocMatrix(REALSXP, nRows, nEndowmentEffects));
		reval = REAL(EVAL);
		int pos = 0;
		for (int i = 0; i < nEvaluationEffects; i++)
		{
				for (int j = 0; j < nRows; j++)
			{
					reval[pos++] = networkChange.
						evaluationEffectContribution(j, i);
// 					Rprintf(" %f\n", networkChange.
// 						evaluationEffectContribution(j, i));
	}
		}
		rendow = REAL(ENDOW);
		pos = 0;
		for (int i = 0; i < nEndowmentEffects; i++)
		{
				for (int j = 0; j < nRows; j++)
			{
					reval[pos++] = networkChange.
						endowmentEffectContribution(j, i);
				}
			}
			SET_VECTOR_ELT(MINISTEP, 9, EVAL);
			SET_VECTOR_ELT(MINISTEP, 10, ENDOW);
			UNPROTECT(2);
		}
	}
	else
	{
		SET_VECTOR_ELT(MINISTEP, 0, mkString("Behavior"));
		const BehaviorChange& behaviorChange =
			dynamic_cast<const BehaviorChange &>(miniStep);
		SET_VECTOR_ELT(MINISTEP, 4, ScalarInteger(0));
		SET_VECTOR_ELT(MINISTEP,
			5,
			ScalarInteger(behaviorChange.difference()));
		if (epochSimulation.pModel()->needChangeContributions())
		{
			nRows = 3;
			PROTECT(EVAL = allocMatrix(REALSXP, nRows, nEvaluationEffects));
			PROTECT(ENDOW = allocMatrix(REALSXP, nRows, nEndowmentEffects));
		reval = REAL(EVAL);
		int pos = 0;
		for (int i = 0; i < nEvaluationEffects; i++)
		{
				for (int j = 0; j < nRows; j++)
			{
				reval[pos++] =
						behaviorChange.evaluationEffectContribution(j, i);
	}
		}
		rendow = REAL(ENDOW);
		pos = 0;
		for (int i = 0; i < nEndowmentEffects; i++)
		{
				for (int j = 0; j < nRows; j++)
			{
					reval[pos++] =
						behaviorChange.endowmentEffectContribution(j, i);
		}
	}
	SET_VECTOR_ELT(MINISTEP, 9, EVAL);
	SET_VECTOR_ELT(MINISTEP, 10, ENDOW);
			UNPROTECT(2);
		}
	}
	SET_VECTOR_ELT(MINISTEP, 1, ScalarInteger(miniStep.variableId()));
	SET_VECTOR_ELT(MINISTEP, 11, ScalarLogical(miniStep.missing(period)));
	SET_VECTOR_ELT(MINISTEP, 2, mkString(miniStep.variableName().c_str()));
	SET_VECTOR_ELT(MINISTEP, 7,	ScalarReal(miniStep.logOptionSetProbability()));
	SET_VECTOR_ELT(MINISTEP, 8, ScalarReal(miniStep.logChoiceProbability()));
	SET_VECTOR_ELT(MINISTEP, 6, ScalarReal(miniStep.reciprocalRate()));

	UNPROTECT(1);
	return MINISTEP;
}
SEXP getChainList(const Chain& chain, const EpochSimulation& epochSimulation)
{
	SEXP ans;

	PROTECT(ans = allocVector(VECSXP, chain.ministepCount() - 1));

	MiniStep *pMiniStep = chain.pFirst()->pNext();
	for (int i = 0; i < chain.ministepCount() - 1; i++)
	{
		SET_VECTOR_ELT(ans, i, getMiniStepList(*pMiniStep, chain.period(),
				epochSimulation));
		pMiniStep = pMiniStep->pNext();
	}
	UNPROTECT(1);
	return ans;
}


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
		//pOneModeNetworkLongitudinalData->calculateProperties();
		//Rprintf("%f %f\n", pOneModeNetworkLongitudinalData->
		//	averageInDegree(), pOneModeNetworkLongitudinalData->
		//	averageOutDegree());
        UNPROTECT(6);
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

		//pNetworkLongitudinalData->calculateProperties();
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
 *  Creates all the interaction effects for one network
 */
	SEXP createInteractionEffects(SEXP EFFECTS, Model *pModel,
		vector<Data *> * pGroupData, const char *networkName,
		int effectCol, int parmCol, int int1Col, int int2Col, int initValCol,
		int typeCol, int groupCol,	int periodCol, int pointerCol, int rateTypeCol,
		int netTypeCol, int intptr1Col, int intptr2Col, int intptr3Col)
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
 *  removes the ML simulation Object.
 */

    SEXP deleteMLSimulation(SEXP RpMLSimulation)
    {
        MLSimulation * pMLSimulation =
			(MLSimulation *) R_ExternalPtrAddr(RpMLSimulation);
		delete pMLSimulation;
		//	Rprintf("deleteModel\n");
		return R_NilValue;
    }
/**
 *  sets up the model options of MAXDEGREE, CONDITIONAL
 */
    SEXP setupModelOptions(SEXP DATAPTR, SEXP MODELPTR, SEXP MAXDEGREE,
		SEXP CONDVAR, SEXP CONDTARGETS, SEXP PROFILEDATA, SEXP PARALLELRUN)
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
 *  retrieves the values of the scores for each of the effects,
 *  for one period. The call will relate to one group only, although all effects
 *  are the same apart from the basic rates.
 */
	void getScores(SEXP EFFECTSLIST, int period, int group, const Data *pData,
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
							(*rderiv)[storederiv++] = pVariable->basicRateDerivative();
						}
					}
					else
					{
						error("Non constant rate effects are not %s",
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
/**
 *  retrieves the values of the candidate parameters for each of the effects,
 *  for one period.
 */
	void getCandidatesAndShapes(SEXP EFFECTSLIST, int period, int group,
		const Data *pData,
		const MLSimulation * pMLSimulation, vector<double> * rcandidates,
		vector<int> * ibayesshapes, int nBatches)
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

		int storeCandidates = 0;

		for (int ii = 0; ii < length(EFFECTSLIST); ii++)
		{
			//	const char * networkName =
			//	CHAR(STRING_ELT(VECTOR_ELT(VECTOR_ELT(EFFECTSLIST, ii),
			//							   nameCol), 0));
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
							//TODO: make this work for more than one variable
							//	const DependentVariable * pVariable =
							//	pMLSimulation->pVariable(networkName);
							for (int batch = 0; batch < nBatches; batch++)
							{
								(*rcandidates)[storeCandidates + batch] =
									pMLSimulation->sampledBasicRates(batch);
								(*ibayesshapes)[storeCandidates + batch] =
									pMLSimulation->
									sampledBasicRatesDistributions(batch);
							}
						}
					}
					else
					{
						error("Non constant rate effects are not %s",
							"implemented for maximum likelihood.");
					}
				}
				else
				{
					const EffectInfo * pEffectInfo = (const EffectInfo *)
						R_ExternalPtrAddr(VECTOR_ELT(VECTOR_ELT(EFFECTS,
									pointerCol), i));
					for (int batch = 0; batch < nBatches; batch++)
					{
						(*rcandidates)[storeCandidates + batch] =
							pMLSimulation->candidates(pEffectInfo, batch);
						(*ibayesshapes)[storeCandidates + batch] = 0;
					}

				}
				storeCandidates += nBatches;
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
/************************************************************************
*******************************************************************************
*************************************************************************/

/**
 *  Does one simulation
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
		int totObservations = 0;
        for (int group = 0; group < nGroups; group++)
            totObservations += (*pGroupData)[group]->observationCount() - 1;

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
		pModel->needChain(returnDependents);

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
						getChainList(*(pEpochSimulation->pChain()),
							*pEpochSimulation);

					SET_VECTOR_ELT(VECTOR_ELT(chains, group), period,
						thisChain);

				}
				if (addChainToStore)
				{
					pModel->chainStore(*(pEpochSimulation->pChain()));
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
        PROTECT(ans = allocVector(VECSXP, 6));

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

    SEXP mlPeriod(SEXP DERIV, SEXP DATAPTR,
		SEXP MODELPTR, SEXP EFFECTSLIST, SEXP MLSIMPTR,
		SEXP THETA, SEXP RETURNDEPS, SEXP GROUP, SEXP PERIOD,
		SEXP NRUNMH)
    {
		/* do some MH steps and return the scores and derivs of the chain
		   at the end */

        /* get hold of the data vector */
		vector<Data *> * pGroupData = (vector<Data *> *)
			R_ExternalPtrAddr(DATAPTR);
		int group = asInteger(GROUP) - 1;

		int period = asInteger(PERIOD) - 1;

		Data * pData = (*pGroupData)[group];

		/* get hold of the model object */
        Model * pModel = (Model *) R_ExternalPtrAddr(MODELPTR);

		/* get hold of the simulation object */
        MLSimulation * pMLSimulation =
			(MLSimulation *) R_ExternalPtrAddr(MLSIMPTR);

		//	int returnDependents = asInteger(RETURNDEPS);

		int deriv = asInteger(DERIV);

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
		pMLSimulation->runEpoch(0);

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
		PROTECT(accepts = allocVector(INTSXP, 6));
		SEXP rejects;
		PROTECT(rejects = allocVector(INTSXP, 6));
		int * iaccepts = INTEGER(accepts);
		int * irejects = INTEGER(rejects);

		for (int i = 0; i < 6; i++)
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

        UNPROTECT(7);
        return(ans);
    }

    SEXP mlMakeChains(SEXP DATAPTR, SEXP MODELPTR, SEXP SIMPLERATES, SEXP PROBS)
    {
		/* set up chains and do the first few MH iters on each */

        /* get hold of the data vector */
		vector<Data *> * pGroupData = (vector<Data *> *)
			R_ExternalPtrAddr(DATAPTR);

        int nGroups = pGroupData->size();
		/* find total number of periods to process */

		int totObservations = 0;
        for (int group = 0; group < nGroups; group++)
            totObservations += (*pGroupData)[group]->observationCount() - 1;

		/* get hold of the model object */
        Model * pModel = (Model *) R_ExternalPtrAddr(MODELPTR);

		SEXP RpChains;
		PROTECT(RpChains = allocVector(VECSXP, totObservations));

		Data * pData = (*pGroupData)[0];
		/* create the simulation object or maybe an array: one period for now*/
        MLSimulation * pMLSimulation = new MLSimulation(pData, pModel);

		/* set simple rates flag */
		pMLSimulation->simpleRates(asInteger(SIMPLERATES));

		/* set probability flags */
		pMLSimulation->insertDiagonalProbability(REAL(PROBS)[0]);
		pMLSimulation->cancelDiagonalProbability(REAL(PROBS)[1]);
		pMLSimulation->permuteProbability(REAL(PROBS)[2]);
		pMLSimulation->insertPermuteProbability(REAL(PROBS)[3]);
		pMLSimulation->deletePermuteProbability(REAL(PROBS)[4]);
		pMLSimulation->randomMissingProbability(REAL(PROBS)[5]);
		pMLSimulation->missingNetworkProbability(REAL(PROBS)[6]);
		pMLSimulation->missingBehaviorProbability(REAL(PROBS)[7]);

		/* initialize the chain: this also initializes the data */
		pMLSimulation->connect(0);
		SEXP ch, ch1;
		PROTECT(ch = getChainDF(*(pMLSimulation->pChain())));
		//	PrintValue(chp);
		/* get the chain up to a reasonable length */
		GetRNGstate();
		pMLSimulation->preburnin();

		/* do some more steps */

		pMLSimulation->setUpProbabilityArray();
		int numSteps=1000;
		for (int i = 0; i < numSteps; i++)
		{
			pMLSimulation->MLStep();
		}

		/* return the pointers */
		for (int i = 0; i < totObservations; i++)
		{
			SET_VECTOR_ELT(RpChains, i,
				R_MakeExternalPtr((void *) pMLSimulation, R_NilValue,
					R_NilValue));
		}
		PROTECT(ch1 = getChainDF(*(pMLSimulation->pChain())));
		//PrintValue(chp);
		//pMLSimulation->pChain()->printConsecutiveCancelingPairs();
		SEXP ans;
		PROTECT(ans = allocVector(VECSXP, 3));
		SET_VECTOR_ELT(ans, 0, RpChains);
		SET_VECTOR_ELT(ans, 1, ch);
		SET_VECTOR_ELT(ans, 2, ch1);
		UNPROTECT(4);
		PutRNGstate();
		return ans;
	}

	SEXP reloadChain(SEXP CHAIN, SEXP MLPTR, SEXP DATAPTR, SEXP PERIOD)
	{
		/* get hold of the data vector */
		vector<Data *> * pGroupData = (vector<Data *> *)
			R_ExternalPtrAddr(DATAPTR);

		int period = asInteger(PERIOD) - 1;
		int group = 0;
		Data * pData = (*pGroupData)[group];

		/* get hold of the ml simulation object */
        MLSimulation * pMLSimulation = (MLSimulation *)
			R_ExternalPtrAddr(MLPTR);

		/* chain is a data.frame: aspect, varname, ego, alter, diff */
		SEXP ASPECT = VECTOR_ELT(CHAIN, 0);
		SEXP VARNAME = VECTOR_ELT(CHAIN, 1);
		SEXP EGO = VECTOR_ELT(CHAIN, 2);
		SEXP ALTER = VECTOR_ELT(CHAIN, 3);
		SEXP DIFF = VECTOR_ELT(CHAIN, 4);
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

		/* attach the chain to the simulation */
		pMLSimulation->pChain(pChain);

		return  R_NilValue;
	}
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

		unsigned numberChains = pModel->chainStore()->size();
		SEXP logprob;
		PROTECT(logprob = allocVector(REALSXP, numberChains));
		double * rlogprob = REAL(logprob);

		/* get the probabilities: */

		for (unsigned i = 0; i < numberChains; i++)
		{
			vector<Chain *> *pChainStore = pModel->chainStore();
			Chain * pChain = (*pChainStore)[i];
			//Rprintf("%x pchain \n", pChain);
			rlogprob[i] = pEpochSimulation->
				calculateChainProbabilities(pChain);
			//Rprintf("%d %d %f logprob\n", i, pChain->ministepCount(), rlogprob[i]);
		}
		UNPROTECT(1);
		delete pEpochSimulation;
		return logprob;
	}
	SEXP clearStoredChains(SEXP MODELPTR)
	{
		/* get hold of the model object */
        Model * pModel = (Model *) R_ExternalPtrAddr(MODELPTR);
		deallocateVector(*(pModel->chainStore()));
		return R_NilValue;
	}
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

	SEXP MCMCcycle(SEXP DATAPTR, SEXP MODELPTR, SEXP MLPTR,
		SEXP EFFECTSLIST, SEXP PERIOD, SEXP GROUP, SEXP SCALEFACTOR,
		SEXP NRUNMH, SEXP NRUNMHBATCHES)
	{
		/* get hold of the data vector */
		vector<Data *> * pGroupData = (vector<Data *> *)
			R_ExternalPtrAddr(DATAPTR);

		int period = asInteger(PERIOD) - 1;
		int group = asInteger(GROUP) - 1;;
		Data * pData = (*pGroupData)[group];

		/* get hold of the model object */
        Model * pModel = (Model *) R_ExternalPtrAddr(MODELPTR);

		/* get hold of the ml simulation object */
        MLSimulation * pMLSimulation = (MLSimulation *)
			R_ExternalPtrAddr(MLPTR);

		int nrunMH = asInteger(NRUNMH);
		int nrunMHBatches = asInteger(NRUNMHBATCHES);
		double scaleFactor = asReal(SCALEFACTOR);
		pModel->numberMHBatches(nrunMHBatches) ;

		pModel->numberMLSteps(nrunMH);

		pModel->BayesianScaleFactor(scaleFactor);

		/* set up storage for statistics for MH accept and reject */
		SEXP accepts;
		PROTECT(accepts = allocMatrix(INTSXP, nrunMHBatches, 6));
		SEXP rejects;
		PROTECT(rejects = allocMatrix(INTSXP, nrunMHBatches, 6));
		int * iaccepts = INTEGER(accepts);
		int * irejects = INTEGER(rejects);

		GetRNGstate();

		pMLSimulation->initializeMCMCcycle();
		for (int batch = 0; batch < nrunMHBatches; batch++)
		{
			pMLSimulation->runEpoch(period);
			// store the MH statistics for this batch
			for (int i = 0; i < 6; i++)
			{
				iaccepts[ i * nrunMHBatches + batch] =
					pMLSimulation->acceptances(i);
				irejects[ i * nrunMHBatches + batch] =
					pMLSimulation->rejections(i);
			}

			pMLSimulation->MHPstep();
		}

		PutRNGstate();

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

		SET_VECTOR_ELT(ans, 0, BayesAccepts);
		SET_VECTOR_ELT(ans, 1, BayesCandidates);
		SET_VECTOR_ELT(ans, 2, BayesShapes);
		SET_VECTOR_ELT(ans, 3, accepts);
		SET_VECTOR_ELT(ans, 4, rejects);
		UNPROTECT(6);
		return  ans;
	}
}
