/******************************************************************************
 * SIENA: Simulation Investigation for Empirical Network Analysis
 *
 * Web: http://www.stats.ox.ac.uk/~snijders/siena/
 *
 * File: sienaPrintUtils.cpp
 *
 * Description: This module contains various utilities, including:
 *
 *  1) printOutData: dump out the data from a single data object for
 *  profiling in sienaProfile.exe. Called from setupModelOptions if
 *  requested by user.
 *  NB. This function is probably not up to date or complete.
 *
 *  2) getMiniStepDF: create a data frame as a SEXP from a ministep
 *
 *  3) getChainDF: create a data frame as a SEXP from a chain
 *
 *  4) getMiniStepList: create a list format SEXP from a ministep
 *
 *  5) getChainList: create a list format SEXP from a chain
 *
 *  SEXP's can be printed within C using PrintValue(SEXP x)
 *****************************************************************************/
/**
 * @file
 * Dump data out for profiling.
 * Convert ministeps and chains to R format objects.
 */
#include <stdexcept>
#include <vector>
#include <cstring>
#include <iostream>
#include <fstream>
#include <Rinternals.h>
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
#include "model/Model.h"
#include "utils/Utils.h"
#include "model/EpochSimulation.h"
#include "model/variables/DependentVariable.h"
#include "model/variables/BehaviorVariable.h"
#include "model/variables/NetworkVariable.h"
#include "model/ml/Chain.h"
#include "model/ml/MiniStep.h"
#include "model/ml/NetworkChange.h"
#include "model/ml/BehaviorChange.h"

using namespace std;
using namespace siena;

//--------------------------------------
// utility functions to process groups
//-------------------------------------

/** Calculate the period number of this group and period, to access
 * correct stored chain.
 */
int periodFromStart(vector<Data *> & pGroupData, int group, int period)
{

	int periodFromStart = 0;
	for (int i = 0; i < group; i++)
	{
		periodFromStart += pGroupData[i]->observationCount() - 1;
	}
	periodFromStart += period;

	return periodFromStart;
}

/** Calculate the total number of periods in all groups, which is the dimension
 * of some returned arrays.
 */

int totalPeriods(vector<Data *> & pGroupData)
{
	int nGroups = pGroupData.size();

	int totObservations = 0;

	for (int group = 0; group < nGroups; group++)
	{
		totObservations += pGroupData[group]->observationCount() - 1;
	}
	return totObservations;

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


/**
 * print out the data for profiling with gprof
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

/** Create an R vector from a behavior variable for a single period 
 *
 */ 
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

/** Create an R matrix from a network variable for a single period 
 *
 */ 
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

/** Create an R 3 column matrix from a network variable for a single period 
 *
 */ 
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
/** Create a data frame with a single row from a ministep. (prints nicely 
 * with PrintValue)
 */
SEXP getMiniStepDF(const MiniStep& miniStep)
{
	SEXP MINISTEP, classname, dimnames, colnames;
	if (miniStep.networkMiniStep() || miniStep.behaviorMiniStep())
	{
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
	else
		return R_NilValue;
}
/** Create a data frame with chain. (prints nicely with PrintValue)
 *
 */
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

/** Create a list from a ministep. Easy to create, but prints untidily!
 *
 */
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

/** Create a list from a chain. Easy to create, but prints untidily!
 *
 */
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
