/******************************************************************************
 * SIENA: Simulation Investigation for Empirical Network Analysis
 *
 * Web: http://www.stats.ox.ac.uk/~snijders/siena/
 *
 * File: GwespFunction.cpp
 *
 * Description: This file contains the implementation of the class
 * GwespFunction.
 *****************************************************************************/
#include "R_ext/Print.h"

#include "GwespFunction.h"
#include "model/tables/NetworkCache.h"
#include "model/tables/EgocentricConfigurationTable.h"
#include <math.h>

namespace siena
{

/**
 * Constructor.
 */
GwespFunction::GwespFunction(string networkName,
	EgocentricConfigurationTable *	(NetworkCache::*pTable)() const,
	double parameter) :
	NetworkAlterFunction(networkName)
{
	this->lparameter = parameter;
	this->lweight = -0.01 * this->lparameter;
	this->lcumulativeWeight[0] = 0;
	for (int i = 1; i < MAX_STATISTIC; i++)
	{
		this->lcumulativeWeight[i] =
			this->lcumulativeWeight[i - 1] + exp(this->lweight * (double) i);
	}
	// And normalize
	for (int i = 1; i < MAX_STATISTIC; i++)
	{
		this->lcumulativeWeight[i] /= this->lcumulativeWeight[MAX_STATISTIC - 1];
	}
	this->lpTable = pTable;
}

/**
 * Initializes this function.
 * @param[in] pData the observed data
 * @param[in] pState the current state of the dependent variables
 * @param[in] period the period of interest
 * @param[in] pCache the cache object to be used to speed up calculations
 */
void GwespFunction::initialize(const Data * pData,
	State * pState,
	int period,
	Cache * pCache)
{
	NetworkAlterFunction::initialize(pData, pState, period, pCache);
	this->lpInitialisedTable = (*this->pNetworkCache().*lpTable)();
}

/**
 * Returns the value of this function for the given alter. It is assumed
 * that the function has been initialized before and pre-processed with
 * respect to a certain ego.
 */
double GwespFunction::value(int alter)
{
	int statistic = lpInitialisedTable->get(alter);
	if (statistic >= MAX_STATISTIC)
	{
		return 1.0;
	}
	else
	{
		return this->lcumulativeWeight[statistic];
	}
}


}
