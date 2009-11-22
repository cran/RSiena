/******************************************************************************
 * SIENA: Simulation Investigation for Empirical Network Analysis
 *
 * Web: http://www.stats.ox.ac.uk/~snijders/siena/
 *
 * File: ReciprocatedTwoPathFunction.cpp
 *
 * Description: This file contains the implementation of the class
 * ReciprocatedTwoPathFunction.
 *****************************************************************************/

#include "ReciprocatedTwoPathFunction.h"
#include "model/tables/NetworkCache.h"
#include "model/tables/EgocentricConfigurationTable.h"

namespace siena
{

/**
 * Constructor.
 * @param[in] networkName the name of the network variable this function is
 * associated with
 */
ReciprocatedTwoPathFunction::ReciprocatedTwoPathFunction(string networkName) :
	NetworkAlterFunction(networkName)
{
	this->lpTable = 0;
}


/**
 * Initializes this function.
 * @param[in] pData the observed data
 * @param[in] pState the current state of the dependent variables
 * @param[in] period the period of interest
 * @param[in] pCache the cache object to be used to speed up calculations
 */
void ReciprocatedTwoPathFunction::initialize(const Data * pData,
	State * pState,
	int period,
	Cache * pCache)
{
	NetworkAlterFunction::initialize(pData, pState, period, pCache);
	this->lpTable = this->pNetworkCache()->pRRTable();
}


/**
 * Returns the value of this function for the given alter. It is assumed
 * that the function has been initialized before and pre-processed with
 * respect to a certain ego.
 */
double ReciprocatedTwoPathFunction::value(int alter)
{
	return this->lpTable->get(alter);
}


/**
 * Returns the value of this function as an integer.
 */
int ReciprocatedTwoPathFunction::intValue(int alter)
{
	return this->lpTable->get(alter);
}

}