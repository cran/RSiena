/******************************************************************************
 * SIENA: Simulation Investigation for Empirical Network Analysis
 *
 * Web: http://www.stats.ox.ac.uk/~snijders/siena/
 *
 * File: MixedNetworkAlterFunction.cpp
 *
 * Description: This file contains the implementation of the class
 * MixedNetworkAlterFunction.
 *****************************************************************************/
#include <R_ext/Print.h>
#include "MixedNetworkAlterFunction.h"
#include "network/Network.h"
#include "model/State.h"
#include "model/tables/Cache.h"

namespace siena
{

MixedNetworkAlterFunction::MixedNetworkAlterFunction(string firstNetworkName,
	string secondNetworkName )
{
	this->lname1 = firstNetworkName;
	this->lname2 = secondNetworkName;
	this->lpFirstNetwork = 0;
	this->lpSecondNetwork = 0;
	this->lpTwoNetworkCache = 0;
}


MixedNetworkAlterFunction::~MixedNetworkAlterFunction()
{
}


/**
 * Initializes this function.
 * @param[in] pData the observed data
 * @param[in] pState the current state of the dependent variables
 * @param[in] period the period of interest
 * @param[in] pCache the cache object to be used to speed up calculations
 */
void MixedNetworkAlterFunction::initialize(const Data * pData,
	State * pState,
	int period,
	Cache * pCache)
{
	AlterFunction::initialize(pData, pState, period, pCache);
	this->lpFirstNetwork = pState->pNetwork(this->lname1);
	this->lpSecondNetwork = pState->pNetwork(this->lname2);
	this->lpTwoNetworkCache = pCache->pTwoNetworkCache(this->lpFirstNetwork,
		this->lpSecondNetwork);
}

}
