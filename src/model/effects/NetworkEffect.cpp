/******************************************************************************
 * SIENA: Simulation Investigation for Empirical Network Analysis
 *
 * Web: http://www.stats.ox.ac.uk/~snijders/siena/
 *
 * File: NetworkEffect.cpp
 *
 * Description: This file contains the implementation of the
 * NetworkEffect class.
 *****************************************************************************/

#include <stdexcept>

#include "NetworkEffect.h"
#include "data/NetworkLongitudinalData.h"
#include "model/State.h"
#include "model/EpochSimulation.h"
#include "model/EffectInfo.h"
#include "model/variables/NetworkVariable.h"
#include "model/tables/Cache.h"
#include "model/tables/NetworkCache.h"
#include "model/tables/EgocentricConfigurationTable.h"

using namespace std;

namespace siena
{

/**
 * Constructor.
 */
NetworkEffect::NetworkEffect(const EffectInfo * pEffectInfo) :
	Effect(pEffectInfo)
{
	this->lpNetwork = 0;
	this->lpNetworkData = 0;
	this->lpNetworkCache = 0;
	this->lpTwoPathTable = 0;
	this->lpReverseTwoPathTable = 0;
	this->lpInStarTable = 0;
	this->lpOutStarTable = 0;
	this->lpCriticalInStarTable = 0;
	this->lpRRTable = 0;
	this->lpRFTable = 0;
	this->lpRBTable = 0;
	this->lpFRTable = 0;
	this->lpBRTable = 0;
}


/**
 * Initializes this effect.
 * @param[in] pData the observed data
 * @param[in] pState the current state of the dependent variables
 * @param[in] period the period of interest
 * @param[in] pCache the cache object to be used to speed up calculations
 */
void NetworkEffect::initialize(const Data * pData,
	State * pState,
	int period,
	Cache * pCache)
{
	Effect::initialize(pData, pState, period, pCache);
	string name = this->pEffectInfo()->variableName();

	this->lpNetworkData = pData->pNetworkData(name);

	if (!this->lpNetworkData)
	{
		throw logic_error("Data for network variable '" + name +"' expected.");
	}

	this->lpNetwork = pState->pNetwork(name);
	this->lpNetworkCache = pCache->pNetworkCache(this->lpNetwork);

	this->lpTwoPathTable = this->lpNetworkCache->pTwoPathTable();
	this->lpReverseTwoPathTable = this->lpNetworkCache->pReverseTwoPathTable();
	this->lpInStarTable = this->lpNetworkCache->pInStarTable();
	this->lpOutStarTable = this->lpNetworkCache->pOutStarTable();
	this->lpCriticalInStarTable = this->lpNetworkCache->pCriticalInStarTable();
	this->lpRRTable = this->lpNetworkCache->pRRTable();
	this->lpRFTable = this->lpNetworkCache->pRFTable();
	this->lpRBTable = this->lpNetworkCache->pRBTable();
	this->lpFRTable = this->lpNetworkCache->pFRTable();
	this->lpBRTable = this->lpNetworkCache->pBRTable();
}


/**
 * Does the necessary preprocessing work for calculating the tie flip
 * contributions for a specific ego. This method must be invoked before
 * calling NetworkEffect::calculateTieFlipContribution(...).
 */
void NetworkEffect::preprocessEgo(int ego)
{
	this->lego = ego;
}


/**
 * Returns if there is a tie from the current ego to the given alter.
 */
bool NetworkEffect::outTieExists(int alter) const
{
	return this->lpNetworkCache->outTieExists(alter);
}


/**
 * Returns if there is a tie from the given alter to the current ego.
 */
bool NetworkEffect::inTieExists(int alter) const
{
	return this->lpNetworkCache->inTieExists(alter);
}


/**
 * Returns the statistic corresponding to this effect as part of
 * the evaluation function.
 */
double NetworkEffect::evaluationStatistic() const
{
	return this->statistic(this->pNetwork());
}


/**
 * Returns the statistic corresponding to this effect as part of
 * the endowment function.
 */
double NetworkEffect::endowmentStatistic(Network * pLostTieNetwork) const
{
	return this->statistic(pLostTieNetwork);
}


/**
 * A convenience method for implementing statistics for both evaluation and
 * endowment function. It assumes that the statistic can be calculated by
 * iterating over ties (i,j) of a network Y and summing up some terms
 * s_{ij}(X) with respect to another network X, namely,
 * s(X,Y) = sum_{(i,j) \in Y} s_{ij}(X).
 * For evaluation function, X = Y.
 * For endowment function, X is the initial network of the period, and Y is the
 * network of ties that have been lost during the network evolution.
 */
double NetworkEffect::statistic(const Network * pSummationTieNetwork) const
{
	// Nothing in the base class.
	return 0;
}

}
