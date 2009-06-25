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
#include "model/EpochSimulation.h"
#include "model/EffectInfo.h"
#include "model/variables/NetworkVariable.h"

using namespace std;

namespace siena
{

/**
 * Constructor.
 */
NetworkEffect::NetworkEffect(const EffectInfo * pEffectInfo) :
	Effect(pEffectInfo)
{
	this->lpVariable = 0;
	this->lpNetworkData = 0;
}


/**
 * Initializes this effect for the use with the given epoch simulation.
 */
void NetworkEffect::initialize(EpochSimulation * pSimulation)
{
	Effect::initialize(pSimulation);

	this->lpVariable =
		dynamic_cast<const NetworkVariable *>(
			pSimulation->pVariable(this->pEffectInfo()->variableName()));

	if (!this->lpVariable)
	{
		throw logic_error("Network variable '" +
			this->pEffectInfo()->variableName() +
			"' expected.");
	}

	this->lpNetworkData =
		dynamic_cast<const NetworkLongitudinalData *>(
			this->lpVariable->pData());
}


/**
 * Initializes this effect for calculating the corresponding statistics.
 * @param[in] pData the observed data
 * @param[in] pState the current state of the dependent variables
 * @param[in] period the period of interest
 */
void NetworkEffect::initialize(const Data * pData, State * pState, int period)
{
	Effect::initialize(pData, pState, period);

	this->lpNetworkData =
		pData->pNetworkData(this->pEffectInfo()->variableName());

	if (!this->lpNetworkData)
	{
		throw logic_error("Data for network variable '" +
			this->pEffectInfo()->variableName() +
			"' expected.");
	}
}


/**
 * Does the necessary preprocessing work for calculating the tie flip
 * contributions for a specific ego. This method must be invoked before
 * calling NetworkEffect::calculateTieFlipContribution(...).
 */
void NetworkEffect::preprocessEgo()
{
	// Empty here. See the derived classes for something more interesting.
}


/**
 * Returns if the given configuration table is used by this effect
 * during the calculation of tie flip contributions. The method has to
 * be overriden by all effects using any of the precalculated configuration
 * tables.
 */
bool NetworkEffect::usesTable(const ConfigurationTable * pTable) const
{
	// No tables used by default.
	return false;
}


/**
 * Returns the statistic corresponding to this effect as part of
 * the evaluation function with respect to the given network.
 */
double NetworkEffect::evaluationStatistic(Network * pNetwork) const
{
	return this->statistic(pNetwork, pNetwork);
}


/**
 * Returns the statistic corresponding to this effect as part of
 * the endowment function with respect to an initial network
 * and a network of lost ties. The current network is implicit as
 * the introduced ties are not relevant for calculating
 * endowment statistics.
 */
double NetworkEffect::endowmentStatistic(Network * pInitialNetwork,
	Network * pLostTieNetwork) const
{
	return this->statistic(pInitialNetwork, pLostTieNetwork);
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
double NetworkEffect::statistic(Network * pNetwork,
	Network * pSummationTieNetwork) const
{
	// Nothing in the base class.
	return 0;
}

}
