/******************************************************************************
 * SIENA: Simulation Investigation for Empirical Network Analysis
 *
 * Web: http://www.stats.ox.ac.uk/~snijders/siena/
 *
 * File: NetworkDependentBehaviorEffect.cpp
 *
 * Description: This file contains the implementation of the
 * NetworkDependentBehaviorEffect class.
 *****************************************************************************/

#include <stdexcept>

#include "NetworkDependentBehaviorEffect.h"
#include "model/EpochSimulation.h"
#include "model/EffectInfo.h"
#include "model/variables/NetworkVariable.h"

namespace siena
{

/**
 * Constructor.
 * @param[in] pEffectInfo the descriptor object of the effect
 */
NetworkDependentBehaviorEffect::NetworkDependentBehaviorEffect(
	const EffectInfo * pEffectInfo) : BehaviorEffect(pEffectInfo)
{
}


/**
 * Destructor.
 */
NetworkDependentBehaviorEffect::~NetworkDependentBehaviorEffect()
{
}


/**
 * Initializes this effect for the use with the given epoch simulation.
 */
void NetworkDependentBehaviorEffect::initialize(EpochSimulation * pSimulation)
{
	BehaviorEffect::initialize(pSimulation);

	this->lpNetworkVariable =
		dynamic_cast<const NetworkVariable *>(
			pSimulation->pVariable(this->pEffectInfo()->interactionName1()));

	if (!this->lpNetworkVariable)
	{
		throw logic_error("Network variable '" +
			this->pEffectInfo()->interactionName1() +
			"' expected.");
	}

	// As long as we don't have behavior effects depending on two-mode
	// networks, we explicitly insist on a one-mode network.

	if (!this->lpNetworkVariable->oneModeNetwork())
	{
		throw logic_error("One-mode network expected.");
	}
}

/**
 * Initializes this effect for calculating the corresponding statistics.
 * @param[in] pData the observed data
 * @param[in] pState the current state of the dependent variables
 * @param[in] period the period of interest
 */
void NetworkDependentBehaviorEffect::initialize(const Data * pData, State * pState, int period)
{
	Effect::initialize(pData, pState, period);

	//TODO
}
}
