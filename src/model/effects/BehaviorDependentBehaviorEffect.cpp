/******************************************************************************
 * SIENA: Simulation Investigation for Empirical Network Analysis
 *
 * Web: http://www.stats.ox.ac.uk/~snijders/siena/
 *
 * File: BehaviorDependentBehaviorEffect.cpp
 *
 * Description: This file contains the implementation of the
 * BehaviorDependentBehaviorEffect class.
 *****************************************************************************/

#include <stdexcept>

#include "BehaviorDependentBehaviorEffect.h"
#include "model/EpochSimulation.h"
#include "model/EffectInfo.h"
#include "model/variables/BehaviorVariable.h"

namespace siena
{

/**
 * Constructor.
 * @param[in] pEffectInfo the descriptor object of the effect
 */
BehaviorDependentBehaviorEffect::BehaviorDependentBehaviorEffect(
	const EffectInfo * pEffectInfo) : BehaviorEffect(pEffectInfo)
{
}


/**
 * Destructor.
 */
BehaviorDependentBehaviorEffect::~BehaviorDependentBehaviorEffect()
{
}


/**
 * Initializes this effect for the use with the given epoch simulation.
 */
void BehaviorDependentBehaviorEffect::initialize(EpochSimulation * pSimulation)
{
	BehaviorEffect::initialize(pSimulation);

	this->lpBehaviorVariable =
		dynamic_cast<const BehaviorVariable *>(
			pSimulation->pVariable(this->pEffectInfo()->interactionName1()));

	if (!this->lpBehaviorVariable)
	{
		throw logic_error("Behavior variable  '" +
			this->pEffectInfo()->interactionName1() +
			"' expected.");
	}
}

/**
 * Initializes this effect for calculating the corresponding statistics.
 * @param[in] pData the observed data
 * @param[in] pState the current state of the dependent variables
 * @param[in] period the period of interest
 */
void BehaviorDependentBehaviorEffect::initialize(const Data * pData, State * pState, int period)
{
	Effect::initialize(pData, pState, period);

	//TODO
}

}
