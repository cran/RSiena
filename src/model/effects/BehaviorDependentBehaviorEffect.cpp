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
#include "data/BehaviorLongitudinalData.h"
#include "model/EpochSimulation.h"
#include "model/EffectInfo.h"
#include "model/variables/BehaviorVariable.h"
#include "model/State.h"

namespace siena
{

/**
 * Constructor.
 * @param[in] pEffectInfo the descriptor object of the effect
 */
BehaviorDependentBehaviorEffect::BehaviorDependentBehaviorEffect(
	const EffectInfo * pEffectInfo) : BehaviorEffect(pEffectInfo)
{
	this->linteractionValues = 0;
	this->lpInteractionBehaviorData = 0;
}


/**
 * Initializes this effect.
 * @param[in] pData the observed data
 * @param[in] pState the current state of the dependent variables
 * @param[in] period the period of interest
 * @param[in] pCache the cache object to be used to speed up calculations
 */
void BehaviorDependentBehaviorEffect::initialize(const Data * pData,
	State * pState,
	int period,
	Cache * pCache)
{
	BehaviorEffect::initialize(pData, pState, period, pCache);
	string name = this->pEffectInfo()->interactionName1();

	this->linteractionValues = pState->behaviorValues(name);
	this->lpInteractionBehaviorData = pData->pBehaviorData(name);

	if (!this->linteractionValues || !this->lpInteractionBehaviorData)
	{
		throw logic_error("Behavior variable  '" + name + "' expected.");
	}
}


/**
 * Returns the interaction value of the given actor
 * centered around the overall mean of all observed values.
 */
double BehaviorDependentBehaviorEffect::centeredInteractionValue(int actor)
	const
{
	return this->linteractionValues[actor] -
		this->lpInteractionBehaviorData->overallMean();
}

}
