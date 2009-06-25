/******************************************************************************
 * SIENA: Simulation Investigation for Empirical Network Analysis
 *
 * Web: http://www.stats.ox.ac.uk/~snijders/siena/
 *
 * File: BehaviorEffect.cpp
 *
 * Description: This file contains the implementation of the class
 * BehaviorEffect.
 *****************************************************************************/

#include <stdexcept>
#include <R.h>

#include "BehaviorEffect.h"
#include "model/EpochSimulation.h"
#include "model/EffectInfo.h"
#include "model/variables/BehaviorVariable.h"

namespace siena
{

/**
 * Constructor.
 */
BehaviorEffect::BehaviorEffect(const EffectInfo * pEffectInfo) :
	Effect(pEffectInfo)
{
}


/**
 * Initializes this effect for the use with the given epoch simulation.
 */
void BehaviorEffect::initialize(EpochSimulation * pSimulation)
{
	Effect::initialize(pSimulation);

	this->lpVariable =
		dynamic_cast<const BehaviorVariable *>(
			pSimulation->pVariable(this->pEffectInfo()->variableName()));

	if (!this->lpVariable)
	{
		throw logic_error("Behavior variable '" +
			this->pEffectInfo()->variableName() +
			"' expected.");
	}
}
/**
 * Initializes this effect for calculating the corresponding statistics.
 * @param[in] pData the observed data
 * @param[in] pState the current state of the dependent variables
 * @param[in] period the period of interest
 */
void BehaviorEffect::initialize(const Data * pData, State * pState, int period)
{
	Effect::initialize(pData, pState, period);

	//TODO
}

}
