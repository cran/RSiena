/******************************************************************************
 * SIENA: Simulation Investigation for Empirical Network Analysis
 *
 * Web: http://www.stats.ox.ac.uk/~snijders/siena/
 *
 * File: ConstantCovariateBehaviorEffect.cpp
 *
 * Description: This file contains the implementation of the
 * ConstantCovariateBehaviorEffect class.
 *****************************************************************************/

#include <stdexcept>

#include "data/ConstantCovariate.h"
#include "ConstantCovariateBehaviorEffect.h"
#include "model/EpochSimulation.h"
#include "model/EffectInfo.h"

namespace siena
{

/**
 * Constructor.
 * @param[in] pEffectInfo the descriptor object of the effect
 */
ConstantCovariateBehaviorEffect::ConstantCovariateBehaviorEffect(
	const EffectInfo * pEffectInfo) : BehaviorEffect(pEffectInfo)
{
}


/**
 * Destructor.
 */
ConstantCovariateBehaviorEffect::~ConstantCovariateBehaviorEffect()
{
}


/**
 * Initializes this effect.
 * @param[in] pData the observed data
 * @param[in] pState the current state of the dependent variables
 * @param[in] period the period of interest
 * @param[in] pCache the cache object to be used to speed up calculations
 */
void ConstantCovariateBehaviorEffect::initialize(const Data * pData,
	State * pState,
	int period,
	Cache * pCache)
{
	BehaviorEffect::initialize(pData, pState, period, pCache);
	string name = this->pEffectInfo()->interactionName1();

	this->lpCovariate = pData->pConstantCovariate(name);

	if (!this->lpCovariate)
	{
		throw logic_error("Constant covariate  '" + name + "' expected.");
	}
}

}
