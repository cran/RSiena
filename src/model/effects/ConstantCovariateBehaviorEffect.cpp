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
 * Initializes this effect for the use with the given epoch simulation.
 */
void ConstantCovariateBehaviorEffect::initialize(EpochSimulation * pSimulation)
{
	BehaviorEffect::initialize(pSimulation);

	this->lpCovariate =
		dynamic_cast<const ConstantCovariate *>(
			pSimulation->pData()->
			pConstantCovariate(this->pEffectInfo()->interactionName1()));

	if (!this->lpCovariate)
	{
		throw logic_error("Constant Covariate  '" +
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
void ConstantCovariateBehaviorEffect::initialize(const Data * pData, State * pState, int period)
{
	Effect::initialize(pData, pState, period);

	//TODO
}
}
