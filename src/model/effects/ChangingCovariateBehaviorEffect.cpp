/******************************************************************************
 * SIENA: Simulation Investigation for Empirical Network Analysis
 *
 * Web: http://www.stats.ox.ac.uk/~snijders/siena/
 *
 * File: ChangingCovariateBehaviorEffect.cpp
 *
 * Description: This file contains the implementation of the
 * ChangingCovariateBehaviorEffect class.
 *****************************************************************************/

#include <stdexcept>

#include "data/ChangingCovariate.h"
#include "ChangingCovariateBehaviorEffect.h"
#include "model/EpochSimulation.h"
#include "model/EffectInfo.h"

namespace siena
{

/**
 * Constructor.
 * @param[in] pEffectInfo the descriptor object of the effect
 */
ChangingCovariateBehaviorEffect::ChangingCovariateBehaviorEffect(
	const EffectInfo * pEffectInfo) : BehaviorEffect(pEffectInfo)
{
}


/**
 * Destructor.
 */
ChangingCovariateBehaviorEffect::~ChangingCovariateBehaviorEffect()
{
}


/**
 * Initializes this effect for the use with the given epoch simulation.
 */
void ChangingCovariateBehaviorEffect::initialize(EpochSimulation * pSimulation)
{
	BehaviorEffect::initialize(pSimulation);

	this->lpCovariate =
		dynamic_cast<const ChangingCovariate *>(
			pSimulation->pData()->
			pChangingCovariate(this->pEffectInfo()->interactionName1()));

	if (!this->lpCovariate)
	{
		throw logic_error("Changing Covariate  '" +
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
void ChangingCovariateBehaviorEffect::initialize(const Data * pData, State * pState, int period)
{
	Effect::initialize(pData, pState, period);

	//TODO
}
}
