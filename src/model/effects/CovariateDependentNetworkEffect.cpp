/******************************************************************************
 * SIENA: Simulation Investigation for Empirical Network Analysis
 *
 * Web: http://www.stats.ox.ac.uk/~snijders/siena/
 *
 * File: CovariateDependentNetworkEffect.cpp
 *
 * Description: This file contains the implementation of the
 * CovariateDependentNetworkEffect class.
 *****************************************************************************/

#include <stdexcept>
#include <R.h>

#include "CovariateDependentNetworkEffect.h"
#include "data/ConstantCovariate.h"
#include "data/ChangingCovariate.h"
#include "data/BehaviorLongitudinalData.h"
#include "model/State.h"
#include "model/EffectInfo.h"
#include "model/EpochSimulation.h"
#include "model/variables/BehaviorVariable.h"

namespace siena
{

/**
 * Constructor.
 */
CovariateDependentNetworkEffect::CovariateDependentNetworkEffect(
	const EffectInfo * pEffectInfo) : NetworkEffect(pEffectInfo)
{
	this->lpConstantCovariate = 0;
	this->lpChangingCovariate = 0;
	this->lpBehaviorVariable = 0;
	this->lpBehaviorData = 0;
}


/**
 * Destructor.
 */
CovariateDependentNetworkEffect::~CovariateDependentNetworkEffect()
{
}


/**
 * Initializes this effect for calculating contributions of changes to
 * evaluation or endowment functions.
 */
void CovariateDependentNetworkEffect::initialize(EpochSimulation * pSimulation)
{
	NetworkEffect::initialize(pSimulation);

	this->lpConstantCovariate =
		pSimulation->pData()->pConstantCovariate(
			this->pEffectInfo()->interactionName1());
	this->lpChangingCovariate =
		pSimulation->pData()->pChangingCovariate(
			this->pEffectInfo()->interactionName1());
	this->lpBehaviorVariable =
		dynamic_cast<const BehaviorVariable *>(
			pSimulation->pVariable(this->pEffectInfo()->interactionName1()));
	this->lpBehaviorData = 0;

	if (!this->lpConstantCovariate &&
		!this->lpChangingCovariate &&
		!this->lpBehaviorVariable)
	{
		throw logic_error("Covariate or dependent behavior variable '" +
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
void CovariateDependentNetworkEffect::initialize(const Data * pData,
	State * pState,
	int period)
{
	NetworkEffect::initialize(pData, pState, period);

	this->lpConstantCovariate =
		pData->pConstantCovariate(this->pEffectInfo()->interactionName1());
	this->lpChangingCovariate =
		pData->pChangingCovariate(this->pEffectInfo()->interactionName1());
	this->lpBehaviorData =
		pData->pBehaviorData(this->pEffectInfo()->interactionName1());
	this->lpBehaviorVariable = 0;

	if (!this->lpConstantCovariate &&
		!this->lpChangingCovariate &&
		!this->lpBehaviorData)
	{
		throw logic_error("Covariate or dependent behavior variable '" +
			this->pEffectInfo()->interactionName1() +
			"' expected.");
	}
}

/**
 * Returns the covariate value for the given actor.
 */
double CovariateDependentNetworkEffect::value(int i) const
{
	double value = 0;

	if (this->lpConstantCovariate)
	{
		value = this->lpConstantCovariate->value(i);
	}
	else if (this->lpChangingCovariate)
	{
		value = this->lpChangingCovariate->value(i, this->period());
	}
	else if (this->lpBehaviorVariable)
	{
		value = this->lpBehaviorVariable->centeredValue(i);
	}
	else
	{
		value = this->lpBehaviorData->value(this->period(), i) -
			this->lpBehaviorData->overallMean();
	}

	return value;
}


/**
 * Returns if the covariate value for the given actor is missing.
 */
bool CovariateDependentNetworkEffect::missing(int i) const
{
	bool missing = false;

	if (this->lpConstantCovariate)
	{
		missing = this->lpConstantCovariate->missing(i);
	}
	else if (this->lpChangingCovariate)
	{
		missing = this->lpChangingCovariate->missing(i, this->period());
	}
	else if (this->lpBehaviorVariable)
	{
		missing = this->lpBehaviorVariable->missingStartValue(i);
	}
	else
	{
		missing = this->lpBehaviorData->missing(this->period(), i);
	}
	return missing;
}


/**
 * Returns the centered similarity of the given actors.
 */
double CovariateDependentNetworkEffect::similarity(int i, int j) const
{
	double similarity = 0;

	if (this->lpConstantCovariate)
	{
		similarity =
			this->lpConstantCovariate->similarity(
				this->lpConstantCovariate->value(i),
				this->lpConstantCovariate->value(j));
	}
	else if (this->lpChangingCovariate)
	{
		similarity =
			this->lpChangingCovariate->similarity(this->value(i),
				this->value(j));
	}
	else if (this->lpBehaviorVariable)
	{
		similarity = this->lpBehaviorVariable->similarity(i, j);
	}
	else
	{
		similarity =
			this->lpBehaviorData->similarity(this->value(i), this->value(j));
	}

	return similarity;
}

}
