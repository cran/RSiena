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
	this->lpBehaviorData = 0;
	this->lvalues = 0;
}


/**
 * Initializes this effect.
 * @param[in] pData the observed data
 * @param[in] pState the current state of the dependent variables
 * @param[in] period the period of interest
 * @param[in] pCache the cache object to be used to speed up calculations
 */
void CovariateDependentNetworkEffect::initialize(const Data * pData,
	State * pState,
	int period,
	Cache * pCache)
{
	NetworkEffect::initialize(pData, pState, period, pCache);
	string name = this->pEffectInfo()->interactionName1();

	this->lpConstantCovariate = pData->pConstantCovariate(name);
	this->lpChangingCovariate = pData->pChangingCovariate(name);
	this->lpBehaviorData = pData->pBehaviorData(name);
	this->lvalues = pState->behaviorValues(name);

	if (!this->lpConstantCovariate &&
		!this->lpChangingCovariate &&
		!(this->lpBehaviorData && this->lvalues))
	{
		throw logic_error("Covariate or dependent behavior variable '" +
			name +
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
	else
	{
		value = this->lvalues[i] - this->lpBehaviorData->overallMean();
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
	else
	{
		similarity =
			this->lpBehaviorData->similarity(this->lvalues[i],
				this->lvalues[j]);
	}

	return similarity;
}

/**
 * Returns the constant covariate associated with this effect.
 */
ConstantCovariate * CovariateDependentNetworkEffect::pConstantCovariate() const
{
		return this->lpConstantCovariate;

}

/**
 * Returns the changing covariate associated with this effect.
 */
ChangingCovariate * CovariateDependentNetworkEffect::pChangingCovariate() const
{
		return this->lpChangingCovariate;

}

/**
 * Returns the changing covariate associated with this effect.
 */
BehaviorLongitudinalData * CovariateDependentNetworkEffect::pBehaviorData() const
{
		return this->lpBehaviorData;

}
}
