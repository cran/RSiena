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

#include "BehaviorEffect.h"
#include "data/Data.h"
#include "data/BehaviorLongitudinalData.h"
#include "model/State.h"
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
void BehaviorEffect::initialize(const Data * pData,
	State * pState,
	int period,
	Cache * pCache)
{
	Effect::initialize(pData, pState, period, pCache);

	string name = this->pEffectInfo()->variableName();

	this->lpBehaviorData = pData->pBehaviorData(name);

	if (!this->lpBehaviorData)
	{
		throw logic_error(
			"Data for behavior variable '" + name +"' expected.");
	}

	this->lvalues = pState->behaviorValues(name);
}


/**
 * Returns the number of actors of the behavior variable associated with this
 * effect.
 */
int BehaviorEffect::n() const
{
	return this->lpBehaviorData->n();
}


/**
 * Returns the behavior value of the given actor.
 */
int BehaviorEffect::value(int actor) const
{
	return this->lvalues[actor];
}


/**
 * Returns the behavior value of the given actor centered around the
 * overall mean of all observed values.
 */
double BehaviorEffect::centeredValue(int actor) const
{
	return this->lvalues[actor] - this->lpBehaviorData->overallMean();
}


/**
 * Returns if the value of the behavioral variable is missing for the given
 * actor at the specified observation.
 */
bool BehaviorEffect::missing(int observation, int actor) const
{
	return this->lpBehaviorData->missing(observation, actor);
}


/**
 * Returns the observed range of the respective behavior variable.
 */
double BehaviorEffect::range() const
{
	return this->lpBehaviorData->range();
}


/**
 * Returns the centered similarity for the given values defined as
 * 1 - |a - b| / range - similarityMean.
 */
double BehaviorEffect::similarity(double a, double b) const
{
	return this->lpBehaviorData->similarity(a, b);
}


/**
 * Returns the similarity mean value over all observations.
 */
double BehaviorEffect::similarityMean() const
{
	return this->lpBehaviorData->similarityMean();
}


/**
 * Returns the statistic corresponding to this effect as part of
 * the evaluation function with respect to the given values of
 * the behavior variable.
 */
double BehaviorEffect::evaluationStatistic(double * currentValues)
{
	double statistic = 0;
	int n = this->n();

	for (int i = 0; i < n; i++)
	{
		this->preprocessEgo(i);
		if (!this->missing(this->period(), i) &&
			!this->missing(this->period() + 1, i))
		{
			statistic += this->egoStatistic(i, currentValues);
		}
	}

	return statistic;
}


/**
 * Returns the statistic corresponding to the given ego with respect to the
 * given values of the behavior variable.
 */
double BehaviorEffect::egoStatistic(int ego, double * currentValues)
{
	throw runtime_error("egoStatistic not implemented for " +
		this->pEffectInfo()->effectName());
}


/**
 * Returns the statistic corresponding to this effect as part of
 * the endowment function with respect to an initial behavior
 * variable and the current state.
 */
double BehaviorEffect::endowmentStatistic(const int * difference,
	double *currentValues)
{
	throw runtime_error("endowmentStatistic not implemented for " +
		this->pEffectInfo()->effectName());
}

/**
 * Does the necessary preprocessing work for calculating the probabilities
 * for a specific ego. This method must be invoked before
 * calling BehaviorEffect::calculateChangeContribution(...).
 */
void BehaviorEffect::preprocessEgo(int ego)
{
	this->lego = ego;
}
}
