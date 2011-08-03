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
 * the endowment function with respect to the given values of
 * the behavior variable.
 */
double BehaviorEffect::endowmentStatistic(const int * difference,
	double * currentValues)
{
	double statistic = 0;
	int n = this->n();

	for (int i = 0; i < n; i++)
	{
		this->preprocessEgo(i);
		if (!this->missing(this->period(), i))
		{
			statistic += this->egoEndowmentStatistic(i, difference,
				currentValues);
		}
	}

	return statistic;
}
/**
 * Returns the statistic corresponding the given ego as part of
 * the endowment function with respect to an initial behavior
 * variable and the current state.
 */
double BehaviorEffect::egoEndowmentStatistic(int i, const int * difference,
	double *currentValues)
{
	throw runtime_error("egoEndowmentStatistic not implemented for " +
		this->pEffectInfo()->effectName());
}


/**
 * Returns the statistic corresponding to this effect as part of
 * the creation function.
 * @param[in] difference an array of differences per each actor where the
 * current value is subtracted from the initial value. Thus positive
 * differences indicate a decrease of actors' behavior, while negative values
 * indicate an increase of actors' behavior.
 * @param[in] currentValues the current state of the behavior variable
 */
double BehaviorEffect::creationStatistic(int * difference,
	double *currentValues)
{
	// Here we use a trick. The creation statistics are very similar to the
	// endowmnent statistics, but instead of summing over all actors with
	// decreasing values, we must now sum over all actors with increasing
	// values. So we just reverse the differences and call the endowment
	// statistic.

	int n = this->n();

	for (int i = 0; i < n; i++)
	{
		difference[i]  = -difference[i];
	}

	double statistic = this->endowmentStatistic(difference, currentValues);

	for (int i = 0; i < n; i++)
	{
		difference[i]  = -difference[i];
	}

	return -statistic;
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
