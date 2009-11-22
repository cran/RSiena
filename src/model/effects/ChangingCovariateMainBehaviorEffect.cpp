/******************************************************************************
 * SIENA: Simulation Investigation for Empirical Network Analysis
 *
 * Web: http://www.stats.ox.ac.uk/~snijders/siena/
 *
 * File: ChangingCovariateMainBehaviorEffect.cpp
 *
 * Description: This file contains the implementation of the
 * ChangingCovariateMainBehaviorEffect class.
 *****************************************************************************/

#include <cmath>

#include "ChangingCovariateMainBehaviorEffect.h"
#include "data/ChangingCovariate.h"
#include "model/variables/BehaviorVariable.h"
#include <R.h>
namespace siena
{

/**
 * Constructor.
 */
ChangingCovariateMainBehaviorEffect::ChangingCovariateMainBehaviorEffect(
	const EffectInfo * pEffectInfo) :
		ChangingCovariateBehaviorEffect(pEffectInfo)
{
}


/**
 * Calculates the change in the statistic corresponding to this effect if
 * the given actor would change his behavior by the given amount.
 */
double ChangingCovariateMainBehaviorEffect::calculateChangeContribution(int
	actor, int difference) const
{
	// The formula for the effect:
	// s_i(x) = v_i *  c_i
	// We need to calculate the change delta in s_i(x), if we changed
	// v_i to v_i + d (d being the given amount of change in v_i).
	// This is d * c_i.

	return difference * this->pCovariate()->value(actor, this->period());
}


/**
 * Returns the statistic corresponding to this effect as part of
 * the evaluation function with respect to the given behavior variable.
 */
double ChangingCovariateMainBehaviorEffect::evaluationStatistic(double * currentValues) const
{
	double statistic = 0;
	int n = this->n();

	const ChangingCovariate * pCovariate = this->pCovariate();

	for (int i = 0; i < n; i++)
	{
		statistic += currentValues[i] * pCovariate->value(i, this->period());
	}

	return statistic;
}


/**
 * Returns the statistic corresponding to this effect as part of
 * the endowment function with respect to the initial values of a
 * behavior variable and the current values.
 */
double ChangingCovariateMainBehaviorEffect::endowmentStatistic(const int * difference,
	double * currentValues) const
{
	double statistic = 0;
	int n = this->n();

	const ChangingCovariate * pCovariate = this->pCovariate();

	for (int i = 0; i < n; i++)
	{
		if (difference[i] > 0)
		{
			statistic +=
				currentValues[i] * pCovariate->value(i, this->period());
		}
	}

	return statistic;
}

}
