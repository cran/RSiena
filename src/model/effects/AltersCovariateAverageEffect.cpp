/******************************************************************************
 * SIENA: Simulation Investigation for Empirical Network Analysis
 *
 * Web: http://www.stats.ox.ac.uk/~snijders/siena/
 *
 * File: AltersCovariateAverageEffect.cpp
 *
 * Description: This file contains the implementation of the
 * AltersCovariateAverageEffect class.
 *****************************************************************************/

#include <stdexcept>

#include "AltersCovariateAverageEffect.h"
#include "data/Data.h"
#include "network/Network.h"
#include "network/IncidentTieIterator.h"
#include "data/ConstantCovariate.h"
#include "data/ChangingCovariate.h"
#include "data/BehaviorLongitudinalData.h"
#include "model/State.h"
#include "model/EffectInfo.h"
#include "model/variables/BehaviorVariable.h"

namespace siena
{

/**
 * Constructor.
 */
AltersCovariateAverageEffect::AltersCovariateAverageEffect(
	const EffectInfo * pEffectInfo) :
	CovariateAndNetworkBehaviorEffect(pEffectInfo)
{
}



/**
 * Calculates the change in the statistic corresponding to this effect if
 * the given actor would change his behavior by the given amount.
 */
double AltersCovariateAverageEffect::calculateChangeContribution(int actor,
	int difference)
{
	return difference * this->averageAlterValue(actor);

}

/**
 * Returns the statistic corresponding to the given ego with respect to the
 * given values of the behavior variable.
 */
double AltersCovariateAverageEffect::egoStatistic(int ego, double * currentValues)
{
	double statistic = 0;

	if (!this->missingDummy(ego))
	{
		statistic = currentValues[ego] * this->averageAlterValue(ego);
	}
	return statistic;
}

/**
 * Returns the statistic corresponding to this effect as part of
 * the endowment function with respect to the initial values of a
 * behavior variable and the current values.
 */
double AltersCovariateAverageEffect::endowmentStatistic(const int * difference,
	double * currentValues)
{
	double statistic = 0;
	int n = this->pNetwork()->n();

	for (int i = 0; i < n; i++)
	{
		if (difference[i] > 0)
		{
			statistic += currentValues[i] * this->averageAlterValue(i);
		}
	}

	return statistic;
}


}
