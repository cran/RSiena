/******************************************************************************
 * SIENA: Simulation Investigation for Empirical Network Analysis
 *
 * Web: http://www.stats.ox.ac.uk/~snijders/siena/
 *
 * File: AverageAlterEffect.cpp
 *
 * Description: This file contains the implementation of the
 * AverageAlterEffect class.
 *****************************************************************************/

#include <cmath>
#include "AverageAlterEffect.h"
#include "network/Network.h"
#include "network/IncidentTieIterator.h"

#include "model/variables/NetworkVariable.h"
#include "model/variables/BehaviorVariable.h"
namespace siena
{

/**
 * Constructor.
 */
AverageAlterEffect::AverageAlterEffect(
	const EffectInfo * pEffectInfo) :
		NetworkDependentBehaviorEffect(pEffectInfo)
{
}


/**
 * Calculates the change in the statistic corresponding to this effect if
 * the given actor would change his behavior by the given amount.
 */
double AverageAlterEffect::calculateChangeContribution(int actor,
	int difference)
{
	double contribution = 0;
	const Network * pNetwork = this->pNetwork();

	if (pNetwork->outDegree(actor) > 0)
	{
		// The formula for the effect:
		// s_i(x) = v_i * avg(v_j) over all neighbors j of i.
		// We need to calculate the change delta in s_i(x), if we changed
		// v_i to v_i + d (d being the given amount of change in v_i).
		// This is d * avg(v_j), the average being taken over all neighbors
		// of i. This is what is calculated below.

		double totalAlterValue = 0;

		for (IncidentTieIterator iter = pNetwork->outTies(actor);
			iter.valid();
			iter.next())
		{
			double alterValue = this->centeredValue(iter.actor());
			totalAlterValue += alterValue;
		}

		contribution = difference * totalAlterValue /
			pNetwork->outDegree(actor);
	}

	return contribution;
}


/**
 * Returns the statistic corresponding to the given ego with respect to the
 * given values of the behavior variable.
 */
double AverageAlterEffect::egoStatistic(int i, double * currentValues)
{
	double statistic = 0;
	const Network * pNetwork = this->pNetwork();
	int neighborCount = 0;

	for (IncidentTieIterator iter = pNetwork->outTies(i);
		 iter.valid();
		 iter.next())
	{
		int j = iter.actor();

		if (!this->missing(this->period(), j) &&
			!this->missing(this->period() + 1, j))
		{
			statistic += currentValues[j];
			neighborCount++;
		}
	}

	if (neighborCount > 0)
	{
		statistic *= currentValues[i] / neighborCount;
	}

	return statistic;
}



/**
 * Returns the statistic corresponding to the given ego as part of
 * the endowment function with respect to the initial values of a
 * behavior variable and the current values.
 */
double AverageAlterEffect::egoEndowmentStatistic(int ego,
	const int * difference,
	double * currentValues)
{
	double statistic = 0;

	const Network * pNetwork = this->pNetwork();

	if (difference[ego] > 0)
	{
		if (pNetwork->outDegree(ego))
		{
			double thisStatistic = 0;
			double previousStatistic = 0;

			for (IncidentTieIterator iter = pNetwork->outTies(ego);
				 iter.valid();
				 iter.next())
			{
				double alterValue = currentValues[iter.actor()];
				double alterPreviousValue = currentValues[iter.actor()];
				// +		difference[iter.actor()];
				thisStatistic += iter.value() * alterValue;
				previousStatistic += iter.value() * alterPreviousValue;
			}

			thisStatistic *= currentValues[ego] / pNetwork->outDegree(ego);
			previousStatistic *=
				(currentValues[ego] + difference[ego]) /
				pNetwork->outDegree(ego);
			statistic = thisStatistic - previousStatistic;
		}
	}

	return statistic;
}

}
