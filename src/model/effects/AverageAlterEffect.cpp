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
#include "data/Network.h"
#include "data/IncidentTieIterator.h"

#include "model/variables/NetworkVariable.h"
#include "model/variables/BehaviorVariable.h"
#include <R.h>
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
	int difference) const
{
	double contribution = 0;
	Network * pNetwork = this->pNetworkVariable()->pNetwork();

	if (pNetwork->outDegree(actor) > 0)
	{
		// The formula for the effect:
		// s_i(x) = v_i * avg(v_j) over all neighbors j of i.
		// We need to calculate the change delta in s_i(x), if we changed
		// v_i to v_i + d (d being the given amount of change in v_i).
		// This is d * avg(v_j), the average being taken over all neighbors
		// of i. This is what is calculated below.

		double averageAlterChange = 0;

		for (IncidentTieIterator iter = pNetwork->outTies(actor);
			iter.valid();
			iter.next())
		{
			double alterValue = this->pVariable()->centeredValue(iter.actor());
			averageAlterChange += alterValue;
		}

		contribution = difference * averageAlterChange /
			pNetwork->outDegree(actor);
	}

	return contribution;
}

/**
 * Returns the statistic corresponding to this effect as part of
 * the evaluation function with respect to the given behavior variable.
 */
double AverageAlterEffect::evaluationStatistic(double * currentValues) const
{
	double statistic = 0;
	int n = this->pVariable()->n();

	for (int i = 0; i < n; i++)
	{
		if (this->pNetworkVariable()->pPredictorNetwork()->outDegree(i))
		{
			double thisStatistic = 0;
			for (IncidentTieIterator iter=this->pNetworkVariable()->
					 pPredictorNetwork()->outTies(i);
				 iter.valid();
				 iter.next())
			{
				double alterValue = currentValues[iter.actor()];
				thisStatistic += iter.value() *	alterValue;
			}
			thisStatistic *= currentValues[i] /
				this->pNetworkVariable()->pPredictorNetwork()->outDegree(i);
			statistic += thisStatistic;
		}
	}
	return statistic;
}


/**
 * Returns the statistic corresponding to this effect as part of
 * the endowment function with respect to the initial values of a
 * behavior variable and the current values.
 */
double AverageAlterEffect::endowmentStatistic(const int * difference,
	double * currentValues) const
{
	double statistic = 0;
	int n = this->pVariable()->n();

	for (int i = n-1; i > -1; i--)
	{
		if (difference[i] > 0)
		{
			if (this->pNetworkVariable()->pPredictorNetwork()->outDegree(i))
			{
				double thisStatistic = 0;
				double previousStatistic = 0;
				for (IncidentTieIterator iter=this->pNetworkVariable()->
						 pPredictorNetwork()->outTies(i);
					 iter.valid();
					 iter.next())
				{
					double alterValue = currentValues[iter.actor()];
					double alterPreviousValue = currentValues[iter.actor()];
					// +		difference[iter.actor()];
					thisStatistic += iter.value() * alterValue;
					previousStatistic += iter.value() * alterPreviousValue;
				}
				thisStatistic *= currentValues[i] /
					this->pNetworkVariable()->pPredictorNetwork()->outDegree(i);
				previousStatistic *= (currentValues[i] + difference[i])/
					this->pNetworkVariable()->pPredictorNetwork()->outDegree(i);
				statistic +=  thisStatistic - previousStatistic;
			}
		}
	}
	return statistic;
}

}
