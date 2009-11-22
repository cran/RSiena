/******************************************************************************
 * SIENA: Simulation Investigation for Empirical Network Analysis
 *
 * Web: http://www.stats.ox.ac.uk/~snijders/siena/
 *
 * File: TotalSimilarityEffect.cpp
 *
 * Description: This file contains the implementation of the
 * TotalSimilarityEffect class.
 *****************************************************************************/

#include <cmath>
#include "TotalSimilarityEffect.h"
#include "network/Network.h"
#include "network/IncidentTieIterator.h"

#include "model/variables/NetworkVariable.h"
#include "model/variables/BehaviorVariable.h"
#include <R.h>
namespace siena
{

/**
 * Constructor.
 */
TotalSimilarityEffect::TotalSimilarityEffect(
	const EffectInfo * pEffectInfo) :
		NetworkDependentBehaviorEffect(pEffectInfo)
{
}


/**
 * Calculates the change in the statistic corresponding to this effect if
 * the given actor would change his behavior by the given amount.
 */
double TotalSimilarityEffect::calculateChangeContribution(int actor,
	int difference) const
{
	double contribution = 0;
	const Network * pNetwork = this->pNetwork();

	if (pNetwork->outDegree(actor) > 0)
	{
		// The formula for the effect:
		// s_i(x) = sum(sim(v_i, v_j) - centeringConstant) over all neighbors
		// j of i.
		// sim(v_i, v_j) = 1.0 - |v_i - v_j| / observedRange
		// We need to calculate the change delta in s_i(x), if we changed
		// v_i to v_i + d (d being the given amount of change in v_i).
		// To this end, we can disregard the centering constant and
		// compute the total change in similarity, namely,
		// sum(sim(v_i + d, v_j) - sim(v_i, v_j)) =
		// sum(1 - |v_i+d-v_j|/range - 1 + |v_i-v_j|/range) =
		// sum(|v_i-v_j| - |v_i+d-v_j|) / range,
		// the sum being taken over all neighbors of i. This is what
		// is calculated below.

		int oldValue = this->value(actor);
		int newValue = oldValue + difference;
		int totalSimilarityChange = 0;

		for (IncidentTieIterator iter = pNetwork->outTies(actor);
			iter.valid();
			iter.next())
		{
			int alterValue = this->value(iter.actor());
			totalSimilarityChange +=
				abs(oldValue - alterValue) - abs(newValue - alterValue);
		}

		contribution = ((double) totalSimilarityChange) / this->range();
	}

	return contribution;
}

/**
 * Returns the statistic corresponding to this effect as part of
 * the evaluation function with respect to the given behavior variable.
 */
double TotalSimilarityEffect::evaluationStatistic(double * currentValues) const
{
	double statistic = 0;
	int n = this->n();

	double similarityMean =  this->similarityMean();

	for (int i = 0; i < n; i++)
	{
		if (this->pNetwork()->outDegree(i))
		{
			double thisStatistic = 0;

			for (IncidentTieIterator iter = this->pNetwork()->outTies(i);
				 iter.valid();
				 iter.next())
			{
				double alterValue = currentValues[iter.actor()];
				double range = this->range();
				thisStatistic += iter.value() *
					(1.0 - fabs(alterValue - currentValues[i]) / range);
				thisStatistic -= similarityMean;
			}

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
double TotalSimilarityEffect::endowmentStatistic(const int * difference,
	double * currentValues) const
{
	double statistic = 0;
	int n = this->n();

	double similarityMean =  this->similarityMean();

	for (int i = 0; i < n; i++)
	{
		if (difference[i] > 0)
		{
			if (this->pNetwork()->outDegree(i))
			{
				double thisStatistic = 0;

				for (IncidentTieIterator iter = this->pNetwork()->outTies(i);
					 iter.valid();
					 iter.next())
				{
					double alterValue = currentValues[iter.actor()];
					double range = this->range();
					thisStatistic += iter.value() *
						(1.0 - fabs(alterValue - currentValues[i]) / range);
					thisStatistic -= similarityMean;
				}

				statistic += 2 * thisStatistic;

				// do the same using the difference in i's value
				// rather than current state and subtract it.
				// not sure whether this is correct.

				thisStatistic = 0;

				for (IncidentTieIterator iter = this->pNetwork()->outTies(i);
					 iter.valid();
					 iter.next())
				{
					double alterValue = currentValues[iter.actor()];
					double range = this->range();
					thisStatistic += iter.value() *
						(1.0 - fabs(alterValue - (difference[i] +
								currentValues[i]))
							/ range);
					thisStatistic -= similarityMean;
				}

				statistic -= thisStatistic;
			}
		}
	}

	return statistic;
}

}
