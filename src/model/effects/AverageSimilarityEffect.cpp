/******************************************************************************
 * SIENA: Simulation Investigation for Empirical Network Analysis
 *
 * Web: http://www.stats.ox.ac.uk/~snijders/siena/
 *
 * File: AverageSimilarityEffect.cpp
 *
 * Description: This file contains the implementation of the
 * AverageSimilarityEffect class.
 *****************************************************************************/

#include <cmath>
#include "AverageSimilarityEffect.h"
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
AverageSimilarityEffect::AverageSimilarityEffect(
	const EffectInfo * pEffectInfo) :
		NetworkDependentBehaviorEffect(pEffectInfo)
{
}


/**
 * Calculates the change in the statistic corresponding to this effect if
 * the given actor would change his behavior by the given amount.
 */
double AverageSimilarityEffect::calculateChangeContribution(int actor,
	int difference) const
{
	double contribution = 0;
	const Network * pNetwork = this->pNetwork();

	if (pNetwork->outDegree(actor) > 0)
	{
		// The formula for the effect:
		// s_i(x) = avg(sim(v_i, v_j) - centeringConstant) over all neighbors
		// j of i.
		// sim(v_i, v_j) = 1.0 - |v_i - v_j| / observedRange
		// We need to calculate the change delta in s_i(x), if we changed
		// v_i to v_i + d (d being the given amount of change in v_i).
		// To this end, we can disregard the centering constant and
		// compute the average change in similarity, namely,
		// avg(sim(v_i + d, v_j) - sim(v_i, v_j)) =
		// avg(1 - |v_i+d-v_j|/range - 1 + |v_i-v_j|/range) =
		// avg(|v_i-v_j| - |v_i+d-v_j|) / range,
		// the average being taken over all neighbors of i. This is what
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

		contribution = ((double) totalSimilarityChange) /
			this->range() /
			pNetwork->outDegree(actor);
	}

	return contribution;
}

/**
 * Returns the statistic corresponding to this effect as part of
 * the evaluation function with respect to the given behavior variable.
 */
double AverageSimilarityEffect::evaluationStatistic(double * currentValues) const
{
	double statistic = 0;
	int n = this->n();
	const Network * pNetwork = this->pNetwork();

	double similarityMean =  this->similarityMean();

	for (int i = 0; i < n; i++)
	{
		if (pNetwork->outDegree(i))
		{
			double thisStatistic = 0;

			for (IncidentTieIterator iter = pNetwork->outTies(i);
				 iter.valid();
				 iter.next())
			{
				double alterValue = currentValues[iter.actor()];
				double range = this->range();
				thisStatistic += iter.value() *
					(1.0 - fabs(alterValue - currentValues[i]) / range);
				thisStatistic -= similarityMean;
			}

			thisStatistic /= pNetwork->outDegree(i);
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
double AverageSimilarityEffect::endowmentStatistic(const int * difference,
	double * currentValues) const
{
	double statistic = 0;
	int n = this->n();
	const Network * pNetwork = this->pNetwork();

	double similarityMean =  this->similarityMean();

	for (int i = 0; i < n; i++)
	{
		if (difference[i] > 0)
		{
			if (pNetwork->outDegree(i))
			{
				double thisStatistic = 0;

				for (IncidentTieIterator iter = pNetwork->outTies(i);
					 iter.valid();
					 iter.next())
				{
					double alterValue = currentValues[iter.actor()];
					double range = this->range();
					thisStatistic += iter.value() *
						(1.0 - fabs(alterValue - currentValues[i]) / range);
					thisStatistic -= similarityMean;
				}

				thisStatistic /= pNetwork->outDegree(i);
				statistic += 2 * thisStatistic;

				// do the same using the difference in i's value
				// rather than current state and subtract it.
				// not sure whether this is correct.

				thisStatistic = 0;

				for (IncidentTieIterator iter = pNetwork->outTies(i);
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

				thisStatistic /= pNetwork->outDegree(i);
				statistic -= thisStatistic;
			}
		}
	}

	return statistic;
}

}
