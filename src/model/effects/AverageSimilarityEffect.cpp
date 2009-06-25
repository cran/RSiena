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
	Network * pNetwork = this->pNetworkVariable()->pNetwork();

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

		int oldValue = this->pVariable()->value(actor);
		int newValue = oldValue + difference;
		int totalSimilarityChange = 0;

		for (IncidentTieIterator iter = pNetwork->outTies(actor);
			iter.valid();
			iter.next())
		{
			int alterValue = this->pVariable()->value(iter.actor());
			totalSimilarityChange +=
				abs(oldValue - alterValue) - abs(newValue - alterValue);
		}

		contribution = ((double) totalSimilarityChange) /
			this->pVariable()->range() /
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
	int n = this->pVariable()->n();

	double similarityMean =  this->pVariable()->similarityMean();
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
				double range = this->pVariable()->range();
				thisStatistic += iter.value() *
					(1.0 - fabs(alterValue - currentValues[i]) / range);
				thisStatistic -= similarityMean;
			}
			thisStatistic /= this->pNetworkVariable()->pPredictorNetwork()->
				outDegree(i);
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
	int n = this->pVariable()->n();

	double similarityMean =  this->pVariable()->similarityMean();
	for (int i = 0; i < n; i++)
	{
		if (difference[i] > 0)
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
					double range = this->pVariable()->range();
					thisStatistic += iter.value() *
						(1.0 - fabs(alterValue - currentValues[i]) / range);
					thisStatistic -= similarityMean;
				}
				thisStatistic /= this->pNetworkVariable()->pPredictorNetwork()->
					outDegree(i);
				statistic += 2 * thisStatistic;


				// do the same using the difference in i's value
				// rather than current state and subtract it.
				// not sure whether this is correct.

				thisStatistic = 0;
				for (IncidentTieIterator iter=this->pNetworkVariable()->
						 pPredictorNetwork()->outTies(i);
					 iter.valid();
					 iter.next())
				{
					double alterValue = currentValues[iter.actor()];
					double range = this->pVariable()->range();
					thisStatistic += iter.value() *
						(1.0 - fabs(alterValue - (difference[i] +
								currentValues[i]))
							/ range);
					thisStatistic -= similarityMean;
				}
				thisStatistic /= this->pNetworkVariable()->pPredictorNetwork()->
					outDegree(i);
				statistic -= thisStatistic;
			}
		}
	}
	return statistic;
}

}
