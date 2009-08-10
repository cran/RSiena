/******************************************************************************
 * SIENA: Simulation Investigation for Empirical Network Analysis
 *
 * Web: http://www.stats.ox.ac.uk/~snijders/siena/
 *
 * File: OutdegreeEffect.cpp
 *
 * Description: This file contains the implementation of the
 * OutdegreeEffect class.
 *****************************************************************************/

#include <cmath>
#include "OutdegreeEffect.h"
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
OutdegreeEffect::OutdegreeEffect(
	const EffectInfo * pEffectInfo) :
		NetworkDependentBehaviorEffect(pEffectInfo)
{
}


/**
 * Calculates the change in the statistic corresponding to this effect if
 * the given actor would change his behavior by the given amount.
 */
double OutdegreeEffect::calculateChangeContribution(int actor,
	int difference) const
{
	Network * pNetwork = this->pNetworkVariable()->pNetwork();

	// The formula for the effect:
	// s_i(x) = v_i * outdegree of i.
	// We need to calculate the change delta in s_i(x), if we changed
	// v_i to v_i + d (d being the given amount of change in v_i).
	// This is  d * outdegree of i. This is what is calculated below.

	return difference * pNetwork->outDegree(actor);
}

/**
 * Returns the statistic corresponding to this effect as part of
 * the evaluation function with respect to the given behavior variable.
 */
double OutdegreeEffect::evaluationStatistic(double * currentValues) const
{
	double statistic = 0;
	int n = this->pVariable()->n();

	for (int i = 0; i < n; i++)
	{
		statistic += currentValues[i] *
			this->pNetworkVariable()->pPredictorNetwork()->outDegree(i);
	}

	return statistic;
}


/**
 * Returns the statistic corresponding to this effect as part of
 * the endowment function with respect to the initial values of a
 * behavior variable and the current values.
 */
double OutdegreeEffect::endowmentStatistic(const int * difference,
	double * currentValues) const
{
	double statistic = 0;
	int n = this->pVariable()->n();

	for (int i = 0; i < n; i++)
	{
		if (difference[i] > 0)
		{
			statistic += currentValues[i] *
				this->pNetworkVariable()->pPredictorNetwork()->outDegree(i);
				//		-(currentValues[i] + difference[i]) *
				//	this->pNetworkVariable()->pPredictorNetwork()->outDegree(i);
		}
	}
	return statistic;
}

}
