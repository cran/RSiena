/******************************************************************************
 * SIENA: Simulation Investigation for Empirical Network Analysis
 *
 * Web: http://www.stats.ox.ac.uk/~snijders/siena/
 *
 * File: CovariateEgoEffect.cpp
 *
 * Description: This file contains the implementation of the
 * CovariateEgoEffect class.
 *****************************************************************************/

#include "CovariateEgoEffect.h"
#include "network/Network.h"
#include "model/variables/NetworkVariable.h"

namespace siena
{

/**
 * Constructor.
 */
CovariateEgoEffect::CovariateEgoEffect(const EffectInfo * pEffectInfo) :
	CovariateDependentNetworkEffect(pEffectInfo)
{
}


/**
 * Calculates the contribution of a tie flip to the given actor.
 */
double CovariateEgoEffect::calculateContribution(int alter) const
{
	return this->value(this->ego());
}


/**
 * See the base class.
 */
double CovariateEgoEffect::statistic(const Network * pSummationTieNetwork)
	const
{
	double statistic = 0;

	for (int i = 0; i < pSummationTieNetwork->n(); i++)
	{
		if (!this->missing(i))
		{
			statistic +=
				pSummationTieNetwork->outDegree(i) * this->value(i);
		}
	}

	return statistic;
}

}
