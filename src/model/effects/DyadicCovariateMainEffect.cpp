/******************************************************************************
 * SIENA: Simulation Investigation for Empirical Network Analysis
 *
 * Web: http://www.stats.ox.ac.uk/~snijders/siena/
 *
 * File: DyadicCovariateMainEffect.cpp
 *
 * Description: This file contains the implementation of the
 * DyadicCovariateMainEffect class.
 *****************************************************************************/
#include <R.h>
#include "DyadicCovariateMainEffect.h"
#include "network/Network.h"
#include "network/TieIterator.h"
#include "model/variables/NetworkVariable.h"

namespace siena
{

/**
 * Constructor.
 */
DyadicCovariateMainEffect::DyadicCovariateMainEffect(
	const EffectInfo * pEffectInfo) :
		DyadicCovariateDependentNetworkEffect(pEffectInfo)
{
}


/**
 * Calculates the contribution of a tie flip to the given actor.
 */
double DyadicCovariateMainEffect::calculateContribution(int alter) const
{
	double change = 0;
	int ego = this->ego();

	if (!this->missing(ego, alter))
	{
		change = this->value(ego, alter);
	}

	return change;
}


/**
 * See the base class for a detailed comment.
 */
double DyadicCovariateMainEffect::statistic(
	const Network * pSummationTieNetwork) const
{
	double statistic = 0;

	for (TieIterator iter = pSummationTieNetwork->ties();
		iter.valid();
		iter.next())
	{
		if (!this->missing(iter.ego(), iter.alter()))
		{
			statistic += this->value(iter.ego(), iter.alter());
		}
	}

	return statistic;
}

}
