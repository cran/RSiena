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
#include "data/Network.h"
#include "data/TieIterator.h"
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
double DyadicCovariateMainEffect::calculateTieFlipContribution(int alter) const
{
	double change = 0;
	int ego = this->pVariable()->ego();

	if (!this->missing(ego, alter))
	{
		change = this->value(ego, alter);
	}

	if (this->pVariable()->outTieExists(alter))
	{
		// The ego would loose the tie and consequently the covariate value
		change = -change;
	}

	return change;
}


/**
 * Returns the statistic corresponding to this effect as part of
 * the evaluation function with respect to the given network.
 */
double DyadicCovariateMainEffect::evaluationStatistic(Network * pNetwork) const
{
	double statistic = 0;

	for (TieIterator iter = pNetwork->ties(); iter.valid(); iter.next())
	{
		if (!this->missing(iter.ego(), iter.alter()))
		{
			statistic += this->value(iter.ego(), iter.alter());
		}
	}

	return statistic;
}


/**
 * Returns the statistic corresponding to this effect as part of
 * the endowment function with respect to an initial network
 * and a network of lost ties. The current network is implicit as
 * the introduced ties are not relevant for calculating
 * endowment statistics.
 */
double DyadicCovariateMainEffect::endowmentStatistic(Network * pInitialNetwork,
	Network * pLostTieNetwork) const
{
	// This is the same as the evaluation statistic computed with respect
	// to the network of lost ties.

	return this->evaluationStatistic(pLostTieNetwork);
}

}
