/******************************************************************************
 * SIENA: Simulation Investigation for Empirical Network Analysis
 *
 * Web: http://www.stats.ox.ac.uk/~snijders/siena/
 *
 * File: DyadicCovariateReciprocityEffect.cpp
 *
 * Description: This file contains the implementation of the
 * DyadicCovariateReciprocityEffect class.
 *****************************************************************************/

#include "DyadicCovariateReciprocityEffect.h"
#include "data/Network.h"
#include "data/CommonNeighborIterator.h"
#include "model/variables/NetworkVariable.h"

namespace siena
{

/**
 * Constructor.
 */
DyadicCovariateReciprocityEffect::DyadicCovariateReciprocityEffect(
	const EffectInfo * pEffectInfo) :
		DyadicCovariateDependentNetworkEffect(pEffectInfo)
{
}


/**
 * Calculates the contribution of a tie flip to the given actor.
 */
double DyadicCovariateReciprocityEffect::calculateTieFlipContribution(
	int alter) const
{
	double change = 0;
	int ego = this->pVariable()->ego();

	if (this->pVariable()->inTieExists(alter) && !this->missing(ego, alter))
	{
		change = this->value(ego, alter);

		if (this->pVariable()->outTieExists(alter))
		{
			// The ego would loose the tie and consequently the covariate value
			change = -change;
		}
	}

	return change;
}


/**
 * Detailed comment in the base class.
 */
double DyadicCovariateReciprocityEffect::statistic(Network * pNetwork,
	Network * pSummationTieNetwork) const
{
	double statistic = 0;
	int n = pNetwork->n();

	// The summation network is a subnetwork of the (main) network.
	// So essentially, we are iterating over ties of the summation network
	// that are reciprocated in the main network, and add up the (non-missing)
	// covariate values for these ties.

	for (int i = 0; i < n; i++)
	{
		CommonNeighborIterator iter(pSummationTieNetwork->outTies(i),
			pNetwork->inTies(i));

		while (iter.valid())
		{
			int j = iter.actor();

			if (!this->missing(i, j))
			{
				statistic += this->value(i, j);
			}

			iter.next();
		}
	}

	return statistic;
}

}
