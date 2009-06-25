/******************************************************************************
 * SIENA: Simulation Investigation for Empirical Network Analysis
 * 
 * Web: http://www.stats.ox.ac.uk/~snijders/siena/
 * 
 * File: CovariateSimilarityEffect.cpp
 * 
 * Description: This file contains the implementation of the
 * CovariateSimilarityEffect class.
 *****************************************************************************/

#include "CovariateSimilarityEffect.h"
#include "data/Network.h"
#include "data/IncidentTieIterator.h"
#include "model/variables/NetworkVariable.h"

namespace siena
{

/**
 * Constructor.
 */
CovariateSimilarityEffect::CovariateSimilarityEffect(
	const EffectInfo * pEffectInfo) :
		CovariateDependentNetworkEffect(pEffectInfo)
{
}


/**
 * Calculates the contribution of a tie flip to the given actor.
 */
double CovariateSimilarityEffect::calculateTieFlipContribution(int alter) const
{
	double change =
		this->similarity(this->pVariable()->ego(), alter);
	
	if (this->pVariable()->outTieExists(alter))
	{
		// The ego would loose the tie, so we take the opposite value
		change = -change;
	}
	
	return change;
}


/**
 * Returns the statistic corresponding to this effect as part of
 * the evaluation function with respect to the given network.
 */
double CovariateSimilarityEffect::evaluationStatistic(Network * pNetwork) const
{
	double statistic = 0;

	for (int i = 0; i < pNetwork->n(); i++)
	{
		if (!this->missing(i))
		{
			for (IncidentTieIterator iter = pNetwork->outTies(i);
				iter.valid();
				iter.next())
			{
				if (!this->missing(iter.actor()))
				{
					statistic +=
						this->similarity(i, iter.actor());
				}
			}
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
double CovariateSimilarityEffect::endowmentStatistic(Network * pInitialNetwork,
	Network * pLostTieNetwork) const
{
	// This is the same as the evaluation statistic computed with respect
	// to the network of lost ties.
	
	return this->evaluationStatistic(pLostTieNetwork);
}

}
