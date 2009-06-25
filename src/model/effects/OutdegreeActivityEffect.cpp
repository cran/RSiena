/******************************************************************************
 * SIENA: Simulation Investigation for Empirical Network Analysis
 * 
 * Web: http://www.stats.ox.ac.uk/~snijders/siena/
 * 
 * File: OutdegreeActivityEffect.cpp
 * 
 * Description: This file contains the implementation of the
 * OutdegreeActivityEffect class.
 *****************************************************************************/

#include "OutdegreeActivityEffect.h"
#include "data/Network.h"
#include "model/variables/NetworkVariable.h"

namespace siena
{

/**
 * Constructor.
 */
OutdegreeActivityEffect::OutdegreeActivityEffect(
	const EffectInfo * pEffectInfo) : NetworkEffect(pEffectInfo)
{
}

	
/**
 * Calculates the contribution of a tie flip to the given actor.
 */
double OutdegreeActivityEffect::calculateTieFlipContribution(int alter) const
{
	double change = 0;
	
	// Current out-degree
	int d =	this->pVariable()->pNetwork()->outDegree(this->pVariable()->ego());
	
	if (this->pVariable()->outTieExists(alter))
	{
		// After a tie withdrawal, the new out-degree would be d-1, and
		// the new effect statistic s_i'=(d-1)^2. The current statistic
		// is s_i=d^2, so the change would be s_i'-s_i = -2d+1.

		change = -2 * d + 1;
	}
	else
	{
		// When introducing a new tie, the new out-degree would be d+1, and
		// the new effect statistic s_i'=(d+1)^2. The current statistic
		// is s_i=d^2, so the change would be s_i'-s_i = 2d+1.
		
		change = 2 * d + 1;
	}
	
	return change;
}


/**
 * Returns the statistic corresponding to this effect as part of
 * the evaluation function with respect to the given network.
 */
double OutdegreeActivityEffect::evaluationStatistic(Network * pNetwork) const
{
	int statistic = 0;

	// Iterate over senders and sum up their squared outdegrees
	
	for (int i = 0; i < pNetwork->n(); i++)
	{
		statistic += pNetwork->outDegree(i) * pNetwork->outDegree(i);
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
double OutdegreeActivityEffect::endowmentStatistic(Network * pInitialNetwork,
	Network * pLostTieNetwork) const
{
	int statistic = 0;

	// Iterate over senders and sum up their outdegrees multiplied by
	// the number of lost outgoing ties.
	
	for (int i = 0; i < pInitialNetwork->n(); i++)
	{
		statistic +=
			pInitialNetwork->outDegree(i) * pLostTieNetwork->outDegree(i);
	}
	
	return statistic;
}

}
