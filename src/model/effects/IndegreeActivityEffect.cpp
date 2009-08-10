/******************************************************************************
 * SIENA: Simulation Investigation for Empirical Network Analysis
 * 
 * Web: http://www.stats.ox.ac.uk/~snijders/siena/
 * 
 * File: IndegreeActivityEffect.cpp
 * 
 * Description: This file contains the implementation of the
 * IndegreeActivityEffect class.
 *****************************************************************************/

#include <cmath>

#include "IndegreeActivityEffect.h"
#include "utils/SqrtTable.h"
#include "data/Network.h"
#include "model/variables/NetworkVariable.h"

namespace siena
{

/**
 * Constructor.
 */
IndegreeActivityEffect::IndegreeActivityEffect(
	const EffectInfo * pEffectInfo, bool root) : NetworkEffect(pEffectInfo)
{
	this->lroot = root;
	this->lsqrtTable = SqrtTable::instance();
}

	
/**
 * Calculates the contribution of a tie flip to the given actor.
 */
double IndegreeActivityEffect::calculateTieFlipContribution(int alter) const
{
	double change =
		this->pVariable()->pNetwork()->inDegree(this->pVariable()->ego());
	
	if (this->lroot)
	{
		change = this->lsqrtTable->sqrt(change);
	}
	
	if (this->pVariable()->outTieExists(alter))
	{
		// If we withdraw the tie, the contribution is negative
		change = -change;
	}
	
	return change;
}


/**
 * Returns the statistic corresponding to this effect as part of
 * the evaluation function with respect to the given network.
 */
double IndegreeActivityEffect::evaluationStatistic(Network * pNetwork) const
{
	// A one-mode network is assumed
	
	double statistic = 0;

	// Iterate over all actors and sum up their indegrees x outdegree
	// products
	
	for (int i = 0; i < pNetwork->n(); i++)
	{
		double indegreeContribution = pNetwork->inDegree(i);
		
		if (this->lroot)
		{
			indegreeContribution =
				this->lsqrtTable->sqrt(indegreeContribution);
		}
		
		statistic += indegreeContribution * pNetwork->outDegree(i);
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
double IndegreeActivityEffect::endowmentStatistic(Network * pInitialNetwork,
	Network * pLostTieNetwork) const
{
	// A one-mode network is assumed

	double statistic = 0;

	// Iterate over all actors and sum up their indegrees multiplied
	// by the number of lost outgoing ties.
	
	for (int i = 0; i < pInitialNetwork->n(); i++)
	{
		double indegreeContribution = pInitialNetwork->inDegree(i);
		
		if (this->lroot)
		{
			indegreeContribution =
				this->lsqrtTable->sqrt(indegreeContribution);
		}
		
		statistic += indegreeContribution * pLostTieNetwork->outDegree(i);
	}
	
	return statistic;
}

}
