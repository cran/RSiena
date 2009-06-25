/******************************************************************************
 * SIENA: Simulation Investigation for Empirical Network Analysis
 * 
 * Web: http://www.stats.ox.ac.uk/~snijders/siena/
 * 
 * File: OutdegreePopularityEffect.cpp
 * 
 * Description: This file contains the implementation of the
 * OutdegreePopularityEffect class.
 *****************************************************************************/

#include <cmath>

#include "OutdegreePopularityEffect.h"
#include "utils/SqrtTable.h"
#include "data/Network.h"
#include "model/variables/NetworkVariable.h"

namespace siena
{

/**
 * Constructor.
 */
OutdegreePopularityEffect::OutdegreePopularityEffect(
	const EffectInfo * pEffectInfo,
	bool root) : NetworkEffect(pEffectInfo)
{
	this->lroot = root;
	this->lsqrtTable = SqrtTable::instance();
}

	
/**
 * Calculates the contribution of a tie flip to the given actor.
 */
double OutdegreePopularityEffect::calculateTieFlipContribution(int alter) const
{
	int degree = this->pVariable()->pNetwork()->outDegree(alter);  
	double change = degree;
	
	if (this->lroot)
	{
		change = this->lsqrtTable->sqrt(degree);
	}
	
	if (this->pVariable()->outTieExists(alter))
	{
		// If we withdraw the tie, we loose the outdegree of the alter
		change = -change;
	}
	
	return change;
}


/**
 * Returns the statistic corresponding to this effect as part of
 * the evaluation function with respect to the given network.
 */
double OutdegreePopularityEffect::evaluationStatistic(Network * pNetwork) const
{
	// We assume a one-mode network because this effect doesn't make sense
	// for two-mode networks.
	
	double statistic = 0;

	// Iterate over all actors and sum up their products of indegrees
	// and outdegrees (or square roots of outdegrees).
	
	for (int i = 0; i < pNetwork->n(); i++)
	{
		int degree = pNetwork->outDegree(i);
		double outDegreeContribution = degree;
		
		if (this->lroot)
		{
			outDegreeContribution = this->lsqrtTable->sqrt(degree);
		}
		
		statistic += pNetwork->inDegree(i) * outDegreeContribution;
	}
	
	if (!this->lroot)
	{
		// TODO: Why do we multiply by the number of senders? Ask Tom.
		statistic *= pNetwork->n();
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
double OutdegreePopularityEffect::endowmentStatistic(Network * pInitialNetwork,
	Network * pLostTieNetwork) const
{
	// We assume a one-mode network because this effect doesn't make sense
	// for two-mode networks.

	double statistic = 0;

	// Iterate over all actors and sum up their outdegrees
	// (or square roots of outdegrees) multiplied by
	// the number of lost incoming ties.
	
	for (int i = 0; i < pInitialNetwork->n(); i++)
	{
		int degree = pInitialNetwork->outDegree(i);
		double outDegreeContribution = degree;
		
		if (this->lroot)
		{
			outDegreeContribution = this->lsqrtTable->sqrt(degree);
		}
		
		statistic += outDegreeContribution * pLostTieNetwork->inDegree(i);
	}
	
	if (!this->lroot)
	{
		// TODO: Why do we multiply by the number of senders? Ask Tom.
		statistic *= pInitialNetwork->n();
	}
	
	return statistic;
}

}
