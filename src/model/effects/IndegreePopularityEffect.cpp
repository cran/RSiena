/******************************************************************************
 * SIENA: Simulation Investigation for Empirical Network Analysis
 * 
 * Web: http://www.stats.ox.ac.uk/~snijders/siena/
 * 
 * File: IndegreePopularityEffect.cpp
 * 
 * Description: This file contains the implementation of the
 * IndegreePopularityEffect class.
 *****************************************************************************/

#include "IndegreePopularityEffect.h"
#include "utils/SqrtTable.h"
#include "data/Network.h"
#include "model/variables/NetworkVariable.h"

namespace siena
{

/**
 * Constructor.
 */
IndegreePopularityEffect::IndegreePopularityEffect(
	const EffectInfo * pEffectInfo, bool root) : NetworkEffect(pEffectInfo)
{
	this->lroot = root;
	this->lsqrtTable = SqrtTable::instance();
}

	
/**
 * Calculates the contribution of a tie flip to the given actor.
 */
double IndegreePopularityEffect::calculateTieFlipContribution(int alter) const
{
	int degree = this->pVariable()->pNetwork()->inDegree(alter);
	
	if (!this->pVariable()->outTieExists(alter))
	{
		// The indegree will increase after introducing the tie
		degree++;
	}
	
	double change = degree;

	if (this->lroot)
	{
		change = this->lsqrtTable->sqrt(degree);
	}
	
	if (this->pVariable()->outTieExists(alter))
	{
		// If we withdraw the tie, we loose the indegree contribution
		// of the alter
		
		change = -change;
	}
	
	return change;
}


/**
 * Returns the statistic corresponding to this effect as part of
 * the evaluation function with respect to the given network.
 */
double IndegreePopularityEffect::evaluationStatistic(Network * pNetwork) const
{
	double statistic = 0;

	// Iterate over receivers and sum up their squared indegrees
	
	for (int i = 0; i < pNetwork->m(); i++)
	{
		int degree = pNetwork->inDegree(i);
		
		if (this->lroot)
		{
			statistic += degree * this->lsqrtTable->sqrt(degree);			
		}
		else
		{
			statistic += degree * degree;
		}
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
double IndegreePopularityEffect::endowmentStatistic(Network * pInitialNetwork,
	Network * pLostTieNetwork) const
{
	double statistic = 0;

	// Iterate over receivers and sum up their indegrees (or roots of
	// indegrees) multiplied by the number of lost incoming ties.
	
	for (int i = 0; i < pInitialNetwork->m(); i++)
	{
		int degree = pInitialNetwork->inDegree(i);
		
		if (this->lroot)
		{
			statistic +=
				pLostTieNetwork->inDegree(i) *
					this->lsqrtTable->sqrt(degree);			
		}
		else
		{
			statistic += pLostTieNetwork->inDegree(i) * degree;
		}
	}
	
	if (!this->lroot)
	{
		// TODO: Why do we multiply by the number of senders? Ask Tom.
		statistic *= pInitialNetwork->n();
	}
	
	return statistic;
}

}
