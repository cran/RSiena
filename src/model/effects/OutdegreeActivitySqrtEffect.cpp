/******************************************************************************
 * SIENA: Simulation Investigation for Empirical Network Analysis
 * 
 * Web: http://www.stats.ox.ac.uk/~snijders/siena/
 * 
 * File: OutdegreeActivitySqrtEffect.cpp
 * 
 * Description: This file contains the implementation of the
 * OutdegreeActivitySqrtEffect class.
 *****************************************************************************/

#include <cmath>

#include "OutdegreeActivitySqrtEffect.h"
#include "utils/SqrtTable.h"
#include "data/Network.h"
#include "model/variables/NetworkVariable.h"

namespace siena
{

/**
 * Constructor.
 */
OutdegreeActivitySqrtEffect::OutdegreeActivitySqrtEffect(
	const EffectInfo * pEffectInfo) : NetworkEffect(pEffectInfo)
{
	this->lsqrtTable = SqrtTable::instance();
}

	
/**
 * Calculates the contribution of a tie flip to the given actor.
 */
double OutdegreeActivitySqrtEffect::calculateTieFlipContribution(int alter)
	const
{
	// Current out-degree
	
	int degree =
		this->pVariable()->pNetwork()->outDegree(this->pVariable()->ego());
	
	int newDegree = degree + 1;
	
	if (this->pVariable()->outTieExists(alter))
	{
		// The tie is going to be withdrawn, so the degree actually decreases
		newDegree = degree - 1;
	}
	
	return newDegree * this->lsqrtTable->sqrt(newDegree) -
		degree * this->lsqrtTable->sqrt(degree);
}


/**
 * Returns the statistic corresponding to this effect as part of
 * the evaluation function with respect to the given network.
 */
double OutdegreeActivitySqrtEffect::evaluationStatistic(Network * pNetwork)
	const
{
	double statistic = 0;

	// Iterate over senders and sum up their outdegrees to the power of 1.5
	
	for (int i = 0; i < pNetwork->n(); i++)
	{
		statistic +=
			pNetwork->outDegree(i) *
				this->lsqrtTable->sqrt(pNetwork->outDegree(i));
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
double OutdegreeActivitySqrtEffect::endowmentStatistic(
	Network * pInitialNetwork,
	Network * pLostTieNetwork) const
{
	double statistic = 0;

	// Iterate over senders and sum up the square roots of their
	// outdegrees multiplied by the number of lost outgoing ties.
	
	for (int i = 0; i < pInitialNetwork->n(); i++)
	{
		statistic +=
			this->lsqrtTable->sqrt(pInitialNetwork->outDegree(i)) *
				pLostTieNetwork->outDegree(i);
	}
	
	return statistic;
}

}
