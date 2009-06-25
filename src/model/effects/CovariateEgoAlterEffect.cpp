/******************************************************************************
 * SIENA: Simulation Investigation for Empirical Network Analysis
 * 
 * Web: http://www.stats.ox.ac.uk/~snijders/siena/
 * 
 * File: CovariateEgoAlterEffect.cpp
 * 
 * Description: This file contains the implementation of the
 * CovariateEgoAlterEffect class.
 *****************************************************************************/

#include "CovariateEgoAlterEffect.h"
#include "data/Network.h"
#include "data/IncidentTieIterator.h"
#include "model/variables/NetworkVariable.h"

namespace siena
{

/**
 * Constructor.
 */
CovariateEgoAlterEffect::CovariateEgoAlterEffect(
	const EffectInfo * pEffectInfo) :
		CovariateDependentNetworkEffect(pEffectInfo)
{
}


/**
 * Calculates the contribution of a tie flip to the given actor.
 */
double CovariateEgoAlterEffect::calculateTieFlipContribution(int alter) const
{
	double change =
		this->value(this->pVariable()->ego()) * this->value(alter);
	
	if (this->pVariable()->outTieExists(alter))
	{
		// The ego would loose the tie, so the contribution is the opposite
		change = -change;
	}
	
	return change;
}


/**
 * Returns the statistic corresponding to this effect as part of
 * the evaluation function with respect to the given network.
 */
double CovariateEgoAlterEffect::evaluationStatistic(Network * pNetwork) const
{
	double statistic = 0;

	for (int i = 0; i < pNetwork->n(); i++)
	{
		if (!this->missing(i))
		{
			double alterValueSum = 0;
			
			for (IncidentTieIterator iter = pNetwork->outTies(i);
				iter.valid();
				iter.next())
			{
				if (!this->missing(iter.actor()))
				{
					alterValueSum +=
						this->value(iter.actor());
				}
			}
			
			statistic += this->value(i) * alterValueSum;
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
double CovariateEgoAlterEffect::endowmentStatistic(Network * pInitialNetwork,
	Network * pLostTieNetwork) const
{
	// This is the same as the evaluation statistic computed with respect
	// to the network of lost ties.
	
	return this->evaluationStatistic(pLostTieNetwork);
}

}
