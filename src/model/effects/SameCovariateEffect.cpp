/******************************************************************************
 * SIENA: Simulation Investigation for Empirical Network Analysis
 * 
 * Web: http://www.stats.ox.ac.uk/~snijders/siena/
 * 
 * File: SameCovariateEffect.cpp
 * 
 * Description: This file contains the implementation of the
 * SameCovariateEffect class.
 *****************************************************************************/

#include <cmath>

#include "SameCovariateEffect.h"
#include "utils/Utils.h"
#include "data/Network.h"
#include "data/IncidentTieIterator.h"
#include "model/variables/NetworkVariable.h"

namespace siena
{

/**
 * Constructor.
 */
SameCovariateEffect::SameCovariateEffect(const EffectInfo * pEffectInfo) :
	CovariateDependentNetworkEffect(pEffectInfo)
{
}


/**
 * Calculates the contribution of a tie flip to the given actor.
 */
double SameCovariateEffect::calculateTieFlipContribution(int alter) const
{
	int change = 0;
	
	if (fabs(this->value(alter) - this->value(this->pVariable()->ego())) <
		EPSILON)
	{
		change = 1;
	}
	
	if (this->pVariable()->outTieExists(alter))
	{
		// The ego would loose the tie, so the change is 0 or -1.
		change = -change;
	}
	
	return change;
}


/**
 * Returns the statistic corresponding to this effect as part of
 * the evaluation function with respect to the given network.
 */
double SameCovariateEffect::evaluationStatistic(Network * pNetwork) const
{
	int statistic = 0;

	for (int i = 0; i < pNetwork->n(); i++)
	{
		if (!this->missing(i))
		{
			double egoValue = this->value(i);
			
			for (IncidentTieIterator iter = pNetwork->outTies(i);
				iter.valid();
				iter.next())
			{
				if (!this->missing(iter.actor()))
				{
					if (fabs(this->value(iter.actor()) - egoValue) < EPSILON)
					{
						statistic++;
					}
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
double SameCovariateEffect::endowmentStatistic(Network * pInitialNetwork,
	Network * pLostTieNetwork) const
{
	// This is the same as the evaluation statistic computed with respect
	// to the network of lost ties.
	
	return this->evaluationStatistic(pLostTieNetwork);
}

}
