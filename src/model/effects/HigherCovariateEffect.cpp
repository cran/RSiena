/******************************************************************************
 * SIENA: Simulation Investigation for Empirical Network Analysis
 *
 * Web: http://www.stats.ox.ac.uk/~snijders/siena/
 *
 * File: HigherCovariateEffect.cpp
 *
 * Description: This file contains the implementation of the
 * HigherCovariateEffect class.
 *****************************************************************************/

#include "HigherCovariateEffect.h"
#include "data/Network.h"
#include "data/IncidentTieIterator.h"
#include "model/variables/NetworkVariable.h"

namespace siena
{

/**
 * Constructor.
 */
HigherCovariateEffect::HigherCovariateEffect(const EffectInfo * pEffectInfo) :
	CovariateDependentNetworkEffect(pEffectInfo)
{
}


/**
 * Calculates the contribution of a tie flip to the given actor.
 */
double HigherCovariateEffect::calculateTieFlipContribution(int alter) const
{
	double change = 0;
	double egoValue = this->value(this->pVariable()->ego());
	double alterValue = this->value(alter);

	if (egoValue > alterValue)
	{
		change = 1;
	}
	else if (egoValue == alterValue)
	{
		change = 0.5;
	}

	if (this->pVariable()->outTieExists(alter))
	{
		// The ego would loose the tie, so the change is negative.
		change = -change;
	}

	return change;
}


/**
 * Detailed comment in the base class.
 */
double HigherCovariateEffect::statistic(Network * pNetwork,
	Network * pSummationTieNetwork) const
{
	double statistic = 0;
	int n = pNetwork->n();

	for (int i = 0; i < n; i++)
	{
		if (this->missing(i))
		{
			// A contribution of 0.5 per each tie.
			statistic += pSummationTieNetwork->outDegree(i) * 0.5;
		}
		else
		{
			double egoValue = this->value(i);

			for (IncidentTieIterator iter = pSummationTieNetwork->outTies(i);
				iter.valid();
				iter.next())
			{
				if (this->missing(iter.actor()))
				{
					statistic += 0.5;
				}
				else
				{
					if (egoValue > this->value(iter.actor()))
					{
						statistic++;
					}
				}
			}
		}
	}

	return statistic;
}

}
