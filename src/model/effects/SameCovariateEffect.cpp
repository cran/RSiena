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
#include "data/CommonNeighborIterator.h"
#include "model/variables/NetworkVariable.h"

namespace siena
{

/**
 * Constructor.
 * @param[in] pEffectInfo the effect descriptor
 * @param[in] reciprocal indicates if only reciprocal ties have to be
 * considered
 */
SameCovariateEffect::SameCovariateEffect(const EffectInfo * pEffectInfo,
	bool reciprocal) :
	CovariateDependentNetworkEffect(pEffectInfo)
{
	this->lreciprocal = reciprocal;
}


/**
 * Calculates the contribution of a tie flip to the given actor.
 */
double SameCovariateEffect::calculateTieFlipContribution(int alter) const
{
	int change = 0;

	if (!this->lreciprocal || this->pVariable()->inTieExists(alter))
	{
		if (fabs(this->value(alter) - this->value(this->pVariable()->ego())) <
			EPSILON)
		{
			change = 1;

			if (this->pVariable()->outTieExists(alter))
			{
				// The ego would loose the tie, so the change is -1.
				change = -1;
			}
		}
	}

	return change;
}


/**
 * Detailed comment in the base class.
 */
double SameCovariateEffect::statistic(Network * pNetwork,
	Network * pSummationTieNetwork) const
{
	double statistic = 0;
	int n = pNetwork->n();

	for (int i = 0; i < n; i++)
	{
		if (!this->missing(i))
		{
			double egoValue = this->value(i);

			// TODO: This is not very elegant. If CommonNeighborIterator and
			// IncidentTieIterator had a common base class, we could join the
			// cycles below.

			if (this->lreciprocal)
			{
				CommonNeighborIterator iter(pSummationTieNetwork->outTies(i),
					pNetwork->inTies(i));

				while (iter.valid())
				{
					if (!this->missing(iter.actor()))
					{
						if (fabs(this->value(iter.actor()) - egoValue) <
							EPSILON)
						{
							statistic++;
						}
					}

					iter.next();
				}
			}
			else
			{
				for (IncidentTieIterator iter =
						pSummationTieNetwork->outTies(i);
					iter.valid();
					iter.next())
				{
					if (!this->missing(iter.actor()))
					{
						if (fabs(this->value(iter.actor()) - egoValue) <
							EPSILON)
						{
							statistic++;
						}
					}
				}
			}
		}
	}

	return statistic;
}

}
