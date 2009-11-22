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
#include "network/Network.h"
#include "network/IncidentTieIterator.h"
#include "network/CommonNeighborIterator.h"
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
double SameCovariateEffect::calculateContribution(int alter) const
{
	int change = 0;

	if (!this->lreciprocal || this->inTieExists(alter))
	{
		if (fabs(this->value(alter) - this->value(this->ego())) < EPSILON)
		{
			change = 1;
		}
	}

	return change;
}


/**
 * Detailed comment in the base class.
 */
double SameCovariateEffect::statistic(const Network * pSummationTieNetwork)
	const
{
	double statistic = 0;
	const Network * pNetwork = this->pNetwork();
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
