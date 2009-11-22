/******************************************************************************
 * SIENA: Simulation Investigation for Empirical Network Analysis
 *
 * Web: http://www.stats.ox.ac.uk/~snijders/siena/
 *
 * File: CovariateSimilarityEffect.cpp
 *
 * Description: This file contains the implementation of the
 * CovariateSimilarityEffect class.
 *****************************************************************************/

#include "CovariateSimilarityEffect.h"
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
CovariateSimilarityEffect::CovariateSimilarityEffect(
	const EffectInfo * pEffectInfo,
	bool reciprocal) :
		CovariateDependentNetworkEffect(pEffectInfo)
{
	this->lreciprocal = reciprocal;
}


/**
 * Calculates the contribution of a tie flip to the given actor.
 */
double CovariateSimilarityEffect::calculateContribution(int alter) const
{
	double change = 0;

	if (!this->lreciprocal || this->inTieExists(alter))
	{
		change = this->similarity(this->ego(), alter);
	}

	return change;
}


/**
 * Detailed comment in the base class.
 */
double CovariateSimilarityEffect::statistic(
	const Network * pSummationTieNetwork) const
{
	double statistic = 0;
	const Network * pNetwork = this->pNetwork();
	int n = pNetwork->n();

	for (int i = 0; i < n; i++)
	{
		if (!this->missing(i))
		{
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
						statistic +=
							this->similarity(i, iter.actor());
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
						statistic +=
							this->similarity(i, iter.actor());
					}
				}
			}
		}
	}

	return statistic;
}

}
