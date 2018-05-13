/******************************************************************************
 * SIENA: Simulation Investigation for Empirical Network Analysis
 *
 * Web: http://www.stats.ox.ac.uk/~snijders/siena/
 *
 * File: SameCovariateActivityEffect.cpp
 *
 * Description: This file contains the implementation of the
 * SameCovariateActivityEffect class.
 *****************************************************************************/

#include <cmath>
#include "SameCovariateActivityEffect.h"
#include "utils/Utils.h"
#include "network/Network.h"
#include "model/variables/NetworkVariable.h"
#include "network/IncidentTieIterator.h"

namespace siena
{

/**
 * Constructor.
 * @param[in] pEffectInfo the effect descriptor
 * @param[in] squared indicates if the covariate values must be squared
 */
SameCovariateActivityEffect::SameCovariateActivityEffect(
		const EffectInfo * pEffectInfo, bool same) :
	CovariateDependentNetworkEffect(pEffectInfo)
{
	this->lsame = same;
}


/**
 * Calculates the contribution of a tie flip to the given actor.
 */
double SameCovariateActivityEffect::calculateContribution(int alter) const
{
	double myvalue = this->value(this->ego());
	double contribution = 0;
	const Network * pNetwork = this->pNetwork();

	if ((lsame) && (fabs(this->value(alter) - myvalue) < EPSILON))
	{
		for (IncidentTieIterator iter = pNetwork->outTies(this->ego());
			iter.valid();
			iter.next())
		{
			// Get the receiver of the outgoing tie.
			int h = iter.actor();
			if (fabs(this->value(h) - myvalue) < EPSILON)
			{
				contribution++;
			}
		}
		if (this->outTieExists(alter))
		{
			contribution--;
		}
		contribution *= 2;
		contribution++;
	}

	if ((!lsame) && (fabs(this->value(alter) - myvalue) >= EPSILON))
	{
		for (IncidentTieIterator iter = pNetwork->outTies(this->ego());
			iter.valid();
			iter.next())
		{
			// Get the receiver of the outgoing tie.
			int h = iter.actor();
			if (fabs(this->value(h) - myvalue) >= EPSILON)
			{
				contribution++;
			}
		}
		if (this->outTieExists(alter))
		{
			contribution--;
		}
		contribution *= 2;
		contribution++;
	}

	return contribution;
}


/**
 * The contribution of the tie from the implicit ego to the given alter
 * to the statistic. It is assumed that preprocessEgo(ego) has been
 * called before.
 */
double SameCovariateActivityEffect::tieStatistic(int alter)
{
	double contribution = 0;
	const Network * pNetwork = this->pNetwork();

	if (!((this->missing(alter)) || (this->missing(this->ego()))))
	{
		double myvalue = this->value(this->ego());

		if (lsame)
		{
			if (fabs(this->value(alter) - myvalue) < EPSILON)
			{
				for (IncidentTieIterator iter = pNetwork->outTies(this->ego());
					iter.valid();
					iter.next())
				{
					// Get the receiver of the outgoing tie.
					int h = iter.actor();
					if ((!this->missing(h)) &&
							(fabs(this->value(h) - myvalue) < EPSILON))
					{
						contribution++;
					}
				}
			}
		}
		else
		{
			if (fabs(this->value(alter) - myvalue) >= EPSILON)
			{
				for (IncidentTieIterator iter = pNetwork->outTies(this->ego());
					iter.valid();
					iter.next())
				{
					// Get the receiver of the outgoing tie.
					int h = iter.actor();
					if ((!this->missing(h)) &&
							(fabs(this->value(h) - myvalue) >= EPSILON))
					{
						contribution++;
					}
				}
			}
		}
	}

	return contribution;
}

}
