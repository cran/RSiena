/******************************************************************************
 * SIENA: Simulation Investigation for Empirical Network Analysis
 *
 * Web: http://www.stats.ox.ac.uk/~snijders/siena/
 *
 * File: CovariateAlterEffect.cpp
 *
 * Description: This file contains the implementation of the
 * CovariateAlterEffect class.
 *****************************************************************************/

#include "CovariateAlterEffect.h"
#include "network/Network.h"
#include "model/variables/NetworkVariable.h"

namespace siena
{

/**
 * Constructor.
 * @param[in] pEffectInfo the effect descriptor
 * @param[in] squared indicates if the covariate values must be squared
 */
CovariateAlterEffect::CovariateAlterEffect(const EffectInfo * pEffectInfo,
	bool squared) :
		CovariateDependentNetworkEffect(pEffectInfo)
{
	this->lsquared = squared;
}


/**
 * Calculates the contribution of a tie flip to the given actor.
 */
double CovariateAlterEffect::calculateContribution(int alter) const
{
	double change = this->value(alter);

	if (this->lsquared)
	{
		change *= change;
	}

	return change;
}


/**
 * The contribution of the tie from the implicit ego to the given alter
 * to the statistic. It is assumed that preprocessEgo(ego) has been
 * called before.
 */
double CovariateAlterEffect::tieStatistic(int alter)
{
	double statistic = 0;

	if (!this->missing(alter))
	{
		statistic = this->value(alter);

		if (this->lsquared)
		{
			statistic *= statistic;
		}
	}

	return statistic;
}

}
