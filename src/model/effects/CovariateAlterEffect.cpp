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
#include "data/Network.h"
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
double CovariateAlterEffect::calculateTieFlipContribution(int alter) const
{
	double change = this->value(alter);

	if (this->lsquared)
	{
		change *= change;
	}

	if (this->pVariable()->outTieExists(alter))
	{
		// The ego would loose the tie and the covariate value of the alter
		change = -change;
	}

	return change;
}


/**
 * Returns the statistic corresponding to this effect as part of
 * the evaluation function with respect to the given network.
 */
double CovariateAlterEffect::evaluationStatistic(Network * pNetwork) const
{
	double statistic = 0;

	for (int i = 0; i < pNetwork->m(); i++)
	{
		if (!this->missing(i))
		{
			double value = this->value(i);

			if (this->lsquared)
			{
				value *= value;
			}

			statistic += pNetwork->inDegree(i) * value;
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
double CovariateAlterEffect::endowmentStatistic(Network * pInitialNetwork,
	Network * pLostTieNetwork) const
{
	// This is the same as the evaluation statistic computed with respect
	// to the network of lost ties.

	return this->evaluationStatistic(pLostTieNetwork);
}

}
