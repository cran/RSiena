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
 * Detailed comment in the base class.
 */
double CovariateAlterEffect::statistic(const Network * pSummationTieNetwork)
	const
{
	double statistic = 0;
	const Network * pNetwork = this->pNetwork();

	for (int i = 0; i < pNetwork->m(); i++)
	{
		if (!this->missing(i))
		{
			double value = this->value(i);

			if (this->lsquared)
			{
				value *= value;
			}

			statistic += pSummationTieNetwork->inDegree(i) * value;
		}
	}

	return statistic;
}

}
