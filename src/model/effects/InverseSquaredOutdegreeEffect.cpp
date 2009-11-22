/******************************************************************************
 * SIENA: Simulation Investigation for Empirical Network Analysis
 *
 * Web: http://www.stats.ox.ac.uk/~snijders/siena/
 *
 * File: InverseSquaredOutdegreeEffect.cpp
 *
 * Description: This file contains the implementation of the class
 * InverseSquaredOutdegreeEffect.
 *****************************************************************************/

#include <stdexcept>
#include "InverseSquaredOutdegreeEffect.h"
#include "network/Network.h"
#include "model/EffectInfo.h"
#include "model/variables/NetworkVariable.h"

namespace siena
{

/**
 * Constructor.
 */
InverseSquaredOutdegreeEffect::InverseSquaredOutdegreeEffect(
	const EffectInfo * pEffectInfo) : NetworkEffect(pEffectInfo)
{
	this->lc = pEffectInfo->internalEffectParameter();

	if (this->lc < 1)
	{
		throw invalid_argument(
			string("InverseSquaredOutdegreeEffect: ") +
			"Parameter value must be at least 1");
	}
}


/**
 * Calculates the contribution of a tie flip to the given actor.
 */
double InverseSquaredOutdegreeEffect::calculateContribution(int alter)
	const
{
	double sum = this->pNetwork()->outDegree(this->ego()) + this->lc;

	if (this->outTieExists(alter))
	{
		// Tie withdrawal
		return -2.0 / ((sum - 1) * sum * (sum + 1));
	}
	else
	{
		// Tie introduction
		return -2.0 / (sum * (sum + 1) * (sum + 2));
	}
}


/**
 * Returns the statistic corresponding to this effect as part of
 * the evaluation function.
 */
double InverseSquaredOutdegreeEffect::evaluationStatistic() const
{
	double statistic = 0;
	const Network * pNetwork = this->pNetwork();
	int n = pNetwork->n();

	for (int i = 0; i < n; i++)
	{
		double sum = pNetwork->outDegree(i) + this->lc;
		statistic += 1.0 / (sum * (sum + 1));
	}

	return statistic;
}


/**
 * Returns the statistic corresponding to this effect as part of
 * the endowment function.
 */
double InverseSquaredOutdegreeEffect::endowmentStatistic(
	Network * pLostTieNetwork) const
{
	throw logic_error(
		"InverseSquaredOutdegreeEffect: Endowment effect not supported.");
}

}
