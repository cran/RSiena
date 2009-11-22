/******************************************************************************
 * SIENA: Simulation Investigation for Empirical Network Analysis
 *
 * Web: http://www.stats.ox.ac.uk/~snijders/siena/
 *
 * File: InverseOutdegreeEffect.cpp
 *
 * Description: This file contains the implementation of the class
 * InverseOutdegreeEffect.
 *****************************************************************************/

#include <stdexcept>
#include "InverseOutdegreeEffect.h"
#include "network/Network.h"
#include "model/EffectInfo.h"
#include "model/variables/NetworkVariable.h"

namespace siena
{

/**
 * Constructor.
 */
InverseOutdegreeEffect::InverseOutdegreeEffect(
	const EffectInfo * pEffectInfo) : NetworkEffect(pEffectInfo)
{
	this->lc = pEffectInfo->internalEffectParameter();

	if (this->lc < 1)
	{
		throw invalid_argument(
			"InverseOutdegreeEffect: Parameter value must be at least 1");
	}
}


/**
 * Calculates the contribution of a tie flip to the given actor.
 */
double InverseOutdegreeEffect::calculateContribution(int alter) const
{
	double sum =
		this->pNetwork()->outDegree(this->ego()) + this->lc;

	if (this->outTieExists(alter))
	{
		// Tie withdrawal
		return -1.0 / ((sum - 1) * sum);
	}
	else
	{
		// Tie introduction
		return -1.0 / ((sum + 1) * sum);
	}
}


/**
 * Returns the statistic corresponding to this effect as part of
 * the evaluation function.
 */
double InverseOutdegreeEffect::evaluationStatistic() const
{
	double statistic = 0;
	const Network * pNetwork = this->pNetwork();
	int n = pNetwork->n();

	for (int i = 0; i < n; i++)
	{
		statistic += 1.0 / (pNetwork->outDegree(i) + this->lc);
	}

	return statistic;
}


/**
 * Returns the statistic corresponding to this effect as part of
 * the endowment function.
 */
double InverseOutdegreeEffect::endowmentStatistic(Network * pLostTieNetwork)
	const
{
	throw logic_error(
		"InverseOutdegreeEffect: Endowment effect not supported.");
}

}
