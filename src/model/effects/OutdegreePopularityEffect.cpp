/******************************************************************************
 * SIENA: Simulation Investigation for Empirical Network Analysis
 *
 * Web: http://www.stats.ox.ac.uk/~snijders/siena/
 *
 * File: OutdegreePopularityEffect.cpp
 *
 * Description: This file contains the implementation of the
 * OutdegreePopularityEffect class.
 *****************************************************************************/

#include <cmath>

#include "OutdegreePopularityEffect.h"
#include "utils/SqrtTable.h"
#include "network/Network.h"
#include "model/variables/NetworkVariable.h"

namespace siena
{

/**
 * Constructor.
 */
OutdegreePopularityEffect::OutdegreePopularityEffect(
	const EffectInfo * pEffectInfo,
	bool root) : NetworkEffect(pEffectInfo)
{
	this->lroot = root;
	this->lsqrtTable = SqrtTable::instance();
}


/**
 * Calculates the contribution of a tie flip to the given actor.
 */
double OutdegreePopularityEffect::calculateContribution(int alter) const
{
	int degree = this->pNetwork()->outDegree(alter);
	double change = degree;

	if (this->lroot)
	{
		change = this->lsqrtTable->sqrt(degree);
	}

	return change;
}


/**
 * See base class.
 */
double OutdegreePopularityEffect::statistic(
	const Network * pSummationTieNetwork) const
{
	// We assume a one-mode network because this effect doesn't make sense
	// for two-mode networks.

	double statistic = 0;
	const Network * pNetwork = this->pNetwork();

	// Iterate over all actors and sum up their outdegrees
	// (or square roots of outdegrees) multiplied by
	// the number of incoming ties in the summation network.

	for (int i = 0; i < pNetwork->n(); i++)
	{
		int degree = pNetwork->outDegree(i);
		double outDegreeContribution = degree;

		if (this->lroot)
		{
			outDegreeContribution = this->lsqrtTable->sqrt(degree);
		}

		statistic += outDegreeContribution * pSummationTieNetwork->inDegree(i);
	}

	if (!this->lroot)
	{
		// TODO: Why do we multiply by the number of senders? Ask Tom.
		statistic *= pNetwork->n();
	}

	return statistic;
}

}
