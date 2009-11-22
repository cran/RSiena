/******************************************************************************
 * SIENA: Simulation Investigation for Empirical Network Analysis
 *
 * Web: http://www.stats.ox.ac.uk/~snijders/siena/
 *
 * File: IndegreeActivityEffect.cpp
 *
 * Description: This file contains the implementation of the
 * IndegreeActivityEffect class.
 *****************************************************************************/

#include <cmath>

#include "IndegreeActivityEffect.h"
#include "utils/SqrtTable.h"
#include "network/Network.h"
#include "model/variables/NetworkVariable.h"

namespace siena
{

/**
 * Constructor.
 */
IndegreeActivityEffect::IndegreeActivityEffect(
	const EffectInfo * pEffectInfo, bool root) : NetworkEffect(pEffectInfo)
{
	this->lroot = root;
	this->lsqrtTable = SqrtTable::instance();
}


/**
 * Calculates the contribution of a tie flip to the given actor.
 */
double IndegreeActivityEffect::calculateContribution(int alter) const
{
	double change =
		this->pNetwork()->inDegree(this->ego());

	if (this->lroot)
	{
		change = this->lsqrtTable->sqrt(change);
	}

	return change;
}


/**
 * See base class for a detailed comment.
 */
double IndegreeActivityEffect::statistic(const Network * pSummationTieNetwork)
	const
{
	// A one-mode network is assumed

	double statistic = 0;
	const Network * pNetwork = this->pNetwork();

	// Iterate over all actors and sum up their indegrees x outdegree
	// products (the outdegrees should be taken from the summation network
	// to support endowment effects).

	for (int i = 0; i < pNetwork->n(); i++)
	{
		double indegreeContribution = pNetwork->inDegree(i);

		if (this->lroot)
		{
			indegreeContribution =
				this->lsqrtTable->sqrt(indegreeContribution);
		}

		statistic += indegreeContribution * pSummationTieNetwork->outDegree(i);
	}

	return statistic;
}

}
