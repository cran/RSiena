/******************************************************************************
 * SIENA: Simulation Investigation for Empirical Network Analysis
 *
 * Web: http://www.stats.ox.ac.uk/~snijders/siena/
 *
 * File: OutdegreeActivitySqrtEffect.cpp
 *
 * Description: This file contains the implementation of the
 * OutdegreeActivitySqrtEffect class.
 *****************************************************************************/

#include <cmath>

#include "OutdegreeActivitySqrtEffect.h"
#include "utils/SqrtTable.h"
#include "network/Network.h"
#include "model/variables/NetworkVariable.h"

namespace siena
{

/**
 * Constructor.
 */
OutdegreeActivitySqrtEffect::OutdegreeActivitySqrtEffect(
	const EffectInfo * pEffectInfo) : NetworkEffect(pEffectInfo)
{
	this->lsqrtTable = SqrtTable::instance();
}


/**
 * Calculates the contribution of a tie flip to the given actor.
 */
double OutdegreeActivitySqrtEffect::calculateContribution(int alter)
	const
{
	// Current out-degree

	int degree = this->pNetwork()->outDegree(this->ego());

	if (this->outTieExists(alter))
	{
		// The tie is going to be withdrawn, so the degree decreases

		return degree * this->lsqrtTable->sqrt(degree) -
			(degree - 1) * this->lsqrtTable->sqrt(degree - 1);
	}
	else
	{
		// The tie is going to be introduced, so the degree increases

		return (degree + 1) * this->lsqrtTable->sqrt(degree + 1) -
			degree * this->lsqrtTable->sqrt(degree);
	}
}


/**
 * See base class.
 */
double OutdegreeActivitySqrtEffect::statistic(
	const Network * pSummationTieNetwork) const
{
	double statistic = 0;
	const Network * pNetwork = this->pNetwork();

	// Iterate over senders and sum up the square roots of their
	// outdegrees multiplied by the number of outgoing ties in the summation
	// network.

	for (int i = 0; i < pNetwork->n(); i++)
	{
		statistic +=
			this->lsqrtTable->sqrt(pNetwork->outDegree(i)) *
				pSummationTieNetwork->outDegree(i);
	}

	return statistic;
}

}
