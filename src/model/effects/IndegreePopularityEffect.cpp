/******************************************************************************
 * SIENA: Simulation Investigation for Empirical Network Analysis
 *
 * Web: http://www.stats.ox.ac.uk/~snijders/siena/
 *
 * File: IndegreePopularityEffect.cpp
 *
 * Description: This file contains the implementation of the
 * IndegreePopularityEffect class.
 *****************************************************************************/

#include "IndegreePopularityEffect.h"
#include "utils/SqrtTable.h"
#include "network/Network.h"
#include "model/variables/NetworkVariable.h"

namespace siena
{

/**
 * Constructor.
 */
IndegreePopularityEffect::IndegreePopularityEffect(
	const EffectInfo * pEffectInfo, bool root) : NetworkEffect(pEffectInfo)
{
	this->lroot = root;
	this->lsqrtTable = SqrtTable::instance();
}


/**
 * Calculates the contribution of a tie flip to the given actor.
 */
double IndegreePopularityEffect::calculateContribution(int alter) const
{
	int degree = this->pNetwork()->inDegree(alter);

	if (!this->outTieExists(alter))
	{
		// The indegree will increase after introducing the tie
		degree++;
	}

	double change = degree;

	if (this->lroot)
	{
		change = this->lsqrtTable->sqrt(degree);
	}

	return change;
}


/**
 * See base class for a detailed comment.
 */
double IndegreePopularityEffect::statistic(
	const Network * pSummationTieNetwork) const
{
	double statistic = 0;
	const Network * pNetwork = this->pNetwork();

	// Iterate over receivers and sum up their indegrees (or roots of
	// indegrees) multiplied by the number of incoming summation ties.

	for (int i = 0; i < pNetwork->m(); i++)
	{
		int degree = pNetwork->inDegree(i);

		if (this->lroot)
		{
			statistic +=
				pSummationTieNetwork->inDegree(i) *
					this->lsqrtTable->sqrt(degree);
		}
		else
		{
			statistic += pSummationTieNetwork->inDegree(i) * degree;
		}
	}

	if (!this->lroot)
	{
		// TODO: Why do we multiply by the number of senders? Ask Tom.
		statistic *= pNetwork->n();
	}

	return statistic;
}

}
