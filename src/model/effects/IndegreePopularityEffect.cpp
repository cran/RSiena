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
 * The contribution of the tie from the implicit ego to the given alter
 * to the statistic. It is assumed that preprocessEgo(ego) has been
 * called before.
 */
double IndegreePopularityEffect::tieStatistic(int alter)
{
	const Network * pNetwork = this->pNetwork();
	int degree = pNetwork->inDegree(alter);
	double statistic;

	if (this->lroot)
	{
		statistic = this->lsqrtTable->sqrt(degree);
	}
	else
	{
		// TODO: Why do we multiply by the number of senders? Ask Tom.
		statistic = degree * pNetwork->n();
	}

	return statistic;
}

}
