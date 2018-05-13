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
		const EffectInfo * pEffectInfo, bool root) :
	NetworkEffect(pEffectInfo)
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
 * The contribution of the tie from the implicit ego to the given alter
 * to the statistic. It is assumed that preprocessEgo(ego) has been
 * called before.
 */
double OutdegreePopularityEffect::tieStatistic(int alter)
{
	double statistic;
	const Network * pNetwork = this->pNetwork();
	int degree = pNetwork->outDegree(alter);

	if (this->lroot)
	{
		statistic = this->lsqrtTable->sqrt(degree);
	}
	else
	{
		// There was a bug here until version 1.1-219
		statistic = degree;
	}

	return statistic;
}

}
