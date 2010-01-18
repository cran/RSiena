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
 * The contribution of the tie from the implicit ego to the given alter
 * to the statistic. It is assumed that preprocessEgo(ego) has been
 * called before.
 */
double IndegreeActivityEffect::tieStatistic(int alter)
{
	const Network * pNetwork = this->pNetwork();
	int inDegree = pNetwork->inDegree(this->ego());
	double statistic = inDegree;

	if (this->lroot)
	{
		statistic = this->lsqrtTable->sqrt(inDegree);
	}

	return statistic;
}

}
