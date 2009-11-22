/******************************************************************************
 * SIENA: Simulation Investigation for Empirical Network Analysis
 *
 * Web: http://www.stats.ox.ac.uk/~snijders/siena/
 *
 * File: InOutDegreeAssortativityEffect.cpp
 *
 * Description: This file contains the implementation of the
 * InOutDegreeAssortativityEffect class.
 *****************************************************************************/

#include <cmath>

#include "InOutDegreeAssortativityEffect.h"
#include "utils/SqrtTable.h"
#include "network/Network.h"
#include "network/TieIterator.h"
#include "model/EffectInfo.h"
#include "model/variables/NetworkVariable.h"

namespace siena
{

/**
 * Constructor.
 */
InOutDegreeAssortativityEffect::InOutDegreeAssortativityEffect(
	const EffectInfo * pEffectInfo) : NetworkEffect(pEffectInfo)
{
	this->lroot = pEffectInfo->internalEffectParameter() == 2;
	this->lsqrtTable = SqrtTable::instance();
}


/**
 * Calculates the contribution of a tie flip to the given actor.
 */
double InOutDegreeAssortativityEffect::calculateContribution(int alter)
	const
{
	double change = 0;
	const Network * pNetwork = this->pNetwork();
	int egoDegree = pNetwork->inDegree(this->ego());
	int alterDegree = pNetwork->outDegree(alter);

	if (this->lroot)
	{
		change =
			this->lsqrtTable->sqrt(egoDegree) *
				this->lsqrtTable->sqrt(alterDegree);
	}
	else
	{
		change = egoDegree * alterDegree;
	}

	return change;
}


/**
 * See base class.
 */
double InOutDegreeAssortativityEffect::statistic(
	const Network * pSummationTieNetwork) const
{
	double statistic = 0;
	const Network * pNetwork = this->pNetwork();

	for (TieIterator iter = pSummationTieNetwork->ties();
		iter.valid();
		iter.next())
	{
		int egoDegree = pNetwork->inDegree(iter.ego());
		int alterDegree = pNetwork->outDegree(iter.alter());

		if (this->lroot)
		{
			statistic +=
				this->lsqrtTable->sqrt(egoDegree) *
					this->lsqrtTable->sqrt(alterDegree);
		}
		else
		{
			statistic += egoDegree * alterDegree;
		}
	}

	return statistic;
}

}
