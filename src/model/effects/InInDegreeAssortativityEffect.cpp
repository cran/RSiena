/******************************************************************************
 * SIENA: Simulation Investigation for Empirical Network Analysis
 *
 * Web: http://www.stats.ox.ac.uk/~snijders/siena/
 *
 * File: InInDegreeAssortativityEffect.cpp
 *
 * Description: This file contains the implementation of the
 * InInDegreeAssortativityEffect class.
 *****************************************************************************/

#include <cmath>

#include "InInDegreeAssortativityEffect.h"
#include "utils/SqrtTable.h"
#include "data/Network.h"
#include "data/TieIterator.h"
#include "model/EffectInfo.h"
#include "model/variables/NetworkVariable.h"

namespace siena
{

/**
 * Constructor.
 */
InInDegreeAssortativityEffect::InInDegreeAssortativityEffect(
	const EffectInfo * pEffectInfo) : NetworkEffect(pEffectInfo)
{
	this->lroot = pEffectInfo->internalEffectParameter() == 2;
	this->lsqrtTable = SqrtTable::instance();
}


/**
 * Calculates the contribution of a tie flip to the given actor.
 */
double InInDegreeAssortativityEffect::calculateTieFlipContribution(int alter)
	const
{
	double change = 0;
	Network * pNetwork = this->pVariable()->pNetwork();
	int egoDegree = pNetwork->inDegree(this->pVariable()->ego());
	int alterDegree = pNetwork->inDegree(alter);

	if (this->pVariable()->outTieExists(alter))
	{
		if (this->lroot)
		{
			change =
				- this->lsqrtTable->sqrt(egoDegree) *
					this->lsqrtTable->sqrt(alterDegree);
		}
		else
		{
			change = - egoDegree * alterDegree;
		}
	}
	else
	{
		if (this->lroot)
		{
			change =
				this->lsqrtTable->sqrt(egoDegree) *
					this->lsqrtTable->sqrt(alterDegree + 1);
		}
		else
		{
			change = egoDegree * (alterDegree + 1);
		}
	}

	return change;
}


/**
 * Returns the statistic corresponding to this effect as part of
 * the evaluation function with respect to the given network.
 */
double InInDegreeAssortativityEffect::evaluationStatistic(Network * pNetwork)
	const
{
	double statistic = 0;

	for (TieIterator iter = pNetwork->ties(); iter.valid(); iter.next())
	{
		int egoDegree = pNetwork->inDegree(iter.ego());
		int alterDegree = pNetwork->inDegree(iter.alter());

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


/**
 * Returns the statistic corresponding to this effect as part of
 * the endowment function with respect to an initial network
 * and a network of lost ties. The current network is implicit as
 * the introduced ties are not relevant for calculating
 * endowment statistics.
 */
double InInDegreeAssortativityEffect::endowmentStatistic(
	Network * pInitialNetwork,
	Network * pLostTieNetwork) const
{
	double statistic = 0;

	for (TieIterator iter = pLostTieNetwork->ties(); iter.valid(); iter.next())
	{
		int egoDegree = pInitialNetwork->inDegree(iter.ego());
		int alterDegree = pInitialNetwork->inDegree(iter.alter());

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
