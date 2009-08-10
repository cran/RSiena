/******************************************************************************
 * SIENA: Simulation Investigation for Empirical Network Analysis
 *
 * Web: http://www.stats.ox.ac.uk/~snijders/siena/
 *
 * File: OutOutDegreeAssortativityEffect.cpp
 *
 * Description: This file contains the implementation of the
 * OutOutDegreeAssortativityEffect class.
 *****************************************************************************/

#include <cmath>
#include <R.h>

#include "OutOutDegreeAssortativityEffect.h"
#include "utils/SqrtTable.h"
#include "data/Network.h"
#include "data/TieIterator.h"
#include "data/IncidentTieIterator.h"
#include "model/EffectInfo.h"
#include "model/variables/NetworkVariable.h"

namespace siena
{

/**
 * Constructor.
 */
OutOutDegreeAssortativityEffect::OutOutDegreeAssortativityEffect(
	const EffectInfo * pEffectInfo) : NetworkEffect(pEffectInfo)
{
	this->lroot = pEffectInfo->internalEffectParameter() == 2;
	this->lsqrtTable = SqrtTable::instance();
}


/**
 * Does the necessary preprocessing work for calculating the tie flip
 * contributions for a specific ego.
 */
void OutOutDegreeAssortativityEffect::preprocessEgo()
{
	int ego = this->pVariable()->ego();
	Network * pNetwork = this->pVariable()->pNetwork();
	this->ldegree = pNetwork->outDegree(ego);

	if (this->lroot)
	{
		this->lsqrtDegree = this->lsqrtTable->sqrt(this->ldegree);
		this->lsqrtDegreePlus = this->lsqrtTable->sqrt(this->ldegree + 1);

		if (this->ldegree >= 1)
		{
			this->lsqrtDegreeMinus = this->lsqrtTable->sqrt(this->ldegree - 1);
		}
	}

	this->lneighborDegreeSum = 0;

	for (IncidentTieIterator iter = pNetwork->outTies(ego);
		iter.valid();
		iter.next())
	{
		int alterDegree = pNetwork->outDegree(iter.actor());

		if (this->lroot)
		{
			this->lneighborDegreeSum += this->lsqrtTable->sqrt(alterDegree);
		}
		else
		{
			this->lneighborDegreeSum += alterDegree;
		}
	}
}


/**
 * Calculates the contribution of a tie flip to the given actor.
 */
double OutOutDegreeAssortativityEffect::calculateTieFlipContribution(int alter)
	const
{
	double change = 0;
	int alterDegree = this->pVariable()->pNetwork()->outDegree(alter);

	if (this->pVariable()->outTieExists(alter))
	{
		if (this->lroot)
		{
			double sqrtAlterDegree = this->lsqrtTable->sqrt(alterDegree);

			change =
				(this->lneighborDegreeSum - sqrtAlterDegree) *
					(this->lsqrtDegreeMinus - this->lsqrtDegree) -
				this->lsqrtDegree * sqrtAlterDegree;
		}
		else
		{
			change =
				- (this->lneighborDegreeSum - alterDegree) -
				this->ldegree * alterDegree;
		}
	}
	else
	{
		if (this->lroot)
		{
			change =
				this->lneighborDegreeSum *
					(this->lsqrtDegreePlus - this->lsqrtDegree) +
				this->lsqrtDegreePlus * this->lsqrtTable->sqrt(alterDegree);
		}
		else
		{
			change =
				this->lneighborDegreeSum + (this->ldegree + 1) * alterDegree;
		}
	}

	return change;
}


/**
 * Returns the statistic corresponding to this effect as part of
 * the evaluation function with respect to the given network.
 */
double OutOutDegreeAssortativityEffect::evaluationStatistic(Network * pNetwork)
	const
{
	double statistic = 0;

	for (TieIterator iter = pNetwork->ties(); iter.valid(); iter.next())
	{
		int egoDegree = pNetwork->outDegree(iter.ego());
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


/**
 * Returns the statistic corresponding to this effect as part of
 * the endowment function with respect to an initial network
 * and a network of lost ties. The current network is implicit as
 * the introduced ties are not relevant for calculating
 * endowment statistics.
 */
double OutOutDegreeAssortativityEffect::endowmentStatistic(
	Network * pInitialNetwork,
	Network * pLostTieNetwork) const
{
	double statistic = 0;

	for (TieIterator iter = pLostTieNetwork->ties(); iter.valid(); iter.next())
	{
		int egoDegree = pInitialNetwork->outDegree(iter.ego());
		int alterDegree = pInitialNetwork->outDegree(iter.alter());

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
