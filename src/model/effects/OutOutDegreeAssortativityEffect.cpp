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
#include "network/Network.h"
#include "network/TieIterator.h"
#include "network/IncidentTieIterator.h"
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
void OutOutDegreeAssortativityEffect::preprocessEgo(int ego)
{
	NetworkEffect::preprocessEgo(ego);

	const Network * pNetwork = this->pNetwork();
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
double OutOutDegreeAssortativityEffect::calculateContribution(int alter)
	const
{
	double change = 0;
	int alterDegree = this->pNetwork()->outDegree(alter);

	if (this->outTieExists(alter))
	{
		if (this->lroot)
		{
			double sqrtAlterDegree = this->lsqrtTable->sqrt(alterDegree);

			change =
				(this->lneighborDegreeSum - sqrtAlterDegree) *
					(this->lsqrtDegree - this->lsqrtDegreeMinus) +
				this->lsqrtDegree * sqrtAlterDegree;
		}
		else
		{
			change =
				(this->lneighborDegreeSum - alterDegree) +
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
 * See base class.
 */
double OutOutDegreeAssortativityEffect::statistic(
	const Network * pSummationTieNetwork) const
{
	double statistic = 0;
	const Network * pNetwork = this->pNetwork();

	for (TieIterator iter = pSummationTieNetwork->ties();
		iter.valid();
		iter.next())
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

}
