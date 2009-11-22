/******************************************************************************
 * SIENA: Simulation Investigation for Empirical Network Analysis
 *
 * Web: http://www.stats.ox.ac.uk/~snijders/siena/
 *
 * File: TransitiveTiesEffect.cpp
 *
 * Description: This file contains the implementation of the class
 * TransitiveTiesEffect.
 *****************************************************************************/

#include "TransitiveTiesEffect.h"
#include "network/OneModeNetwork.h"
#include "network/TieIterator.h"
#include "model/variables/NetworkVariable.h"
#include "model/tables/ConfigurationTable.h"

namespace siena
{

/**
 * Constructor.
 */
TransitiveTiesEffect::TransitiveTiesEffect(const EffectInfo * pEffectInfo) :
	NetworkEffect(pEffectInfo)
{
}


/**
 * Calculates the contribution of a tie flip to the given actor.
 */
double TransitiveTiesEffect::calculateContribution(int alter) const
{
	// Suppose we introduce the tie from the ego i to the alter j, which causes
	// another tie (i,h) to become transitive. It means that there was no
	// two-path from i to h before and that i -> j -> h is the only two-path
	// from i to h now. Hence <(i,h),(j,h)> is one of the critical in-stars
	// between i and j.

	double change =
		this->pCriticalInStarTable()->get(alter);

	// Test if the tie (i,j) is transitive itself (the efficiency could be
	// improved a bit, if we tested only for the existence of a two-path).

	if (this->pTwoPathTable()->get(alter) > 0)
	{
		change++;
	}

	return change;
}


/**
 * See base class.
 */
double TransitiveTiesEffect::statistic(const Network * pSummationTieNetwork)
	const
{
	OneModeNetwork * pOneModeNetwork = (OneModeNetwork *) this->pNetwork();
	int statistic = 0;

	for (TieIterator iter = pSummationTieNetwork->ties();
		iter.valid();
		iter.next())
	{
		if (pOneModeNetwork->twoPathCount(iter.ego(), iter.alter()) > 0)
		{
			statistic++;
		}
	}

	// TODO: Shouldn't we divide by 2 for symmetric networks?
	return statistic;
}

}
