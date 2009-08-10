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
#include "data/OneModeNetwork.h"
#include "data/TieIterator.h"
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
double TransitiveTiesEffect::calculateTieFlipContribution(int alter) const
{
	// Suppose we introduce the tie from the ego i to the alter j, which causes
	// another tie (i,h) to become transitive. It means that there was no
	// two-path from i to h before and that i -> j -> h is the only two-path
	// from i to h now. Hence <(i,h),(j,h)> is one of the critical in-stars
	// between i and j.
	
	double change =
		this->pVariable()->pCriticalInStarTable()->get(alter);
	
	// Test if the tie (i,j) is transitive itself (the efficiency could be
	// improved a bit, if we tested only for the existence of a two-path).
	
	if (this->pVariable()->pTwoPathTable()->get(alter) > 0)
	{
		change++;
	}
	
	// If we are actually withdrawing the tie, then we are loosing this
	// many transitive ties.
	
	if (this->pVariable()->outTieExists(alter))
	{
		change = -change;
	}
	
	return change;
}


/**
 * Returns if the given configuration table is used by this effect
 * during the calculation of tie flip contributions.
 */
bool TransitiveTiesEffect::usesTable(const ConfigurationTable * pTable) const
{
	return pTable == this->pVariable()->pCriticalInStarTable() ||
		pTable == this->pVariable()->pTwoPathTable();
}


/**
 * Returns the statistic corresponding to this effect as part of
 * the evaluation function with respect to the given network.
 */
double TransitiveTiesEffect::evaluationStatistic(Network * pNetwork) const
{
	OneModeNetwork * pOneModeNetwork = (OneModeNetwork *) pNetwork;
	int statistic = 0;
	
	for (TieIterator iter = pOneModeNetwork->ties(); iter.valid(); iter.next())
	{
		if (pOneModeNetwork->twoPathCount(iter.ego(), iter.alter()) > 0)
		{
			statistic++;
		}
	}
	
	// TODO: Shouldn't we divide by 2 for symmetric networks?
	return statistic;
}


/**
 * Returns the statistic corresponding to this effect as part of
 * the endowment function with respect to an initial network
 * and a network of lost ties. The current network is implicit as
 * the introduced ties are not relevant for calculating
 * endowment statistics.
 */
double TransitiveTiesEffect::endowmentStatistic(Network * pInitialNetwork,
	Network * pLostTieNetwork) const
{
	OneModeNetwork * pOneModeNetwork = (OneModeNetwork *) pInitialNetwork;
	int statistic = 0;
	
	for (TieIterator iter = pLostTieNetwork->ties(); iter.valid(); iter.next())
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
