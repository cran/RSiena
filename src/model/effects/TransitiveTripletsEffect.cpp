/******************************************************************************
 * SIENA: Simulation Investigation for Empirical Network Analysis
 * 
 * Web: http://www.stats.ox.ac.uk/~snijders/siena/
 * 
 * File: TransitiveTripletsEffect.cpp
 * 
 * Description: This file contains the implementation of the class
 * TransitiveTripletsEffect.
 *****************************************************************************/

#include "TransitiveTripletsEffect.h"
#include "data/OneModeNetwork.h"
#include "data/TieIterator.h"
#include "model/variables/NetworkVariable.h"
#include "model/tables/ConfigurationTable.h"

namespace siena
{

/**
 * Constructor.
 */
TransitiveTripletsEffect::TransitiveTripletsEffect(
	const EffectInfo * pEffectInfo) : NetworkEffect(pEffectInfo)
{
}


/**
 * Calculates the contribution of a tie flip to the given actor.
 */
double TransitiveTripletsEffect::calculateTieFlipContribution(int alter) const
{
	// If we are introducing a tie from the ego i to the alter j, then each
	// two-path from i to j contributes one transitive triplet, just as each
	// in-star between i and j.
	
	double change =
		this->pVariable()->pTwoPathTable()->get(alter) +
		this->pVariable()->pInStarTable()->get(alter);
	
	// Withdrawals contribute negatively.
	
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
bool TransitiveTripletsEffect::usesTable(const ConfigurationTable * pTable)
	const
{
	return pTable == this->pVariable()->pTwoPathTable() ||
		pTable == this->pVariable()->pInStarTable();
}


/**
 * Returns the statistic corresponding to this effect as part of
 * the evaluation function with respect to the given network.
 */
double TransitiveTripletsEffect::evaluationStatistic(Network * pNetwork) const
{
	OneModeNetwork * pOneModeNetwork = (OneModeNetwork *) pNetwork;
	int transitiveTripletCount = 0;
	
	for (TieIterator iter = pOneModeNetwork->ties(); iter.valid(); iter.next())
	{
		transitiveTripletCount +=
			pOneModeNetwork->twoPathCount(iter.ego(), iter.alter());
	}
	
	return transitiveTripletCount;
}


/**
 * Returns the statistic corresponding to this effect as part of
 * the endowment function with respect to an initial network
 * and a network of lost ties. The current network is implicit as
 * the introduced ties are not relevant for calculating
 * endowment statistics.
 */
double TransitiveTripletsEffect::endowmentStatistic(Network * pInitialNetwork,
	Network * pLostTieNetwork) const
{
	OneModeNetwork * pOneModeNetwork = (OneModeNetwork *) pInitialNetwork;
	int counter = 0;
	
	for (TieIterator iter = pLostTieNetwork->ties(); iter.valid(); iter.next())
	{
		counter +=
			pOneModeNetwork->twoPathCount(iter.ego(), iter.alter());
	}
	
	return counter;
}

}
