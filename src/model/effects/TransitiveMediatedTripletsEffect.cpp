/******************************************************************************
 * SIENA: Simulation Investigation for Empirical Network Analysis
 * 
 * Web: http://www.stats.ox.ac.uk/~snijders/siena/
 * 
 * File: TransitiveMediatedTripletsEffect.cpp
 * 
 * Description: This file contains the implementation of the class
 * TransitiveMediatedTripletsEffect.
 *****************************************************************************/

#include "TransitiveMediatedTripletsEffect.h"
#include "data/OneModeNetwork.h"
#include "data/TieIterator.h"
#include "model/variables/NetworkVariable.h"
#include "model/tables/ConfigurationTable.h"

namespace siena
{

/**
 * Constructor.
 */
TransitiveMediatedTripletsEffect::TransitiveMediatedTripletsEffect(
	const EffectInfo * pEffectInfo) : NetworkEffect(pEffectInfo)
{
}


/**
 * Calculates the contribution of a tie flip to the given actor.
 */
double TransitiveMediatedTripletsEffect::calculateTieFlipContribution(
	int alter) const
{
	// The number of out-stars to the ego and the given alter is the amount
	// of change for this tie flip.
	
	double change = this->pVariable()->pOutStarTable()->get(alter);
	
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
bool TransitiveMediatedTripletsEffect::usesTable(
	const ConfigurationTable * pTable) const
{
	return pTable == this->pVariable()->pOutStarTable();
}


/**
 * Returns the statistic corresponding to this effect as part of
 * the evaluation function with respect to the given network.
 */
double TransitiveMediatedTripletsEffect::evaluationStatistic(
	Network * pNetwork) const
{
	OneModeNetwork * pOneModeNetwork = (OneModeNetwork *) pNetwork;
	int statistic = 0;
	
	for (TieIterator iter = pOneModeNetwork->ties(); iter.valid(); iter.next())
	{
		statistic +=
			pOneModeNetwork->outTwoStarCount(iter.ego(), iter.alter());
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
double TransitiveMediatedTripletsEffect::endowmentStatistic(
	Network * pInitialNetwork,
	Network * pLostTieNetwork) const
{
	OneModeNetwork * pOneModeNetwork = (OneModeNetwork *) pInitialNetwork;
	int statistic = 0;
	
	for (TieIterator iter = pLostTieNetwork->ties(); iter.valid(); iter.next())
	{
		statistic +=
			pOneModeNetwork->outTwoStarCount(iter.ego(), iter.alter());
	}
	
	return statistic;
}

}
