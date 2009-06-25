/******************************************************************************
 * SIENA: Simulation Investigation for Empirical Network Analysis
 * 
 * Web: http://www.stats.ox.ac.uk/~snijders/siena/
 * 
 * File: TransitiveTriadsEffect.cpp
 * 
 * Description: This file contains the implementation of the class
 * TransitiveTriadsEffect.
 *****************************************************************************/

#include "TransitiveTriadsEffect.h"
#include "data/OneModeNetwork.h"
#include "data/TieIterator.h"
#include "model/variables/NetworkVariable.h"
#include "model/tables/ConfigurationTable.h"

namespace siena
{

/**
 * Constructor.
 */
TransitiveTriadsEffect::TransitiveTriadsEffect(
	const EffectInfo * pEffectInfo) : NetworkEffect(pEffectInfo)
{
}


/**
 * Calculates the contribution of a tie flip to the given actor.
 */
double TransitiveTriadsEffect::calculateTieFlipContribution(int alter) const
{
	double change = this->pVariable()->pTwoPathTable()->get(alter);
	
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
bool TransitiveTriadsEffect::usesTable(const ConfigurationTable * pTable) const
{
	return pTable == this->pVariable()->pTwoPathTable();
}


/**
 * Returns the statistic corresponding to this effect as part of
 * the evaluation function with respect to the given network.
 */
double TransitiveTriadsEffect::evaluationStatistic(Network * pNetwork) const
{
	OneModeNetwork * pOneModeNetwork = (OneModeNetwork *) pNetwork;
	int transitiveTriadCount = 0;
	
	for (TieIterator iter = pOneModeNetwork->ties(); iter.valid(); iter.next())
	{
		transitiveTriadCount +=
			pOneModeNetwork->twoPathCount(iter.ego(), iter.alter());
	}
	
	return transitiveTriadCount / 6.0;
}


/**
 * Returns the statistic corresponding to this effect as part of
 * the endowment function with respect to an initial network
 * and a network of lost ties. The current network is implicit as
 * the introduced ties are not relevant for calculating
 * endowment statistics.
 */
double TransitiveTriadsEffect::endowmentStatistic(Network * pInitialNetwork,
	Network * pLostTieNetwork) const
{
	OneModeNetwork * pOneModeNetwork = (OneModeNetwork *) pInitialNetwork;
	int counter = 0;
	
	for (TieIterator iter = pLostTieNetwork->ties(); iter.valid(); iter.next())
	{
		counter +=
			pOneModeNetwork->twoPathCount(iter.ego(), iter.alter());
	}
	
	return counter / 6.0;
}

}
