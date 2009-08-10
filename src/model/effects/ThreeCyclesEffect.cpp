/******************************************************************************
 * SIENA: Simulation Investigation for Empirical Network Analysis
 * 
 * Web: http://www.stats.ox.ac.uk/~snijders/siena/
 * 
 * File: ThreeCyclesEffect.cpp
 * 
 * Description: This file contains the implementation of the class
 * ThreeCyclesEffect.
 *****************************************************************************/

#include "ThreeCyclesEffect.h"
#include "data/OneModeNetwork.h"
#include "data/TieIterator.h"
#include "model/variables/NetworkVariable.h"
#include "model/tables/ConfigurationTable.h"

namespace siena
{

/**
 * Constructor.
 */
ThreeCyclesEffect::ThreeCyclesEffect(const EffectInfo * pEffectInfo) :
	NetworkEffect(pEffectInfo)
{
}


/**
 * Calculates the contribution of a tie flip to the given actor.
 */
double ThreeCyclesEffect::calculateTieFlipContribution(int alter) const
{
	// The absolute value of the contribution is going to be the number
	// of two-paths from the alter to the ego -- a number, which is stored
	// in the table of reverse two paths.
	
	double change =
		this->pVariable()->pReverseTwoPathTable()->get(alter);
	
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
bool ThreeCyclesEffect::usesTable(const ConfigurationTable * pTable) const
{
	return pTable == this->pVariable()->pReverseTwoPathTable();
}


/**
 * Returns the statistic corresponding to this effect as part of
 * the evaluation function with respect to the given network.
 */
double ThreeCyclesEffect::evaluationStatistic(Network * pNetwork) const
{
	OneModeNetwork * pOneModeNetwork = (OneModeNetwork *) pNetwork;
	int statistic = 0;
	
	for (TieIterator iter = pOneModeNetwork->ties(); iter.valid(); iter.next())
	{
		statistic +=
			pOneModeNetwork->twoPathCount(iter.alter(), iter.ego());
	}
	
	return statistic / 3.0;
}


/**
 * Returns the statistic corresponding to this effect as part of
 * the endowment function with respect to an initial network
 * and a network of lost ties. The current network is implicit as
 * the introduced ties are not relevant for calculating
 * endowment statistics.
 */
double ThreeCyclesEffect::endowmentStatistic(Network * pInitialNetwork,
	Network * pLostTieNetwork) const
{
	OneModeNetwork * pOneModeNetwork = (OneModeNetwork *) pInitialNetwork;
	int statistic = 0;
	
	for (TieIterator iter = pLostTieNetwork->ties(); iter.valid(); iter.next())
	{
		statistic +=
			pOneModeNetwork->twoPathCount(iter.alter(), iter.ego());
	}
	
	return statistic / 3.0;
}

}
