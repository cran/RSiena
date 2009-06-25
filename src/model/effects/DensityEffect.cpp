/******************************************************************************
 * SIENA: Simulation Investigation for Empirical Network Analysis
 * 
 * Web: http://www.stats.ox.ac.uk/~snijders/siena/
 * 
 * File: DensityEffect.cpp
 * 
 * Description: This file contains the implementation of the class
 * DensityEffect.
 *****************************************************************************/

#include "DensityEffect.h"
#include "data/Network.h"
#include "model/variables/NetworkVariable.h"

namespace siena
{

/**
 * Constructor.
 */
DensityEffect::DensityEffect(const EffectInfo * pEffectInfo) :
	NetworkEffect(pEffectInfo)
{
}


/**
 * Calculates the contribution of a tie flip to the given actor.
 */
double DensityEffect::calculateTieFlipContribution(int alter) const
{
	double change = 1;
	
	if (this->pVariable()->outTieExists(alter))
	{
		// The ego would loose one tie
		change = -1;
	}
	
	return change;
}


/**
 * Returns the statistic corresponding to this effect as part of
 * the evaluation function with respect to the given network.
 */
double DensityEffect::evaluationStatistic(Network * pNetwork) const
{
	return pNetwork->tieCount();
}


/**
 * Returns the statistic corresponding to this effect as part of
 * the endowment function with respect to an initial network
 * and a network of lost ties. The current network is implicit as
 * the introduced ties are not relevant for calculating
 * endowment statistics.
 */
double DensityEffect::endowmentStatistic(Network * pInitialNetwork,
	Network * pLostTieNetwork) const
{
	return pLostTieNetwork->tieCount();
}

}
