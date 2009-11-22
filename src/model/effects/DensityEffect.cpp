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
#include <stdexcept>

#include "DensityEffect.h"
#include "network/Network.h"
#include "model/variables/NetworkVariable.h"
#include "data/OneModeNetworkLongitudinalData.h"

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
double DensityEffect::calculateContribution(int alter) const
{
	return 1;
}


/**
 * See base class for a detailed comment.
 */
double DensityEffect::statistic(const Network * pSummationTieNetwork) const
{
	double statistic = pSummationTieNetwork->tieCount();
	const OneModeNetworkLongitudinalData * pData =
		dynamic_cast<const OneModeNetworkLongitudinalData *>(this->pData());

	if (pData)
	{
		if (pData->symmetric())
		{
			statistic /= 2;
		}
	}

	return statistic;
}

}
