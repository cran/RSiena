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
#include "network/OneModeNetwork.h"
#include "network/TieIterator.h"
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
double TransitiveMediatedTripletsEffect::calculateContribution(
	int alter) const
{
	// The number of out-stars to the ego and the given alter is the amount
	// of change for this tie flip.

	return this->pOutStarTable()->get(alter);
}


/**
 * See base class.
 */
double TransitiveMediatedTripletsEffect::statistic(
	const Network * pSummationTieNetwork) const
{
	OneModeNetwork * pOneModeNetwork = (OneModeNetwork *) this->pNetwork();
	int statistic = 0;

	for (TieIterator iter = pSummationTieNetwork->ties();
		iter.valid();
		iter.next())
	{
		statistic +=
			pOneModeNetwork->outTwoStarCount(iter.ego(), iter.alter());
	}

	return statistic;
}

}
