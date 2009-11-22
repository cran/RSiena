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
#include "network/OneModeNetwork.h"
#include "network/TieIterator.h"
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
double TransitiveTriadsEffect::calculateContribution(int alter) const
{
	return this->pTwoPathTable()->get(alter);
}


/**
 * See base class.
 */
double TransitiveTriadsEffect::statistic(const Network * pSummationTieNetwork)
	const
{
	OneModeNetwork * pOneModeNetwork = (OneModeNetwork *) this->pNetwork();
	int counter = 0;

	for (TieIterator iter = pSummationTieNetwork->ties();
		iter.valid();
		iter.next())
	{
		counter +=
			pOneModeNetwork->twoPathCount(iter.ego(), iter.alter());
	}

	return counter / 6.0;
}

}
