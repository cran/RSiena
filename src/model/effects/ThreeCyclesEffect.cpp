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
#include "network/OneModeNetwork.h"
#include "network/TieIterator.h"
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
double ThreeCyclesEffect::calculateContribution(int alter) const
{
	// The absolute value of the contribution is going to be the number
	// of two-paths from the alter to the ego -- a number, which is stored
	// in the table of reverse two paths.

	return this->pReverseTwoPathTable()->get(alter);
}


/**
 * See base class.
 */
double ThreeCyclesEffect::statistic(const Network * pSummationTieNetwork) const
{
	OneModeNetwork * pOneModeNetwork = (OneModeNetwork *) this->pNetwork();
	int statistic = 0;

	for (TieIterator iter = pSummationTieNetwork->ties();
		iter.valid();
		iter.next())
	{
		statistic +=
			pOneModeNetwork->twoPathCount(iter.alter(), iter.ego());
	}

	return statistic / 3.0;
}

}
