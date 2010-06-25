/******************************************************************************
 * SIENA: Simulation Investigation for Empirical Network Analysis
 *
 * Web: http://www.stats.ox.ac.uk/~snijders/siena/
 *
 * File: CovariateDistance2AlterEffect.cpp
 *
 * Description: This file contains the implementation of the
 * CovariateDistance2AlterEffect class.
 *****************************************************************************/
#include "CovariateDistance2AlterEffect.h"
#include "network/Network.h"
#include "model/variables/NetworkVariable.h"

namespace siena
{

/**
 * Constructor.
 * @param[in] pEffectInfo the effect descriptor
 */
CovariateDistance2AlterEffect::CovariateDistance2AlterEffect(const EffectInfo * pEffectInfo) :
		CovariateDistance2NetworkEffect(pEffectInfo)
{
}


/**
 * Calculates the contribution of a tie flip to the given actor.
 */
double CovariateDistance2AlterEffect::calculateContribution(int alter) const
{
	double change = this->averageAlterValue(alter);
	return change;
}


/**
 * The contribution of the tie from the implicit ego to the given alter
 * to the statistic. It is assumed that preprocessEgo(ego) has been
 * called before.
 */
double CovariateDistance2AlterEffect::tieStatistic(int alter)
{
	double statistic = 0;

	if (!this->missingDummy(alter))
	{
		statistic = this->averageAlterValue(alter);

	}
	return statistic;
}

}
