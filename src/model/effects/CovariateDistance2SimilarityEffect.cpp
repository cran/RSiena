/******************************************************************************
 * SIENA: Simulation Investigation for Empirical Network Analysis
 *
 * Web: http://www.stats.ox.ac.uk/~snijders/siena/
 *
 * File: CovariateDistance2SimilarityEffect.cpp
 *
 * Description: This file contains the implementation of the
 * CovariateDistance2SimilarityEffect class.
 *****************************************************************************/


#include "CovariateDistance2SimilarityEffect.h"
#include "data/NetworkLongitudinalData.h"
#include "network/Network.h"
#include "network/IncidentTieIterator.h"
#include "network/CommonNeighborIterator.h"
#include "model/variables/NetworkVariable.h"

namespace siena
{

/**
 * Constructor.
 * @param[in] pEffectInfo the effect descriptor
 */
CovariateDistance2SimilarityEffect::CovariateDistance2SimilarityEffect(
	const EffectInfo * pEffectInfo) :
		CovariateDistance2NetworkEffect(pEffectInfo)
{
}


/**
 * Calculates the contribution of a tie flip to the given actor.
 */
double CovariateDistance2SimilarityEffect::calculateContribution(int alter) const
{
	double change = 0;

	change = this->similarityNetwork(this->ego(), alter,
		this->pData()->name());

	return change;
}


/**
 * The contribution of the tie from the implicit ego to the given alter
 * to the statistic. It is assumed that preprocessEgo(ego) has been
 * called before.
 */
double CovariateDistance2SimilarityEffect::tieStatistic(int alter)
{
	double statistic = 0;

	if (!this->missingDummy(this->ego()) && !this->missingDummy(alter) )
	{
		statistic = this->similarityNetwork(this->ego(), alter,
			this->pData()->name());
	}
	return statistic;
}

}
