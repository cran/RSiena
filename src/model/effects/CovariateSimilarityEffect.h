/******************************************************************************
 * SIENA: Simulation Investigation for Empirical Network Analysis
 * 
 * Web: http://www.stats.ox.ac.uk/~snijders/siena/
 * 
 * File: CovariateSimilarityEffect.h
 * 
 * Description: This file contains the definition of the
 * CovariateSimilarityEffect class.
 *****************************************************************************/

#ifndef COVARIATESIMILARITYEFFECT_H_
#define COVARIATESIMILARITYEFFECT_H_

#include "CovariateDependentNetworkEffect.h"

namespace siena
{

/**
 * Covariate-related similarity effect (see manual).
 */
class CovariateSimilarityEffect : public CovariateDependentNetworkEffect
{
public:
	CovariateSimilarityEffect(const EffectInfo * pEffectInfo);
	
	virtual double calculateTieFlipContribution(int alter) const;
	virtual double evaluationStatistic(Network * pNetwork) const;
	virtual double endowmentStatistic(Network * pInitialNetwork,
		Network * pLostTieNetwork) const;
};

}

#endif /*COVARIATESIMILARITYEFFECT_H_*/
