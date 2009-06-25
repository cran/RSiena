/******************************************************************************
 * SIENA: Simulation Investigation for Empirical Network Analysis
 * 
 * Web: http://www.stats.ox.ac.uk/~snijders/siena/
 * 
 * File: SameCovariateEffect.h
 * 
 * Description: This file contains the definition of the
 * SameCovariateEffect class.
 *****************************************************************************/

#ifndef SAMECOVARIATEEFFECT_H_
#define SAMECOVARIATEEFFECT_H_

#include "CovariateDependentNetworkEffect.h"

namespace siena
{

/**
 * Same covariate effect (see manual).
 */
class SameCovariateEffect : public CovariateDependentNetworkEffect
{
public:
	SameCovariateEffect(const EffectInfo * pEffectInfo);
	
	virtual double calculateTieFlipContribution(int alter) const;
	virtual double evaluationStatistic(Network * pNetwork) const;
	virtual double endowmentStatistic(Network * pInitialNetwork,
		Network * pLostTieNetwork) const;
};

}

#endif /*SAMECOVARIATEEFFECT_H_*/
