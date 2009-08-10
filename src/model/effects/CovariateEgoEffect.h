/******************************************************************************
 * SIENA: Simulation Investigation for Empirical Network Analysis
 * 
 * Web: http://www.stats.ox.ac.uk/~snijders/siena/
 * 
 * File: CovariateEgoEffect.h
 * 
 * Description: This file contains the definition of the
 * CovariateEgoEffect class.
 *****************************************************************************/

#ifndef COVARIATEEGOEFFECT_H_
#define COVARIATEEGOEFFECT_H_

#include "CovariateDependentNetworkEffect.h"

namespace siena
{

/**
 * Covariate-ego or covariate-related activity effect (see manual).
 */
class CovariateEgoEffect : public CovariateDependentNetworkEffect
{
public:
	CovariateEgoEffect(const EffectInfo * pEffectInfo);
	
	virtual double calculateTieFlipContribution(int alter) const;
	virtual double evaluationStatistic(Network * pNetwork) const;
	virtual double endowmentStatistic(Network * pInitialNetwork,
		Network * pLostTieNetwork) const;
};

}

#endif /*COVARIATEEGOEFFECT_H_*/
