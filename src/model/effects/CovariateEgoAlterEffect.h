/******************************************************************************
 * SIENA: Simulation Investigation for Empirical Network Analysis
 * 
 * Web: http://www.stats.ox.ac.uk/~snijders/siena/
 * 
 * File: CovariateEgoAlterEffect.h
 * 
 * Description: This file contains the definition of the
 * CovariateEgoAlterEffect class.
 *****************************************************************************/

#ifndef COVARIATEEGOALTEREFFECT_H_
#define COVARIATEEGOALTEREFFECT_H_

#include "CovariateDependentNetworkEffect.h"

namespace siena
{

/**
 * Covariate-ego x alter effect (see manual).
 */
class CovariateEgoAlterEffect : public CovariateDependentNetworkEffect
{
public:
	CovariateEgoAlterEffect(const EffectInfo * pEffectInfo);
	
	virtual double calculateTieFlipContribution(int alter) const;
	virtual double evaluationStatistic(Network * pNetwork) const;
	virtual double endowmentStatistic(Network * pInitialNetwork,
		Network * pLostTieNetwork) const;
};

}

#endif /*COVARIATEEGOALTEREFFECT_H_*/
