/******************************************************************************
 * SIENA: Simulation Investigation for Empirical Network Analysis
 * 
 * Web: http://www.stats.ox.ac.uk/~snijders/siena/
 * 
 * File: CovariateAlterEffect.h
 * 
 * Description: This file contains the definition of the
 * CovariateAlterEffect class.
 *****************************************************************************/

#ifndef COVARIATEALTEREFFECT_H_
#define COVARIATEALTEREFFECT_H_

#include "CovariateDependentNetworkEffect.h"

namespace siena
{

/**
 * Covariate-alter or covariate-related popularity effect (see manual).
 */
class CovariateAlterEffect : public CovariateDependentNetworkEffect
{
public:
	CovariateAlterEffect(const EffectInfo * pEffectInfo);
	
	virtual double calculateTieFlipContribution(int alter) const;
	virtual double evaluationStatistic(Network * pNetwork) const;
	virtual double endowmentStatistic(Network * pInitialNetwork,
		Network * pLostTieNetwork) const;
};

}

#endif /*COVARIATEALTEREFFECT_H_*/
