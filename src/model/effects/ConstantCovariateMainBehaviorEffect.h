/******************************************************************************
 * SIENA: Simulation Investigation for Empirical Network Analysis
 *
 * Web: http://www.stats.ox.ac.uk/~snijders/siena/
 *
 * File: ConstantCovariateMainBehaviorEffect.h
 *
 * Description: This file contains the definition of the
 * ConstantCovariateMainBehaviorEffect class.
 *****************************************************************************/

#ifndef CONSTANTCOVARIATEMAINBEHAVIOREFFECT_H_
#define CONSTANTCOVARIATEMAINBEHAVIOREFFECT_H_

#include "ConstantCovariateBehaviorEffect.h"

namespace siena
{

/**
 * Constant covariate main behavior effect defined as the product of the ego 
 * with the covariate.
 */
	class ConstantCovariateMainBehaviorEffect : public ConstantCovariateBehaviorEffect
{
public:
	ConstantCovariateMainBehaviorEffect(const EffectInfo * pEffectInfo);

	virtual double calculateChangeContribution(int actor,
		int difference) const;
	virtual double evaluationStatistic(double * currentValues) const;
	virtual double endowmentStatistic(const int * difference,
		double * currentValues) const;
};

}

#endif /*CONSTANTCOVARIATEMAINBEHAVIOREFFECT_H_*/
