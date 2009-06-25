/******************************************************************************
 * SIENA: Simulation Investigation for Empirical Network Analysis
 *
 * Web: http://www.stats.ox.ac.uk/~snijders/siena/
 *
 * File: ChangingCovariateMainBehaviorEffect.h
 *
 * Description: This file contains the definition of the
 * ChangingCovariateMainBehaviorEffect class.
 *****************************************************************************/

#ifndef CHANGINGCOVARIATEMAINBEHAVIOREFFECT_H_
#define CHANGINGCOVARIATEMAINBEHAVIOREFFECT_H_

#include "ChangingCovariateBehaviorEffect.h"

namespace siena
{

/**
 * Changing covariate main behavior effect defined as the product of the ego 
 * with the covariate.
 */
	class ChangingCovariateMainBehaviorEffect : public ChangingCovariateBehaviorEffect
{
public:
	ChangingCovariateMainBehaviorEffect(const EffectInfo * pEffectInfo);

	virtual double calculateChangeContribution(int actor,
		int difference) const;
	virtual double evaluationStatistic(double * currentValues) const;
	virtual double endowmentStatistic(const int * difference,
		double * currentValues) const;
};

}

#endif /*CHANGINGCOVARIATEMAINBEHAVIOREFFECT_H_*/
