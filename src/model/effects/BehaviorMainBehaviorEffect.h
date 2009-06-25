/******************************************************************************
 * SIENA: Simulation Investigation for Empirical Network Analysis
 *
 * Web: http://www.stats.ox.ac.uk/~snijders/siena/
 *
 * File: BehaviorMainBehaviorEffect.h
 *
 * Description: This file contains the definition of the
 * BehaviorMainBehaviorEffect class.
 *****************************************************************************/

#ifndef BEHAVIORMAINBEHAVIOREFFECT_H_
#define BEHAVIORMAINBEHAVIOREFFECT_H_

#include "BehaviorDependentBehaviorEffect.h"

namespace siena
{

/**
 * Behavior variable main behavior effect defined as the product of the ego 
 * with the covariate.
 */
	class BehaviorMainBehaviorEffect : public BehaviorDependentBehaviorEffect
{
public:
	BehaviorMainBehaviorEffect(const EffectInfo * pEffectInfo);

	virtual double calculateChangeContribution(int actor,
		int difference) const;
	virtual double evaluationStatistic(double * currentValues) const;
	virtual double endowmentStatistic(const int * difference,
		double * currentValues) const;
};

}

#endif /*BEHAVIORMAINBEHAVIOREFFECT_H_*/
