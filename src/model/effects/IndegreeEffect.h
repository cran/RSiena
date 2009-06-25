/******************************************************************************
 * SIENA: Simulation Investigation for Empirical Network Analysis
 *
 * Web: http://www.stats.ox.ac.uk/~snijders/siena/
 *
 * File: IndegreerEffect.h
 *
 * Description: This file contains the definition of the
 * IndegreeEffect class.
 *****************************************************************************/

#ifndef INDEGREEEFFECT_H_
#define INDEGREEEFFECT_H_

#include "NetworkDependentBehaviorEffect.h"

namespace siena
{

/**
 * Indegree effect defined as the product of the ego with its
 * indegree (with respect to a certain network).
 */
	class IndegreeEffect : public NetworkDependentBehaviorEffect
{
public:
	IndegreeEffect(const EffectInfo * pEffectInfo);

	virtual double calculateChangeContribution(int actor,
		int difference) const;
	virtual double evaluationStatistic(double * currentValues) const;
	virtual double endowmentStatistic(const int * difference,
		double * currentValues) const;
};

}

#endif /*INDEGREEEFFECT_H_*/

