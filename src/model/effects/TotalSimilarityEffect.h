/******************************************************************************
 * SIENA: Simulation Investigation for Empirical Network Analysis
 *
 * Web: http://www.stats.ox.ac.uk/~snijders/siena/
 *
 * File: TotalSimilarityEffect.h
 *
 * Description: This file contains the definition of the
 * TotalSimilarityEffect class.
 *****************************************************************************/

#ifndef TOTALSIMILARITYEFFECT_H_
#define TOTALSIMILARITYEFFECT_H_

#include "NetworkDependentBehaviorEffect.h"

namespace siena
{

/**
 * Total similarity effect defined as the sum of centered similarity
 * of the ego with each of its neighbors (with respect to a certain network).
 */
	class TotalSimilarityEffect : public NetworkDependentBehaviorEffect
{
public:
	TotalSimilarityEffect(const EffectInfo * pEffectInfo);

	virtual double calculateChangeContribution(int actor,
		int difference) const;
	virtual double evaluationStatistic(double * currentValues) const;
	virtual double endowmentStatistic(const int * difference,
		double * currentValues) const;
};

}

#endif /*TOTALSIMILARITYEFFECT_H_*/
