/******************************************************************************
 * SIENA: Simulation Investigation for Empirical Network Analysis
 *
 * Web: http://www.stats.ox.ac.uk/~snijders/siena/
 *
 * File: AverageSimilarityEffect.h
 *
 * Description: This file contains the definition of the
 * AverageSimilarityEffect class.
 *****************************************************************************/

#ifndef AVERAGESIMILARITYEFFECT_H_
#define AVERAGESIMILARITYEFFECT_H_

#include "NetworkDependentBehaviorEffect.h"

namespace siena
{

/**
 * Average similarity effect defined as the average centered similarity
 * of the ego with each of its neighbors (with respect to a certain network).
 */
class AverageSimilarityEffect : public NetworkDependentBehaviorEffect
{
public:
	AverageSimilarityEffect(const EffectInfo * pEffectInfo);

	virtual double calculateChangeContribution(int actor,
		int difference) const;
	virtual double evaluationStatistic(double * currentValues) const;
	virtual double endowmentStatistic(const int * difference,
		double * currentValues) const;
};

}

#endif /*AVERAGESIMILARITYEFFECT_H_*/
