/******************************************************************************
 * SIENA: Simulation Investigation for Empirical Network Analysis
 *
 * Web: http://www.stats.ox.ac.uk/~snijders/siena/
 *
 * File: CovariateSimilarityEffect.h
 *
 * Description: This file contains the definition of the
 * CovariateSimilarityEffect class.
 *****************************************************************************/

#ifndef COVARIATESIMILARITYEFFECT_H_
#define COVARIATESIMILARITYEFFECT_H_

#include "CovariateDependentNetworkEffect.h"

namespace siena
{

/**
 * Covariate-related similarity and covariate-related similarity x reciprocity
 * effects (see manual).
 */
class CovariateSimilarityEffect : public CovariateDependentNetworkEffect
{
public:
	CovariateSimilarityEffect(const EffectInfo * pEffectInfo,
		bool reciprocal);

	virtual double calculateTieFlipContribution(int alter) const;

protected:
	virtual double statistic(Network * pNetwork,
		Network * pSummationTieNetwork) const;

private:
	// Indicates if the reciprocal version of the similarity effect is
	// required

	bool lreciprocal;
};

}

#endif /*COVARIATESIMILARITYEFFECT_H_*/
