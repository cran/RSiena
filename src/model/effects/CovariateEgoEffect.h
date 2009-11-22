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

	virtual double calculateContribution(int alter) const;

protected:
	virtual double statistic(const Network * pSummationTieNetwork) const;
};

}

#endif /*COVARIATEEGOEFFECT_H_*/
