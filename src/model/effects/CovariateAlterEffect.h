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
 * Covariate-alter and covariate squared-alter effects (see manual).
 */
class CovariateAlterEffect : public CovariateDependentNetworkEffect
{
public:
	CovariateAlterEffect(const EffectInfo * pEffectInfo, bool squared);

	virtual double calculateTieFlipContribution(int alter) const;
	virtual double evaluationStatistic(Network * pNetwork) const;
	virtual double endowmentStatistic(Network * pInitialNetwork,
		Network * pLostTieNetwork) const;

private:
	bool lsquared;
};

}

#endif /*COVARIATEALTEREFFECT_H_*/
