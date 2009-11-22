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

	virtual double calculateContribution(int alter) const;

protected:
	virtual double statistic(const Network * pSummationTieNetwork) const;

private:
	bool lsquared;
};

}

#endif /*COVARIATEALTEREFFECT_H_*/
