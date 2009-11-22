/******************************************************************************
 * SIENA: Simulation Investigation for Empirical Network Analysis
 *
 * Web: http://www.stats.ox.ac.uk/~snijders/siena/
 *
 * File: CovariateIndirectTiesEffect.h
 *
 * Description: This file contains the declaration of the class
 * CovariateIndirectTiesEffect.
 *****************************************************************************/

#ifndef COVARIATEINDIRECTTIESEFFECT_H_
#define COVARIATEINDIRECTTIESEFFECT_H_

#include "model/effects/CovariateDependentNetworkEffect.h"

namespace siena
{

/**
 * Covariate of indirect ties effect (see manual).
 */
class CovariateIndirectTiesEffect : public CovariateDependentNetworkEffect
{
public:
	CovariateIndirectTiesEffect(const EffectInfo * pEffectInfo);

	virtual double calculateContribution(int alter) const;
	virtual double evaluationStatistic() const;
	virtual double endowmentStatistic(Network * pLostTieNetwork) const;
};

}

#endif /*COVARIATEINDIRECTTIESEFFECT_H_*/
