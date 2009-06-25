/******************************************************************************
 * SIENA: Simulation Investigation for Empirical Network Analysis
 * 
 * Web: http://www.stats.ox.ac.uk/~snijders/siena/
 * 
 * File: DyadicCovariateMainEffect.h
 * 
 * Description: This file contains the definition of the
 * DyadicCovariateMainEffect class.
 *****************************************************************************/

#ifndef DYADICCOVARIATEMAINEFFECT_H_
#define DYADICCOVARIATEMAINEFFECT_H_

#include "DyadicCovariateDependentNetworkEffect.h"

namespace siena
{

/**
 * Dyadic covariate main effect (see manual).
 */
class DyadicCovariateMainEffect : public DyadicCovariateDependentNetworkEffect
{
public:
	DyadicCovariateMainEffect(const EffectInfo * pEffectInfo);
	
	virtual double calculateTieFlipContribution(int alter) const;
	virtual double evaluationStatistic(Network * pNetwork) const;
	virtual double endowmentStatistic(Network * pInitialNetwork,
		Network * pLostTieNetwork) const;
};

}

#endif /*DYADICCOVARIATEMAINEFFECT_H_*/
