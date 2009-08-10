/******************************************************************************
 * SIENA: Simulation Investigation for Empirical Network Analysis
 *
 * Web: http://www.stats.ox.ac.uk/~snijders/siena/
 *
 * File: DyadicCovariateReciprocityEffect.h
 *
 * Description: This file contains the definition of the
 * DyadicCovariateReciprocityEffect class.
 *****************************************************************************/

#ifndef DYADICCOVARIATERECIPROCITYEFFECT_H_
#define DYADICCOVARIATERECIPROCITYEFFECT_H_

#include "DyadicCovariateDependentNetworkEffect.h"

namespace siena
{

/**
 * Dyadic covariate x reciprocity effect (see manual).
 */
class DyadicCovariateReciprocityEffect :
	public DyadicCovariateDependentNetworkEffect
{
public:
	DyadicCovariateReciprocityEffect(const EffectInfo * pEffectInfo);

	virtual double calculateTieFlipContribution(int alter) const;

protected:
	virtual double statistic(Network * pNetwork,
		Network * pSummationTieNetwork) const;
};

}

#endif /*DYADICCOVARIATERECIPROCITYEFFECT_H_*/
