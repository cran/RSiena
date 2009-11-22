/******************************************************************************
 * SIENA: Simulation Investigation for Empirical Network Analysis
 *
 * Web: http://www.stats.ox.ac.uk/~snijders/siena/
 *
 * File: WWXClosureEffect.h
 *
 * Description: This file contains the definition of the
 * WWXClosureEffect class.
 *****************************************************************************/

#ifndef WWXCLOSUREEFFECT_H_
#define WWXCLOSUREEFFECT_H_

#include "DyadicCovariateDependentNetworkEffect.h"

namespace siena
{

/**
 * WW => X closure of covariate effect (see manual).
 */
class WWXClosureEffect : public DyadicCovariateDependentNetworkEffect
{
public:
	WWXClosureEffect(const EffectInfo * pEffectInfo);
	virtual ~WWXClosureEffect();

	virtual void initialize(const Data * pData,
		State * pState,
		int period,
		Cache * pCache);

	virtual void preprocessEgo(int ego);
	virtual double calculateContribution(int alter) const;

protected:
	virtual double statistic(const Network * pSummationTieNetwork) const;

private:
	void calculateSums(int i, const Network * pNetwork, double * sums) const;

	// For a fixed i, this variable stores the value of sum_h w_{ih} w_{hj} for
	// each j.

	double * lsums;
};

}

#endif /*WWXCLOSUREEFFECT_H_*/
