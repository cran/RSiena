/******************************************************************************
 * SIENA: Simulation Investigation for Empirical Network Analysis
 *
 * Web: http://www.stats.ox.ac.uk/~snijders/siena/
 *
 * File: XWXClosureEffect.h
 *
 * Description: This file contains the definition of the
 * XWXClosureEffect class.
 *****************************************************************************/

#ifndef XWXCLOSUREEFFECT_H_
#define XWXCLOSUREEFFECT_H_

#include "DyadicCovariateDependentNetworkEffect.h"

namespace siena
{

/**
 * XW => X closure of covariate effect (see manual).
 */
class XWXClosureEffect : public DyadicCovariateDependentNetworkEffect
{
public:
	XWXClosureEffect(const EffectInfo * pEffectInfo);
	virtual ~XWXClosureEffect();

	virtual void initialize(const Data * pData,
		State * pState,
		int period,
		Cache * pCache);

	virtual void preprocessEgo(int ego);
	virtual double calculateContribution(int alter) const;

protected:
	virtual double statistic(const Network * pSummationTieNetwork) const;

private:
	void calculateTwoPathSums(int i,
		const Network * pNetwork,
		double * sums) const;
	void calculateInStarSums(int i,
		const Network * pNetwork,
		double * sums) const;

	// For a fixed i, this variable stores the value of sum_h x_{ih} w_{hj} for
	// each j.

	double * ltwoPathSums;

	// For a fixed i, this variable stores the value of sum_h x_{ih} w_{jh} for
	// each j.

	double * linStarSums;
};

}

#endif /*XWXCLOSUREEFFECT_H_*/
