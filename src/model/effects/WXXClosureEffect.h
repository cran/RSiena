/******************************************************************************
 * SIENA: Simulation Investigation for Empirical Network Analysis
 *
 * Web: http://www.stats.ox.ac.uk/~snijders/siena/
 *
 * File: WXXClosureEffect.h
 *
 * Description: This file contains the definition of the
 * WXXClosureEffect class.
 *****************************************************************************/

#ifndef WXXCLOSUREEFFECT_H_
#define WXXCLOSUREEFFECT_H_

#include "DyadicCovariateDependentNetworkEffect.h"

namespace siena
{

/**
 * WX => X closure of covariate effect (see manual).
 */
class WXXClosureEffect : public DyadicCovariateDependentNetworkEffect
{
public:
	WXXClosureEffect(const EffectInfo * pEffectInfo);
	virtual ~WXXClosureEffect();

	virtual void initialize(EpochSimulation * pSimulation);
	virtual void initialize(const Data * pData, State * pState, int period);

	virtual void preprocessEgo();
	virtual double calculateTieFlipContribution(int alter) const;

protected:
	virtual double statistic(Network * pNetwork,
		Network * pSummationTieNetwork) const;

private:
	void calculateSums(int i, Network * pNetwork, double * sums) const;

	// For a fixed i, this variable stores the value of sum_h w_{ih} x_{hj} for
	// each j.

	double * lsums;
};

}

#endif /*WXXCLOSUREEFFECT_H_*/
