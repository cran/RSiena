/******************************************************************************
 * SIENA: Simulation Investigation for Empirical Network Analysis
 *
 * Web: http://www.stats.ox.ac.uk/~snijders/siena/
 *
 * File: OutTruncEffect.h
 *
 * Description: This file contains the definition of the
 * OutTruncEffect class.
 *****************************************************************************/

#ifndef OUTTRUNCEFFECT_H_
#define OUTTRUNCEFFECT_H_

#include "NetworkEffect.h"

namespace siena
{

/**
 * This class defines the outdegree activity effect defined by
 * s_i(x)= min{x_{i+},c}. The corresponding statistic is
 * the sum of outdegrees truncated at c over all actors.
 */
class OutTruncEffect : public NetworkEffect
{
public:
	OutTruncEffect(const EffectInfo * pEffectInfo);

	virtual double calculateContribution(int alter) const;
	virtual double endowmentStatistic(Network * pLostTieNetwork);

protected:
	virtual double egoStatistic(int ego,
		const Network * pSummationTieNetwork);

private:
	double lc;
};

}

#endif /*OUTTRUNCEFFECT_H_*/
