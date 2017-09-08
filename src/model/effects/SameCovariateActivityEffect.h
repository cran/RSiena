/******************************************************************************
 * SIENA: Simulation Investigation for Empirical Network Analysis
 *
 * Web: http://www.stats.ox.ac.uk/~snijders/siena/
 *
 * File: SameCovariateActivityEffect.h
 *
 * Description: This file contains the definition of the
 * SameCovariateActivityEffect class.
 *****************************************************************************/

#ifndef SAMECOVARIATEACTIVITYEFFECT_H_
#define SAMECOVARIATEACTIVITYEFFECT_H_

#include "CovariateDependentNetworkEffect.h"

namespace siena
{

/**
 * Same and different covariate activity effects (see manual).
 */
class SameCovariateActivityEffect : public CovariateDependentNetworkEffect
{
public:
	SameCovariateActivityEffect(const EffectInfo * pEffectInfo, bool same);

	virtual double calculateContribution(int alter) const;

protected:
	virtual double tieStatistic(int alter);

private:
	bool lsame;
};

}

#endif /*SAMECOVARIATEACTIVITYEFFECT_H_*/
