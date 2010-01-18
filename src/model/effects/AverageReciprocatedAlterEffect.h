/******************************************************************************
 * SIENA: Simulation Investigation for Empirical Network Analysis
 *
 * Web: http://www.stats.ox.ac.uk/~snijders/siena/
 *
 * File: AverageReciprocatedAlterEffect.h
 *
 * Description: This file contains the definition of the
 * AverageReciprocatedAlterEffect class.
 *****************************************************************************/

#ifndef AVERAGERECIPROCATEDALTEREFFECT_H_
#define AVERAGERECIPROCATEDALTEREFFECT_H_

#include "NetworkDependentBehaviorEffect.h"

namespace siena
{

/**
 * Average reciprocated alter effect (see manual).
 */
class AverageReciprocatedAlterEffect : public NetworkDependentBehaviorEffect
{
public:
	AverageReciprocatedAlterEffect(const EffectInfo * pEffectInfo);

	virtual double calculateChangeContribution(int actor,
		int difference);
	virtual double egoStatistic(int ego, double * currentValues);
};

}

#endif /*AVERAGERECIPROCATEDALTEREFFECT_H_*/
