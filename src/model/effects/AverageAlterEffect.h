/******************************************************************************
 * SIENA: Simulation Investigation for Empirical Network Analysis
 *
 * Web: http://www.stats.ox.ac.uk/~snijders/siena/
 *
 * File: AverageAlterEffect.h
 *
 * Description: This file contains the definition of the
 * AverageAlterEffect class.
 *****************************************************************************/

#ifndef AVERAGEALTEREFFECT_H_
#define AVERAGEALTEREFFECT_H_

#include "NetworkDependentBehaviorEffect.h"

namespace siena
{

/**
 * Average alter effect defined as the product of the ego with the average
 * of its neighbors (with respect to a certain network).
 */
class AverageAlterEffect : public NetworkDependentBehaviorEffect
{
public:
	AverageAlterEffect(const EffectInfo * pEffectInfo);

	virtual double calculateChangeContribution(int actor,
		int difference);
	virtual double endowmentStatistic(const int * difference,
		double * currentValues);
	virtual double egoStatistic(int ego, double * currentValues);
};

}

#endif /*AVERAGEALTEREFFECT_H_*/
