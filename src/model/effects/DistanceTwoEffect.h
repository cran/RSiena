/******************************************************************************
 * SIENA: Simulation Investigation for Empirical Network Analysis
 *
 * Web: http://www.stats.ox.ac.uk/~snijders/siena/
 *
 * File: DistanceTwoEffect.h
 *
 * Description: This file contains the declaration of the class
 * DistanceTwoEffect.
 *****************************************************************************/

#ifndef DISTANCETWOEFFECT_H_
#define DISTANCETWOEFFECT_H_

#include "model/effects/NetworkEffect.h"

namespace siena
{

/**
 * This class defines a generalization of the "number of distances two" effect,
 * where one can specify the number of two-paths required between a pair of
 * actors to qualify for a distance two pair. The manual defines the effect for
 * parameter values 1 and 2.
 */
class DistanceTwoEffect : public NetworkEffect
{
public:
	DistanceTwoEffect(const EffectInfo * pEffectInfo,
		int requiredTwoPathCount);

	virtual double calculateTieFlipContribution(int alter) const;
	virtual bool usesTable(const ConfigurationTable * pTable) const;
	virtual double evaluationStatistic(Network * pNetwork) const;
	virtual double endowmentStatistic(Network * pInitialNetwork,
		Network * pLostTieNetwork) const;

private:
	int lrequiredTwoPathCount;
};

}

#endif /*DISTANCETWOEFFECT_H_*/
