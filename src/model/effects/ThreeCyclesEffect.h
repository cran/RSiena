/******************************************************************************
 * SIENA: Simulation Investigation for Empirical Network Analysis
 * 
 * Web: http://www.stats.ox.ac.uk/~snijders/siena/
 * 
 * File: ThreeCyclesEffect.h
 * 
 * Description: This file contains the declaration of the class
 * ThreeCyclesEffect.
 *****************************************************************************/

#ifndef THREECYCLESEFFECT_H_
#define THREECYCLESEFFECT_H_

#include "model/effects/NetworkEffect.h"

namespace siena
{

/**
 * This class defines the 3-cycles effect.
 */
class ThreeCyclesEffect : public NetworkEffect
{
public:
	ThreeCyclesEffect(const EffectInfo * pEffectInfo);

	virtual double calculateTieFlipContribution(int alter) const;
	virtual bool usesTable(const ConfigurationTable * pTable) const;
	virtual double evaluationStatistic(Network * pNetwork) const;
	virtual double endowmentStatistic(Network * pInitialNetwork,
		Network * pLostTieNetwork) const;
};

}

#endif /*THREECYCLESEFFECT_H_*/
