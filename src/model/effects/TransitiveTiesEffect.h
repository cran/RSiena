/******************************************************************************
 * SIENA: Simulation Investigation for Empirical Network Analysis
 * 
 * Web: http://www.stats.ox.ac.uk/~snijders/siena/
 * 
 * File: TransitiveTiesEffect.h
 * 
 * Description: This file contains the declaration of the class
 * TransitiveTiesEffect.
 *****************************************************************************/

#ifndef TRANSITIVETIESEFFECT_H_
#define TRANSITIVETIESEFFECT_H_

#include "model/effects/NetworkEffect.h"

namespace siena
{

/**
 * This class defines the transitive ties effect.
 */
class TransitiveTiesEffect : public NetworkEffect
{
public:
	TransitiveTiesEffect(const EffectInfo * pEffectInfo);

	virtual double calculateTieFlipContribution(int alter) const;
	virtual bool usesTable(const ConfigurationTable * pTable) const;
	virtual double evaluationStatistic(Network * pNetwork) const;
	virtual double endowmentStatistic(Network * pInitialNetwork,
		Network * pLostTieNetwork) const;
};

}

#endif /*TRANSITIVETIESEFFECT_H_*/
