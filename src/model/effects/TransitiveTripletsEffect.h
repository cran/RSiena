/******************************************************************************
 * SIENA: Simulation Investigation for Empirical Network Analysis
 * 
 * Web: http://www.stats.ox.ac.uk/~snijders/siena/
 * 
 * File: TransitiveTripletsEffect.h
 * 
 * Description: This file contains the declaration of the class
 * TransitiveTripletsEffect.
 *****************************************************************************/

#ifndef TRANSITIVETRIPLETSEFFECT_H_
#define TRANSITIVETRIPLETSEFFECT_H_

#include "model/effects/NetworkEffect.h"

namespace siena
{

/**
 * This class defines the transitive triplets effect for non-symmetric
 * networks.
 */
class TransitiveTripletsEffect : public NetworkEffect
{
public:
	TransitiveTripletsEffect(const EffectInfo * pEffectInfo);

	virtual double calculateTieFlipContribution(int alter) const;
	virtual bool usesTable(const ConfigurationTable * pTable) const;
	virtual double evaluationStatistic(Network * pNetwork) const;
	virtual double endowmentStatistic(Network * pInitialNetwork,
		Network * pLostTieNetwork) const;
};

}

#endif /*TRANSITIVETRIPLETSEFFECT_H_*/
