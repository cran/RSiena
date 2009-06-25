/******************************************************************************
 * SIENA: Simulation Investigation for Empirical Network Analysis
 * 
 * Web: http://www.stats.ox.ac.uk/~snijders/siena/
 * 
 * File: ReciprocityEffect.h
 * 
 * Description: This file contains the declaration of the class
 * ReciprocityEffect.
 *****************************************************************************/

#ifndef RECIPROCITYEFFECT_H_
#define RECIPROCITYEFFECT_H_

#include "model/effects/NetworkEffect.h"

namespace siena
{

/**
 * This class defines the reciprocity effect. This effect makes sense for
 * non-symmetric networks only.
 */
class ReciprocityEffect : public NetworkEffect
{
public:
	ReciprocityEffect(const EffectInfo * pEffectInfo);

	virtual double calculateTieFlipContribution(int alter) const;
	virtual double evaluationStatistic(Network * pNetwork) const;
	virtual double endowmentStatistic(Network * pInitialNetwork,
		Network * pLostTieNetwork) const;
};

}

#endif /*RECIPROCITYEFFECT_H_*/
