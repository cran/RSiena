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

	virtual double calculateContribution(int alter) const;

protected:
	virtual double tieStatistic(int alter);
};

}

#endif /*TRANSITIVETRIPLETSEFFECT_H_*/
