/******************************************************************************
 * SIENA: Simulation Investigation for Empirical Network Analysis
 *
 * Web: http://www.stats.ox.ac.uk/~snijders/siena/
 *
 * File: SameCovariateTransitiveTripletsEffect.h
 *
 * Description: This file contains the declaration of the class
 * SameCovariateTransitiveTripletsEffect.
 *****************************************************************************/

#ifndef SAMECOVARIATETRANSITIVETRIPLETSEFFECT_H_
#define SAMECOVARIATETRANSITIVETRIPLETSEFFECT_H_

#include "model/effects/NetworkEffect.h"
#include "CovariateDependentNetworkEffect.h"

namespace siena
{

class SameCovariateTransitiveTripletsEffect : 
                   public CovariateDependentNetworkEffect
{
public:
	SameCovariateTransitiveTripletsEffect(
                   const EffectInfo * pEffectInfo, bool reciprocal);
	virtual double calculateContribution(int alter) const;

protected:
	virtual double tieStatistic(int alter);

private:
	// Indicates if the reciprocal version of the 
   // same covariate transitive triplets effect is required;
   // currently not implemented.

	bool lreciprocal;
};

}

#endif /*SAMECOVARIATETRANSITIVETRIPLETSEFFECT_H_*/
