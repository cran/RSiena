/******************************************************************************
 * SIENA: Simulation Investigation for Empirical Network Analysis
 *
 * Web: http://www.stats.ox.ac.uk/~snijders/siena/
 *
 * File: CovariateDistance2SimilarityEffect.h
 *
 * Description: This file contains the definition of the
 * CovariateDistance2SimilarityEffect class.
 *****************************************************************************/

#ifndef COVARIATEDISTANCE2SIMILARITYEFFECT_H_
#define COVARIATEDISTANCE2SIMILARITYEFFECT_H_

#include "CovariateDistance2NetworkEffect.h"

namespace siena
{

/**
 * Covariate-related distance2 similarity effect (see manual).
 */
class CovariateDistance2SimilarityEffect : public CovariateDistance2NetworkEffect
{
public:
	CovariateDistance2SimilarityEffect(const EffectInfo * pEffectInfo);

	virtual double calculateContribution(int alter) const;

protected:
	virtual double tieStatistic(int alter);

private:
};

}

#endif /*COVARIATDISTANCE2ESIMILARITYEFFECT_H_*/
