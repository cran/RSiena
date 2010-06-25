/******************************************************************************
 * SIENA: Simulation Investigation for Empirical Network Analysis
 *
 * Web: http://www.stats.ox.ac.uk/~snijders/siena/
 *
 * File: CovariateDistance2AlterEffect.h
 *
 * Description: This file contains the definition of the
 * CovariateDistance2AlterEffect class.
 *****************************************************************************/

#ifndef COVARIATEDISTANCE2ALTEREFFECT_H_
#define COVARIATEDISTANCE2ALTEREFFECT_H_

#include "CovariateDistance2NetworkEffect.h"

namespace siena
{

/**
 * Covariate - alter at distance 2 effect (see manual).
 */
class CovariateDistance2AlterEffect : public CovariateDistance2NetworkEffect
{
public:
	CovariateDistance2AlterEffect(const EffectInfo * pEffectInfo);

	virtual double calculateContribution(int alter) const;

protected:
	virtual double tieStatistic(int alter);

private:

};

}

#endif /*COVARIATEDISTANCE2ALTEREFFECT_H_*/
