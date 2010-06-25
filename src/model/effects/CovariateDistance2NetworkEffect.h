/******************************************************************************
 * SIENA: Simulation Investigation for Empirical Network Analysis
 *
 * Web: http://www.stats.ox.ac.uk/~snijders/siena/
 *
 * File: CovariateDistance2NetworkEffect.h
 *
 * Description: This file contains the definition of the
 * CovariateDistance2NetworkEffect class.
 *****************************************************************************/

#ifndef COVARIATEDISTANCE2NETWORKEFFECT_H_
#define COVARIATEDISTANCE2NETWORKEFFECT_H_
#include <string>
#include "CovariateDependentNetworkEffect.h"
#include "NetworkEffect.h"

namespace siena
{

// ----------------------------------------------------------------------------
// Section: Forward declarations
// ----------------------------------------------------------------------------

class ConstantCovariate;
class ChangingCovariate;
class BehaviorVariable;
class BehaviorLongitudinalData;


// ----------------------------------------------------------------------------
// Section: Class definition
// ----------------------------------------------------------------------------

/**
 * The base class for network effects depending on an individual
 * covariate (constant, changing, or dependent behavior variable) at distance 2.
 */
class CovariateDistance2NetworkEffect : public CovariateDependentNetworkEffect
{
public:
	CovariateDistance2NetworkEffect(const EffectInfo * pEffectInfo);
	virtual ~CovariateDistance2NetworkEffect();

	virtual void initialize(const Data * pData,
		State * pState,
		int period,
		Cache * pCache);

	virtual void preprocessEgo(int ego);


protected:
	bool missingDummy(int i) const;
	double averageAlterValue(int i) const;
	double similarityNetwork(int i, int j, std::string networkName) const;

private:
	double * laverageAlterValues;
	bool * laverageAlterMissing;
};

}

#endif /*COVARIATEDISTANCE2NETWORKEFFECT_H_*/
