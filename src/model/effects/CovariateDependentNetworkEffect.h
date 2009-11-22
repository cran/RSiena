/******************************************************************************
 * SIENA: Simulation Investigation for Empirical Network Analysis
 *
 * Web: http://www.stats.ox.ac.uk/~snijders/siena/
 *
 * File: CovariateDependentNetworkEffect.h
 *
 * Description: This file contains the definition of the
 * CovariateDependentNetworkEffect class.
 *****************************************************************************/

#ifndef COVARIATEDEPENDENTNETWORKEFFECT_H_
#define COVARIATEDEPENDENTNETWORKEFFECT_H_

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
 * covariate (constant, changing, or dependent behavior variable).
 */
class CovariateDependentNetworkEffect : public NetworkEffect
{
public:
	CovariateDependentNetworkEffect(const EffectInfo * pEffectInfo);

	virtual void initialize(const Data * pData,
		State * pState,
		int period,
		Cache * pCache);

protected:
	double value(int i) const;
	bool missing(int i) const;
	double similarity(int i, int j) const;

private:
	ConstantCovariate * lpConstantCovariate;
	ChangingCovariate * lpChangingCovariate;
	BehaviorLongitudinalData * lpBehaviorData;

	// The current value of a behavior variable per each actor.
	// This array is 0 for covariate-based effects.

	const int * lvalues;
};

}

#endif /*COVARIATEDEPENDENTNETWORKEFFECT_H_*/
