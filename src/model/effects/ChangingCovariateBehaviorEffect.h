/******************************************************************************
 * SIENA: Simulation Investigation for Empirical Network Analysis
 *
 * Web: http://www.stats.ox.ac.uk/~snijders/siena/
 *
 * File: ChangingCovariateBehaviorEffect.h
 *
 * Description: This file contains the definition of the
 * ChangingCovariateBehaviorEffect class.
 *****************************************************************************/

#ifndef CHANGINGCOVARIATEBEHAVIOREFFECT_H_
#define CHANGINGCOVARIATEBEHAVIOREFFECT_H_

#include "BehaviorEffect.h"

namespace siena
{

// ----------------------------------------------------------------------------
// Section: Forward declarations
// ----------------------------------------------------------------------------

class ChangingCovariate;


// ----------------------------------------------------------------------------
// Section: ChangingCovariateBehaviorEffect class
// ----------------------------------------------------------------------------

/**
 * The base class for all behavior effects depending on some changing covariate.
 */
class ChangingCovariateBehaviorEffect : public BehaviorEffect
{
public:
	ChangingCovariateBehaviorEffect(const EffectInfo * pEffectInfo);

	virtual void initialize(const Data * pData,
		State * pState,
		int period,
		Cache * pCache);

protected:
	inline const ChangingCovariate * pCovariate() const;

private:
	// The covariate this effect is interacting with
	const ChangingCovariate * lpCovariate;
};


// ----------------------------------------------------------------------------
// Section: Inline methods
// ----------------------------------------------------------------------------

/**
 * Returns the covariate this effect is interacting with.
 */
const ChangingCovariate * ChangingCovariateBehaviorEffect::pCovariate() const
{
	return this->lpCovariate;
}

}

#endif /*CHANGINGCOVARIATEBEHAVIOREFFECT_H_*/
