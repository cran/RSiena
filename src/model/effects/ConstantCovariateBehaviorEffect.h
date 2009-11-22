/******************************************************************************
 * SIENA: Simulation Investigation for Empirical Network Analysis
 *
 * Web: http://www.stats.ox.ac.uk/~snijders/siena/
 *
 * File: ConstantCovariateBehaviorEffect.h
 *
 * Description: This file contains the definition of the
 * ConstantCovariateBehaviorEffect class.
 *****************************************************************************/

#ifndef CONSTANTCOVARIATEBEHAVIOREFFECT_H_
#define CONSTANTCOVARIATEBEHAVIOREFFECT_H_

#include "BehaviorEffect.h"

namespace siena
{

// ----------------------------------------------------------------------------
// Section: Forward declarations
// ----------------------------------------------------------------------------

class ConstantCovariate;


// ----------------------------------------------------------------------------
// Section: ConstantCovariateBehaviorEffect class
// ----------------------------------------------------------------------------

/**
 * The base class for all behavior effects depending on some covariate.
 */
class ConstantCovariateBehaviorEffect : public BehaviorEffect
{
public:
	ConstantCovariateBehaviorEffect(const EffectInfo * pEffectInfo);
	virtual ~ConstantCovariateBehaviorEffect();

	virtual void initialize(const Data * pData,
		State * pState,
		int period,
		Cache * pCache);

protected:
	inline const ConstantCovariate * pCovariate() const;

private:
	// The covariate this effect is interacting with
	const ConstantCovariate * lpCovariate;
};


// ----------------------------------------------------------------------------
// Section: Inline methods
// ----------------------------------------------------------------------------

/**
 * Returns the covariate this effect is interacting with.
 */
const ConstantCovariate * ConstantCovariateBehaviorEffect::pCovariate()
	const
{
	return this->lpCovariate;
}

}

#endif /*CONSTANTCOVARIATEBEHAVIOREFFECT_H_*/
