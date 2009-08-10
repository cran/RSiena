/******************************************************************************
 * SIENA: Simulation Investigation for Empirical Network Analysis
 *
 * Web: http://www.stats.ox.ac.uk/~snijders/siena/
 *
 * File: BehaviorDependentBehaviorEffect.h
 *
 * Description: This file contains the definition of the
 * BehaviorVariableBehaviorEffect class.
 *****************************************************************************/

#ifndef BEHAVIORDEPENDENTBEHAVIOREFFECT_H_
#define BEHAVIORDEPENDENTBEHAVIOREFFECT_H_

#include "BehaviorEffect.h"

namespace siena
{

// ----------------------------------------------------------------------------
// Section: Forward declarations
// ----------------------------------------------------------------------------

class BehaviorVariable;


// ----------------------------------------------------------------------------
// Section: BehaviorDependentBehaviorEffect class
// ----------------------------------------------------------------------------

/**
 * The base class for all behavior effects depending on some other behavior
 * variable.
 */
class BehaviorDependentBehaviorEffect : public BehaviorEffect
{
public:
	BehaviorDependentBehaviorEffect(const EffectInfo * pEffectInfo);
	virtual ~BehaviorDependentBehaviorEffect();

	virtual void initialize(EpochSimulation * pSimulation);
	virtual void initialize(const Data * pData, State * pState, int period);

protected:
	inline const BehaviorVariable * pBehaviorVariable() const;

private:
	// The behavior variable this effect is interacting with
	const BehaviorVariable * lpBehaviorVariable;
};


// ----------------------------------------------------------------------------
// Section: Inline methods
// ----------------------------------------------------------------------------

/**
 * Returns the behavior variable this effect is interacting with.
 */
const BehaviorVariable * BehaviorDependentBehaviorEffect::pBehaviorVariable()
	const
{
	return this->lpBehaviorVariable;
}

}

#endif /*BEHAVIORDEPENDENTBEHAVIOREFFECT_H_*/
