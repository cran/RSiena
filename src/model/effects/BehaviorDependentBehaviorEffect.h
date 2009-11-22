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

	virtual void initialize(const Data * pData,
		State * pState,
		int period,
		Cache * pCache);

protected:
	inline int interactionValue(int actor) const;
	double centeredInteractionValue(int actor) const;

private:
	// The values of the behavior variable this effect is interacting with
	const int * linteractionValues;

	// The observed data for the behavior variable this effect is
	// interacting with

	BehaviorLongitudinalData * lpInteractionBehaviorData;
};


// ----------------------------------------------------------------------------
// Section: Inline methods
// ----------------------------------------------------------------------------

/**
 * Returns the value of the behavior variable this effect is interacting with.
 */
int BehaviorDependentBehaviorEffect::interactionValue(int actor) const
{
	return this->linteractionValues[actor];
}

}

#endif /*BEHAVIORDEPENDENTBEHAVIOREFFECT_H_*/
