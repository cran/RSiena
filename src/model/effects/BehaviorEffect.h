/******************************************************************************
 * SIENA: Simulation Investigation for Empirical Network Analysis
 *
 * Web: http://www.stats.ox.ac.uk/~snijders/siena/
 *
 * File: BehaviorEffect.h
 *
 * Description: This file contains the definition of the
 * BehaviorEffect class.
 *****************************************************************************/

#ifndef BEHAVIOREFFECT_H_
#define BEHAVIOREFFECT_H_

#include "Effect.h"

namespace siena
{

// ----------------------------------------------------------------------------
// Section: Forward declarations
// ----------------------------------------------------------------------------

class BehaviorLongitudinalData;


// ----------------------------------------------------------------------------
// Section: BehaviorEffect class
// ----------------------------------------------------------------------------

/**
 * The base class for all behavior effects.
 */
class BehaviorEffect : public Effect
{
public:
	BehaviorEffect(const EffectInfo * pEffectInfo);

	virtual void initialize(const Data * pData,
		State * pState,
		int period,
		Cache * pCache);

	/**
	 * Calculates the change in the statistic corresponding to this effect if
	 * the given actor would change his behavior by the given amount.
	 */
	virtual double calculateChangeContribution(int actor,
		int difference) const = 0;

	/**
	 * Returns the statistic corresponding to this effect as part of
	 * the evaluation function with respect to the behavior variable.
	 */
	virtual double evaluationStatistic(double * currentValues) const = 0;

	/**
	 * Returns the statistic corresponding to this effect as part of
	 * the endowment function with respect to an initial behavior
	 * variable and the current state.
	 */
	virtual double endowmentStatistic(const int * difference,
		double *currentValues) const = 0;

protected:
	int n() const;
	int value(int actor) const;
	double centeredValue(int actor) const;
	double range() const;
	double similarityMean() const;

private:
	BehaviorLongitudinalData * lpBehaviorData;
	const int * lvalues;
};

}

#endif /*BEHAVIOREFFECT_H_*/
