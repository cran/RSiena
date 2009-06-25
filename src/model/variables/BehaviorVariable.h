/******************************************************************************
 * SIENA: Simulation Investigation for Empirical Network Analysis
 *
 * Web: http://www.stats.ox.ac.uk/~snijders/siena/
 *
 * File: BehaviorVariable.h
 *
 * Description: This file contains the definition of the
 * BehaviorVariable class.
 *****************************************************************************/

#ifndef BEHAVIORVARIABLE_H_
#define BEHAVIORVARIABLE_H_

#include "DependentVariable.h"

namespace siena
{

// ----------------------------------------------------------------------------
// Section: Forward declarations
// ----------------------------------------------------------------------------

class BehaviorLongitudinalData;


// ----------------------------------------------------------------------------
// Section: BehaviorVariable class
// ----------------------------------------------------------------------------

/**
 * This class represents the state of a behavioral dependent variable.
 * @see DependentVariable
 */
class BehaviorVariable : public DependentVariable
{
public:
	BehaviorVariable(BehaviorLongitudinalData * pData,
		EpochSimulation * pSimulation);
	virtual ~BehaviorVariable();

	virtual int m() const;
	virtual LongitudinalData * pData() const;
	virtual void initialize(int period);
	virtual void actOnJoiner(const ActorSet * pActorSet, int actor);
	virtual void actOnLeaver(const ActorSet * pActorSet, int actor);
	virtual void setLeaverBack(const ActorSet * pActorSet, int actor);

	virtual void makeChange(int actor);

	bool missingStartValue(int actor) const;
	int value(int actor) const;
	double centeredValue(int actor) const;
	double similarity(int i, int j) const;
	int predictorValue(int actor) const;
	void predictorValue(int actor, int value);
	double centeredPredictorValue(int actor) const;
	const int * values() const;
	int range() const;
	double similarityMean() const;

private:
	double totalEvaluationContribution(int actor,
		int difference) const;
	double totalEndowmentContribution(int actor,
		int difference) const;
	void accumulateScores(int difference) const;

	// The observed data for this behavioral variable
	BehaviorLongitudinalData * lpData;

	// The underlying set of actors.
	const ActorSet * lpActorSet;

	// The current value of the variable per each actor
	int * lvalues;

	// The mean value of the variable over all active actors
	double lmean;

	// A two-dimensional array of change contributions to effects, where
	// rows correspond to differences and columns correspond to effects in the
	// evaluation function.

	double ** levaluationEffectContribution;

	// A two-dimensional array of change contributions to effects, where
	// rows correspond to differences and columns correspond to effects in the
	// endowment function.

	double ** lendowmentEffectContribution;

	// Selection probability per each difference

	double * lprobabilities;

	// Predictor values of the variable: starting values with missings zeroed.

	int * lpredictorValues;
};

}

#endif /*BEHAVIORVARIABLE_H_*/
