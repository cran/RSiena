/******************************************************************************
 * SIENA: Simulation Investigation for Empirical Network Analysis
 *
 * Web: http://www.stats.ox.ac.uk/~snijders/siena/
 *
 * File: BehaviorChange.h
 *
 * Description: This file contains the definition of the BehaviorChange class.
 *****************************************************************************/


#ifndef BEHAVIORCHANGE_H_
#define BEHAVIORCHANGE_H_

#include <vector>

#include "MiniStep.h"

namespace siena
{

// ----------------------------------------------------------------------------
// Section: Forward declarations
// ----------------------------------------------------------------------------

class BehaviorLongitudinalData;


// ----------------------------------------------------------------------------
// Section: Class definition
// ----------------------------------------------------------------------------

/**
 * Defines a ministep changing a behavior variable.
 */
class BehaviorChange: public MiniStep
{
public:
	BehaviorChange(BehaviorLongitudinalData * pData,
		int ego,
		int difference);
	virtual ~BehaviorChange();

	inline int difference() const;

	virtual bool behaviorMiniStep() const;
	virtual void makeChange(DependentVariable * pVariable);
	virtual bool missing(int period) const;
	virtual bool missingStart(int period) const;
	virtual bool missingEnd(int period) const;
	virtual MiniStep * createReverseMiniStep() const;
	virtual MiniStep * createCopyMiniStep() const;
	virtual bool firstOfConsecutiveCancelingPair() const;
	double evaluationEffectContribution(int difference, int effect) const;
	double endowmentEffectContribution(int difference, int effect) const;
	double creationEffectContribution(int difference, int effect) const;
	void allocateEffectContributionArrays(int nEvaluationEffects,
		int nEndowmentEffects,
		int nCreationEffects);
	void evaluationEffectContribution(double value, int difference, int effect);
	void endowmentEffectContribution(double value, int difference, int effect);
	void creationEffectContribution(double value, int difference, int effect);

private:
	// The longitudinal data object for the corresponding behavior variable
	BehaviorLongitudinalData * lpData;

	// The amount of change
	int ldifference;

	// A vector of change contributions to effects, where one element,
	// of length the number of effects in the evaluation function,
	// corresponds to each possible difference.

	vector<vector <double> > levaluationEffectContribution;

	// A vector of change contributions to effects, where one element,
	// of length the number of effects in the endowment function,
	// corresponds to each possible difference.

	vector<vector <double> > lendowmentEffectContribution;

	// A vector of change contributions to effects, where one element,
	// of length the number of effects in the tie creation function,
	// corresponds to each possible difference.

	vector<vector <double> > lcreationEffectContribution;
};


// ----------------------------------------------------------------------------
// Section: Inline methods
// ----------------------------------------------------------------------------

/**
 * Returns the amount of change in this ministep.
 */
int BehaviorChange::difference() const
{
	return this->ldifference;
}


}

#endif /* BEHAVIORCHANGE_H_ */
