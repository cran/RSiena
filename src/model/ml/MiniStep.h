/******************************************************************************
 * SIENA: Simulation Investigation for Empirical Network Analysis
 *
 * Web: http://www.stats.ox.ac.uk/~snijders/siena/
 *
 * File: MiniStep.h
 *
 * Description: This file contains the definition of the MiniStep class.
 *****************************************************************************/


#ifndef MINISTEP_H_
#define MINISTEP_H_

#include <string>

using namespace std;

namespace siena
{

// ----------------------------------------------------------------------------
// Section: Forward declarations
// ----------------------------------------------------------------------------

class DependentVariable;


// ----------------------------------------------------------------------------
// Section: Class definition
// ----------------------------------------------------------------------------

/**
 * Defines a single ministep as part of a Chain object in the
 * Maximum-Likelihood calculations.
 */
class MiniStep
{
public:
	MiniStep(int ego, string variableName, int difference);
	virtual ~MiniStep();

	inline int ego() const;
	string variableName() const;
	inline int difference() const;

	inline double logProbability() const;
	void logProbability(double probability);

	inline double reciprocalRate() const;
	void reciprocalRate(double value);

	inline MiniStep * pPrevious() const;
	inline MiniStep * pNext() const;
	void pPrevious(MiniStep * pMiniStep);
	void pNext(MiniStep * pMiniStep);

	virtual void makeChange(DependentVariable * pVariable);

private:
	// The actor making the change
	int lego;

	// The name of the dependent variable to be changed
	string lvariableName;

	// The amount of change (+1,0,-1 for dichotomous variables)
	int ldifference;

	// Log probability of making this ministep
	double llogProbability;

	// Reciprocal of aggregate (summed) rate function immediately
	// before this ministep.

	double lreciprocalRate;

	// Points to the previous ministep in the same chain
	MiniStep * lpPrevious;

	// Points to the next ministep in the same chain
	MiniStep * lpNext;
};


// ----------------------------------------------------------------------------
// Section: Inline methods
// ----------------------------------------------------------------------------

/**
 * Returns the ego of this ministep.
 */
int MiniStep::ego() const
{
	return this->lego;
}


/**
 * Returns the amount of change in this ministep
 * (+1,0,-1 for dichotomous variables).
 */
int MiniStep::difference() const
{
	return this->ldifference;
}


/**
 * Returns the log probability of making this ministep.
 */
double MiniStep::logProbability() const
{
	return this->llogProbability;
}


/**
 * Returns the reciprocal of aggregate rate function immediately before
 * this ministep.
 */
double MiniStep::reciprocalRate() const
{
	return this->lreciprocalRate;
}


/**
 * Returns the previous ministep in the same chain.
 */
MiniStep * MiniStep::pPrevious() const
{
	return this->lpPrevious;
}


/**
 * Returns the next ministep in the same chain.
 */
MiniStep * MiniStep::pNext() const
{
	return this->lpNext;
}

}

#endif /* MINISTEP_H_ */
