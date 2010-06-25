/******************************************************************************
 * SIENA: Simulation Investigation for Empirical Network Analysis
 *
 * Web: http://www.stats.ox.ac.uk/~snijders/siena/
 *
 * File: NetworkChange.h
 *
 * Description: This file contains the definition of the NetworkChange class.
 *****************************************************************************/

#ifndef NETWORKCHANGE_H_
#define NETWORKCHANGE_H_

#include <vector>

#include "MiniStep.h"

namespace siena
{

// ----------------------------------------------------------------------------
// Section: Forward declarations
// ----------------------------------------------------------------------------

class NetworkLongitudinalData;


// ----------------------------------------------------------------------------
// Section: Class definition
// ----------------------------------------------------------------------------

/**
 * Defines a ministep changing a network variable.
 */
class NetworkChange: public MiniStep
{
public:
	NetworkChange(NetworkLongitudinalData * pData,
		int ego,
		int alter);
	virtual ~NetworkChange();

	virtual bool networkMiniStep() const;
	inline int alter() const;
	virtual void makeChange(DependentVariable * pVariable);
	virtual bool diagonal() const;
	virtual bool missing(int period) const;
	virtual MiniStep * createReverseMiniStep() const;
	virtual MiniStep * createCopyMiniStep() const;
	double evaluationEffectContribution(int alter, int effect) const;
	double endowmentEffectContribution(int alter, int effect) const;
	void allocateEffectContributionArrays(int nEvaluationEffects, 
		int nEndowmentEffects, int m);
	void evaluationEffectContribution(double value, int alter, int effect);
	void endowmentEffectContribution(double value, int alter, int effect);

private:
	// The longitudinal data object for the corresponding network variable
	NetworkLongitudinalData * lpData;

	// The alter whose incoming tie is changed
	int lalter;

	// A vector of change contributions to effects, where one element, 
	// of length the number of effects in the evaluation function, 
	// corresponds to each receiver.

	vector<vector <double> > levaluationEffectContribution;

	// A vector of change contributions to effects, where one element, 
	// of length the number of effects in the endowment function, 
	// corresponds to each receiver.

	vector< vector <double> > lendowmentEffectContribution;
};


// ----------------------------------------------------------------------------
// Section: Inline methods
// ----------------------------------------------------------------------------

/**
 * Returns the alter of this ministep.
 */
int NetworkChange::alter() const
{
	return this->lalter;
}

}

#endif /* NETWORKCHANGE_H_ */
