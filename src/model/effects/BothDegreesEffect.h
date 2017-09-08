/******************************************************************************
 * SIENA: Simulation Investigation for Empirical Network Analysis
 *
 * Web: http://www.stats.ox.ac.uk/~snijders/siena/
 *
 * File: BothDegreesEffect.h
 *
 * Description: This file contains the definition of the
 * BothDegreesEffect class.
 *****************************************************************************/

#ifndef BOTHDEGREESEFFECT_H_
#define BOTHDEGREESEFFECT_H_

#include "NetworkEffect.h"

namespace siena
{

// ----------------------------------------------------------------------------
// Section: Forward declarations
// ----------------------------------------------------------------------------

class SqrtTable;

// ----------------------------------------------------------------------------
// Section: Class definition
// ----------------------------------------------------------------------------

/**
 * This class defines the both degrees effect defined as
 * the sum of the indegree popularity and outdegree activity effects.
 */
class BothDegreesEffect : public NetworkEffect
{
public:
	BothDegreesEffect(const EffectInfo * pEffectInfo);
	virtual double calculateContribution(int alter) const;

protected:
	virtual double tieStatistic(int alter);

private:
	// Indicates if the square root of degrees must be used
	bool lroot;
	// Lookup table for fast square root calculations
	SqrtTable * lsqrtTable;
};

}

#endif /*BOTHDEGREESEFFECT_H_*/
