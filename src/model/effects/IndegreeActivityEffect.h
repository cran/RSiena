/******************************************************************************
 * SIENA: Simulation Investigation for Empirical Network Analysis
 *
 * Web: http://www.stats.ox.ac.uk/~snijders/siena/
 *
 * File: IndegreeActivityEffect.h
 *
 * Description: This file contains the definition of the
 * IndegreeActivityEffect class.
 *****************************************************************************/

#ifndef INDEGREEACTIVITYEFFECT_H_
#define INDEGREEACTIVITYEFFECT_H_

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
 * This class defines the indegree activity effect defined by
 * s_i(x)= x_{+i} x_{i+} or
 * s_i(x)= sqrt(x_{+i}) x_{i+}, depending on the parameter of the
 * constructor.
 * The corresponding statistic is
 * the sum of indegree x outdegree products (or sqrt(indegree) x outdegree)
 * over all actors.
 */
class IndegreeActivityEffect : public NetworkEffect
{
public:
	IndegreeActivityEffect(const EffectInfo * pEffectInfo, bool root);

	virtual double calculateContribution(int alter) const;

protected:
	virtual double statistic(const Network * pSummationTieNetwork) const;

private:
	// Indicates if the square root of indegrees must be used
	bool lroot;

	// Lookup table for fast square root calculations
	SqrtTable * lsqrtTable;
};

}

#endif /*INDEGREEACTIVITYEFFECT_H_*/
