/******************************************************************************
 * SIENA: Simulation Investigation for Empirical Network Analysis
 * 
 * Web: http://www.stats.ox.ac.uk/~snijders/siena/
 * 
 * File: OutdegreeActivitySqrtEffect.h
 * 
 * Description: This file contains the definition of the
 * OutdegreeActivitySqrtEffect class.
 *****************************************************************************/

#ifndef OUTDEGREEACTIVITYSQRTEFFECT_H_
#define OUTDEGREEACTIVITYSQRTEFFECT_H_

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
 * This class defines the outdegree activity (sqrt) effect defined by
 * s_i(x)= x_{i+}^1.5. The corresponding statistic is sum_i x_{i+}^1.5.
 */
class OutdegreeActivitySqrtEffect : public NetworkEffect
{
public:
	OutdegreeActivitySqrtEffect(const EffectInfo * pEffectInfo);
	
	virtual double calculateTieFlipContribution(int alter) const;
	virtual double evaluationStatistic(Network * pNetwork) const;
	virtual double endowmentStatistic(Network * pInitialNetwork,
		Network * pLostTieNetwork) const;
	
private:
	// Lookup table for fast square root calculations
	SqrtTable * lsqrtTable;
};

}

#endif /*OUTDEGREEACTIVITYSQRTEFFECT_H_*/
