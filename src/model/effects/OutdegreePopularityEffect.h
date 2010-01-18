/******************************************************************************
 * SIENA: Simulation Investigation for Empirical Network Analysis
 *
 * Web: http://www.stats.ox.ac.uk/~snijders/siena/
 *
 * File: OutdegreePopularityEffect.h
 *
 * Description: This file contains the definition of the
 * OutdegreePopularityEffect class.
 *****************************************************************************/

#ifndef OUTDEGREEPOPULARITYEFFECT_H_
#define OUTDEGREEPOPULARITYEFFECT_H_

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
 * This class defines the indegree popularity effect defined by
 * s_i(x)=sum_j x_{ij} x_{j+} or
 * s_i(x)=sum_j x_{ij} sqrt(x_{j+}), depending on the parameters.
 * The corresponding statistic is
 * the sum of indegree x outdegree products (or indegree x sqrt(outdegree))
 * over all actors.
 */
class OutdegreePopularityEffect : public NetworkEffect
{
public:
	OutdegreePopularityEffect(const EffectInfo * pEffectInfo,
		bool root);

	virtual double calculateContribution(int alter) const;

protected:
	virtual double tieStatistic(int alter);

private:
	// Indicates if the square root of outdegrees must be used
	bool lroot;

	// Lookup table for fast square root calculations
	SqrtTable * lsqrtTable;
};

}

#endif /*OUTDEGREEPOPULARITYEFFECT_H_*/
