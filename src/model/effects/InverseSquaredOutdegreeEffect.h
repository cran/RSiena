/******************************************************************************
 * SIENA: Simulation Investigation for Empirical Network Analysis
 *
 * Web: http://www.stats.ox.ac.uk/~snijders/siena/
 *
 * File: InverseSquaredOutdegreeEffect.h
 *
 * Description: This file contains the declaration of the class
 * InverseSquaredOutdegreeEffect.
 *****************************************************************************/

#ifndef INVERSESQUAREDOUTDEGREEEFFECT_H_
#define INVERSESQUAREDOUTDEGREEEFFECT_H_

#include "NetworkEffect.h"

namespace siena
{

/**
 * This class defines the inverse squared outdegree effect defined as
 * s_i = 1/[(outdegree(i) + c) * (outdegree(i) + c + 1)], where c is
 * a parameter. See the manual for effect definitions.
 */
class InverseSquaredOutdegreeEffect : public NetworkEffect
{
public:
	InverseSquaredOutdegreeEffect(const EffectInfo * pEffectInfo);

	virtual double calculateContribution(int alter) const;
	virtual double evaluationStatistic() const;
	virtual double endowmentStatistic(Network * pLostTieNetwork) const;

private:
	double lc;
};

}

#endif /*INVERSESQUAREDOUTDEGREEEFFECT_H_*/
