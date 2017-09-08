/******************************************************************************
 * SIENA: Simulation Investigation for Empirical Network Analysis
 *
 * Web: http://www.stats.ox.ac.uk/~snijders/siena/
 *
 * File: RecipdegreeActivityEffect.cpp
 *
 * Description: This file contains the implementation of the
 * RecipdegreeActivityEffect class.
 *****************************************************************************/

#include <string>
#include <stdexcept>
#include "RecipdegreeActivityEffect.h"
#include "network/OneModeNetwork.h"
#include "model/variables/NetworkVariable.h"

namespace siena
{

/**
 * Constructor.
 */
RecipdegreeActivityEffect::RecipdegreeActivityEffect(
	const EffectInfo * pEffectInfo) : NetworkEffect(pEffectInfo)
{
}

/**
 * Calculates the contribution of a tie flip to the given actor.
 */
double RecipdegreeActivityEffect::calculateContribution(int alter) const
{
	const OneModeNetwork * pONetwork =
		dynamic_cast<const OneModeNetwork *>(this->pNetwork());

	if (!pONetwork)
	{
		throw runtime_error(
			"One-mode network expected in ReciprocalDegreeBehaviorEffect");
	}
	
	double degree = pONetwork->reciprocalDegree(this->ego());

	if (this->inTieExists(alter))
	{
		degree += this->pNetwork()->outDegree(this->ego());
//		degree += pONetwork->outDegree(this->ego()); should be the same
		if (this->outTieExists(alter))
		{
			degree --;
		}
		else
		{
			degree ++;
		}
	}
	
	double change = degree;
	return change;
}


/**
 * The contribution of the tie from the implicit ego to the given alter
 * to the statistic. It is assumed that preprocessEgo(ego) has been
 * called before.
 */
double RecipdegreeActivityEffect::tieStatistic(int alter)
{
	const OneModeNetwork * pONetwork =
		dynamic_cast<const OneModeNetwork *>(this->pNetwork());

	if (!pONetwork)
	{
		throw runtime_error(
			"One-mode network expected in ReciprocalDegreeBehaviorEffect");
	}
	return pONetwork->reciprocalDegree(this->ego());
}

}
