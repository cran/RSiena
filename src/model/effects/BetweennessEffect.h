/******************************************************************************
 * SIENA: Simulation Investigation for Empirical Network Analysis
 *
 * Web: http://www.stats.ox.ac.uk/~snijders/siena/
 *
 * File: BetweennessEffect.h
 *
 * Description: This file contains the declaration of the class
 * BetweennessEffect.
 *****************************************************************************/

#ifndef BETWEENNESSEFFECT_H_
#define BETWEENNESSEFFECT_H_

#include "NetworkEffect.h"

namespace siena
{

/**
 * This class defines the betweenness effect.
 */
class BetweennessEffect : public NetworkEffect
{
public:
	BetweennessEffect(const EffectInfo * pEffectInfo);

	virtual double calculateTieFlipContribution(int alter) const;
	virtual bool usesTable(const ConfigurationTable * pTable) const;

protected:
	virtual double statistic(Network * pNetwork,
		Network * pSummationTieNetwork) const;
};

}

#endif /*BETWEENNESSEFFECT_H_*/
