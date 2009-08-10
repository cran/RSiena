/******************************************************************************
 * SIENA: Simulation Investigation for Empirical Network Analysis
 *
 * Web: http://www.stats.ox.ac.uk/~snijders/siena/
 *
 * File: DenseTriadsEffect.h
 *
 * Description: This file contains the declaration of the class
 * DenseTriadsEffect.
 *****************************************************************************/

#ifndef DENSETRIADSEFFECT_H_
#define DENSETRIADSEFFECT_H_

#include "NetworkEffect.h"

namespace siena
{

/**
 * This class defines the dense triads effect. See the manual for effect
 * definitions.
 */
class DenseTriadsEffect : public NetworkEffect
{
public:
	DenseTriadsEffect(const EffectInfo * pEffectInfo);

	virtual double calculateTieFlipContribution(int alter) const;
	virtual bool usesTable(const ConfigurationTable * pTable) const;

protected:
	virtual double statistic(Network * pNetwork,
		Network * pSummationTieNetwork) const;

private:
	int ldensity;
};

}

#endif /*DENSETRIADSEFFECT_H_*/
