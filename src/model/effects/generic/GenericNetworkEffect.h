/******************************************************************************
 * SIENA: Simulation Investigation for Empirical Network Analysis
 *
 * Web: http://www.stats.ox.ac.uk/~snijders/siena/
 *
 * File: GenericNetworkEffect.h
 *
 * Description: This file contains the definition of the
 * GenericNetworkEffect class.
 *****************************************************************************/


#ifndef GENERICNETWORKEFFECT_H_
#define GENERICNETWORKEFFECT_H_

#include "model/effects/NetworkEffect.h"

namespace siena
{

class AlterFunction;


class GenericNetworkEffect: public NetworkEffect
{
public:
	GenericNetworkEffect(const EffectInfo * pEffectInfo,
		AlterFunction * pFunction);
	GenericNetworkEffect(const EffectInfo * pEffectInfo,
		AlterFunction * pEffectFunction,
		AlterFunction * pStatisticFunction);
	virtual ~GenericNetworkEffect();

	virtual void initialize(const Data * pData,
		State * pState,
		int period,
		Cache * pCache);
	virtual void preprocessEgo(int ego);

	virtual double calculateContribution(int alter) const;

protected:
	virtual double tieStatistic(int alter);

private:
	AlterFunction * lpEffectFunction;
	AlterFunction * lpStatisticFunction;
};

}

#endif /* GENERICNETWORKEFFECT_H_ */
