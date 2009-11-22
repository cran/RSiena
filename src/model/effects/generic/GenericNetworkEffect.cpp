/******************************************************************************
 * SIENA: Simulation Investigation for Empirical Network Analysis
 *
 * Web: http://www.stats.ox.ac.uk/~snijders/siena/
 *
 * File: GenericNetworkEffect.cpp
 *
 * Description: This file contains the implementation of the class
 * GenericNetworkEffect.
 *****************************************************************************/

#include "GenericNetworkEffect.h"
#include "network/Network.h"
#include "network/IncidentTieIterator.h"
#include "model/effects/generic/AlterFunction.h"
#include "model/tables/Cache.h"

namespace siena
{

GenericNetworkEffect::GenericNetworkEffect(const EffectInfo * pEffectInfo,
	AlterFunction * pFunction) : NetworkEffect(pEffectInfo)
{
	this->lpFunction = pFunction;
}


GenericNetworkEffect::~GenericNetworkEffect()
{
}


/**
 * Initializes this effect.
 * @param[in] pData the observed data
 * @param[in] pState the current state of the dependent variables
 * @param[in] period the period of interest
 * @param[in] pCache the cache object to be used to speed up calculations
 */
void GenericNetworkEffect::initialize(const Data * pData,
	State * pState,
	int period,
	Cache * pCache)
{
	NetworkEffect::initialize(pData, pState, period, pCache);
	this->lpFunction->initialize(pData, pState, period, pCache);
}


/**
 * Does the necessary preprocessing work for calculating the tie flip
 * contributions for a specific ego. This method must be invoked before
 * calling NetworkEffect::calculateTieFlipContribution(...).
 */
void GenericNetworkEffect::preprocessEgo(int ego)
{
	NetworkEffect::preprocessEgo(ego);
	this->lpFunction->preprocessEgo(ego);
}


/**
 * Assuming that the ego would flip the tie to the given actor,
 * this method calculates the change in the statistic corresponding
 * to this effect. The method has to be overriden by all concrete
 * effect classes.
 */
double GenericNetworkEffect::calculateContribution(int alter) const
{
	return this->lpFunction->value(alter);
}


/**
 * See base class.
 */
double GenericNetworkEffect::statistic(const Network * pSummationTieNetwork)
	const
{
	int n = pSummationTieNetwork->n();
	Cache * pCache = this->pCache();
	double statistic = 0;

	for (int i = 0; i < n; i++)
	{
		pCache->initialize(i);
		this->lpFunction->preprocessEgo(i);

		for (IncidentTieIterator iter = pSummationTieNetwork->outTies(i);
			iter.valid();
			iter.next())
		{
			statistic += this->lpFunction->value(iter.actor());
		}
	}

	return statistic;
}

}
