/******************************************************************************
 * SIENA: Simulation Investigation for Empirical Network Analysis
 *
 * Web: http://www.stats.ox.ac.uk/~snijders/siena/
 *
 * File: BehaviorEffect.cpp
 *
 * Description: This file contains the implementation of the class
 * BehaviorEffect.
 *****************************************************************************/

#include <stdexcept>
#include <R.h>

#include "BehaviorEffect.h"
#include "data/Data.h"
#include "data/BehaviorLongitudinalData.h"
#include "model/State.h"
#include "model/EffectInfo.h"
#include "model/variables/BehaviorVariable.h"

namespace siena
{

/**
 * Constructor.
 */
BehaviorEffect::BehaviorEffect(const EffectInfo * pEffectInfo) :
	Effect(pEffectInfo)
{
	this->lpBehaviorData = 0;
	this->lvalues = 0;
}


/**
 * Initializes this effect.
 * @param[in] pData the observed data
 * @param[in] pState the current state of the dependent variables
 * @param[in] period the period of interest
 * @param[in] pCache the cache object to be used to speed up calculations
 */
void BehaviorEffect::initialize(const Data * pData,
	State * pState,
	int period,
	Cache * pCache)
{
	Effect::initialize(pData, pState, period, pCache);

	string name = this->pEffectInfo()->variableName();

	this->lpBehaviorData = pData->pBehaviorData(name);

	if (!this->lpBehaviorData)
	{
		throw logic_error(
			"Data for behavior variable '" + name +"' expected.");
	}

	this->lvalues = pState->behaviorValues(name);
}


/**
 * Returns the number of actors of the behavior variable associated with this
 * effect.
 */
int BehaviorEffect::n() const
{
	return this->lpBehaviorData->n();
}


/**
 * Returns the behavior value of the given actor.
 */
int BehaviorEffect::value(int actor) const
{
	return this->lvalues[actor];
}


/**
 * Returns the behavior value of the given actor centered around the
 * overall mean of all observed values.
 */
double BehaviorEffect::centeredValue(int actor) const
{
	return this->lvalues[actor] - this->lpBehaviorData->overallMean();
}


/**
 * Returns the observed range of the respective behavior variable.
 */
double BehaviorEffect::range() const
{
	return this->lpBehaviorData->range();
}


/**
 * Returns the similarity mean value over all observations.
 */
double BehaviorEffect::similarityMean() const
{
	return this->lpBehaviorData->similarityMean();
}

}
