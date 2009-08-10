/******************************************************************************
 * SIENA: Simulation Investigation for Empirical Network Analysis
 * 
 * Web: http://www.stats.ox.ac.uk/~snijders/siena/
 * 
 * File: Effect.cpp
 * 
 * Description: This file contains the implementation of the class Effect.
 *****************************************************************************/

#include "Effect.h"
#include "model/EffectInfo.h"

namespace siena
{

// ----------------------------------------------------------------------------
// Section: Construction, destruction
// ----------------------------------------------------------------------------

/**
 * Constructor.
 */
Effect::Effect(const EffectInfo * pEffectInfo)
{
	this->lpEffectInfo = pEffectInfo;
	this->lweight = pEffectInfo->parameter();
	this->lperiod = 0;
}


/**
 * Destroys this effect.
 */
Effect::~Effect()
{
}


// ----------------------------------------------------------------------------
// Section: Accessors
// ----------------------------------------------------------------------------

/**
 * Stores the weight of this effect to be used when calculating the owner
 * function.
 */
void Effect::weight(double weight)
{
	this->lweight = weight;
}


/**
 * Returns the effect info object this effect is based on.
 */
const EffectInfo * Effect::pEffectInfo() const
{
	return this->lpEffectInfo;
}


// ----------------------------------------------------------------------------
// Section: Initialization
// ----------------------------------------------------------------------------

/**
 * Initializes this effect for calculating contributions of changes to
 * evaluation or endowment functions.
 */
void Effect::initialize(EpochSimulation * pSimulation)
{	
}


/**
 * Initializes this effect before the simulation of the given period.
 */
void Effect::initializeBeforeSimulation(int period)
{
	this->lperiod = period;
}


/**
 * Initializes this effect for calculating the corresponding statistics.
 * @param[in] pData the observed data
 * @param[in] pState the current state of the dependent variables
 * @param[in] period the period of interest
 */
void Effect::initialize(const Data * pData, State * pState, int period)
{
	this->lperiod = period;
}

}
