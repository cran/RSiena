/******************************************************************************
 * SIENA: Simulation Investigation for Empirical Network Analysis
 *
 * Web: http://www.stats.ox.ac.uk/~snijders/siena/
 *
 * File: MiniStep.cpp
 *
 * Description: This file contains the implementation of the class MiniStep.
 *****************************************************************************/

#include "MiniStep.h"
#include "model/variables/DependentVariable.h"

namespace siena
{

/**
 * Constructs a new ministep.
 * @param[in] ego the actor making the change
 * @param[in] variableName the name of the dependent variable to be changed
 * @param[in] difference the amount of change
 * (-1,0,+1 for dichotomous variables)
 */
MiniStep::MiniStep(int ego, string variableName, int difference)
{
	this->lego = ego;
	this->lvariableName = variableName;
	this->ldifference = difference;
	this->llogProbability = 0;
	this->lreciprocalRate = 0;
	this->lpPrevious = 0;
	this->lpNext = 0;
}


/**
 * Deallocates this ministep.
 */
MiniStep::~MiniStep()
{
}


/**
 * Returns the name of the dependent variable that this ministep is changing.
 */
string MiniStep::variableName() const
{
	return this->lvariableName;
}


/**
 * Stores the log probability of making this ministep.
 */
void MiniStep::logProbability(double probability)
{
	this->llogProbability = probability;
}


/**
 * Stores the reciprocal of aggregate rate function immediately before
 * this ministep.
 */
void MiniStep::reciprocalRate(double value)
{
	this->lreciprocalRate = value;
}


/**
 * Stores the pointer to the previous ministep in the same chain.
 */
void MiniStep::pPrevious(MiniStep * pMiniStep)
{
	this->lpPrevious = pMiniStep;
}


/**
 * Stores the pointer to the next ministep in the same chain.
 */
void MiniStep::pNext(MiniStep * pMiniStep)
{
	this->lpNext = pMiniStep;
}


/**
 * Changes the given dependent variable according to this ministep.
 */
void MiniStep::makeChange(DependentVariable * pVariable)
{
	// Nothing in the base class.
}

}
