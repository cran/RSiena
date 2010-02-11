/******************************************************************************
 * SIENA: Simulation Investigation for Empirical Network Analysis
 *
 * Web: http://www.stats.ox.ac.uk/~snijders/siena/
 *
 * File: BehaviorChange.cpp
 *
 * Description: This file contains the implementation of the class
 * BehaviorChange.
 *****************************************************************************/

#include "BehaviorChange.h"
#include "model/variables/BehaviorVariable.h"

namespace siena
{

/**
 * Constructs a new behavior ministep.
 * @param[in] ego the actor making the change
 * @param[in] variableName the name of the dependent variable to be changed
 * @param[in] difference the amount of change
 * (-1,0,+1 for dichotomous variables)
 */
BehaviorChange::BehaviorChange(int ego,
	string variableName,
	int difference) : MiniStep(ego, variableName, difference)
{
}


/**
 * Deallocates this ministep.
 */
BehaviorChange::~BehaviorChange()
{
}


/**
 * Changes the given behavior variable according to this ministep.
 */
void BehaviorChange::makeChange(DependentVariable * pVariable)
{
	MiniStep::makeChange(pVariable);

	if (this->difference() != 0)
	{
		BehaviorVariable * pBehaviorVariable =
			dynamic_cast<BehaviorVariable *>(pVariable);
		int oldValue = pBehaviorVariable->value(this->ego());
		pBehaviorVariable->value(this->ego(), oldValue + this->difference());
	}
}

}
