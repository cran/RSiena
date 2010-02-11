/******************************************************************************
 * SIENA: Simulation Investigation for Empirical Network Analysis
 *
 * Web: http://www.stats.ox.ac.uk/~snijders/siena/
 *
 * File: BehaviorChange.h
 *
 * Description: This file contains the definition of the BehaviorChange class.
 *****************************************************************************/


#ifndef BEHAVIORCHANGE_H_
#define BEHAVIORCHANGE_H_

#include "MiniStep.h"

namespace siena
{

/**
 * Defines a ministep changing a behavior variable.
 */
class BehaviorChange: public MiniStep
{
public:
	BehaviorChange(int ego,
		string variableName,
		int difference);
	virtual ~BehaviorChange();

	virtual void makeChange(DependentVariable * pVariable);
};

}

#endif /* BEHAVIORCHANGE_H_ */
