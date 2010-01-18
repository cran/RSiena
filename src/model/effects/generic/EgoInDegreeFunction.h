/******************************************************************************
 * SIENA: Simulation Investigation for Empirical Network Analysis
 *
 * Web: http://www.stats.ox.ac.uk/~snijders/siena/
 *
 * File: EgoInDegreeFunction.h
 *
 * Description: This file contains the definition of the
 * EgoInDegreeFunction class.
 *****************************************************************************/

#ifndef EGOINDEGREEFUNCTION_H_
#define EGOINDEGREEFUNCTION_H_

#include "NetworkAlterFunction.h"
#include "IntAlterFunction.h"

namespace siena
{

/**
 * Defines a function that returns the in-degree of the ego regardless
 * of the alter.
 */
class EgoInDegreeFunction: public NetworkAlterFunction, IntAlterFunction
{
public:
	EgoInDegreeFunction(string networkName);

	virtual double value(int alter);
	virtual int intValue(int alter);
};

}

#endif /* EGOINDEGREEFUNCTION_H_ */
