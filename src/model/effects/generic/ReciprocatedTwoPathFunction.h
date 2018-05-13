/******************************************************************************
 * SIENA: Simulation Investigation for Empirical Network Analysis
 *
 * Web: http://www.stats.ox.ac.uk/~snijders/siena/
 *
 * File: ReciprocatedTwoPathFunction.h
 *
 * Description: This file contains the definition of the
 * ReciprocatedTwoPathFunction class.
 *****************************************************************************/

#ifndef RECIPROCATEDTWOPATHFUNCTION_H_
#define RECIPROCATEDTWOPATHFUNCTION_H_

#include "OneModeNetworkAlterFunction.h"
#include "IntAlterFunction.h"

namespace siena
{

class ConfigurationTable;


/**
 * Defines a function that returns the number of reciprocated two-paths between
 * the ego and alters in a network of the given name. A two-path i -> j -> h is
 * reciprocated if there exists the two-path h -> j -> i.
 */
class ReciprocatedTwoPathFunction:
	public OneModeNetworkAlterFunction, IntAlterFunction
{
public:
	ReciprocatedTwoPathFunction(std::string networkName);

	virtual void initialize(const Data * pData,
		State * pState,
		int period,
		Cache * pCache);

	virtual double value(int alter);
	virtual int intValue(int alter);

private:
	ConfigurationTable * lpTable;
};

}

#endif /* RECIPROCATEDTWOPATHFUNCTION_H_ */
