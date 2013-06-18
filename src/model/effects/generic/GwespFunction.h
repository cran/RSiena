/******************************************************************************
 * SIENA: Simulation Investigation for Empirical Network Analysis
 *
 * Web: http://www.stats.ox.ac.uk/~snijders/siena/
 *
 * File: GwespFunction.h
 *
 * Description: This file contains the definition of the
 * GwespFunction class.
 *****************************************************************************/

#ifndef GWESPFUNCTION_H_
#define GWESPFUNCTION_H_
#define MAX_STATISTIC 1000

#include "NetworkAlterFunction.h"

namespace siena
{
// ----------------------------------------------------------------------------
// Section: Forward declarations
// ----------------------------------------------------------------------------

class EgocentricConfigurationTable;
class NetworkCache;
// ----------------------------------------------------------------------------
// Section: Class definition
// ----------------------------------------------------------------------------

/**
 * Defines a function that returns ???.
 */
class GwespFunction: public NetworkAlterFunction
{
public:
	GwespFunction(string networkName,
		EgocentricConfigurationTable * (NetworkCache::*pTable)() const,
		double parameter) ;
	virtual void initialize(const Data * pData,
		State * pState,
		int period,
		Cache * pCache);

	virtual double value(int alter);
private:
	EgocentricConfigurationTable * (NetworkCache::*lpTable)() const;
	double lparameter;
	double lweight;
	double lcumulativeWeight[MAX_STATISTIC];
	EgocentricConfigurationTable *lpInitialisedTable;

};
}

#endif /* GWESPFUNCTION_H_ */
