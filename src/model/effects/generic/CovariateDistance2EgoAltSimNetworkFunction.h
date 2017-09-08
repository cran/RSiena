/******************************************************************************
 * SIENA: Simulation Investigation for Empirical Network Analysis
 *
 * Web: http://www.stats.ox.ac.uk/~snijders/siena/
 *
 * File: CovariateDistance2EgoAltSimNetworkFunction.h
 *
 * Description: This file contains the definition of the
 * CovariateDistance2EgoAltSimNetworkFunction class.
 *****************************************************************************/

#ifndef COVARIATEDISTANCE2EGOALTSIMNETWORKFUNCTION_H_
#define COVARIATEDISTANCE2EGOALTSIMNETWORKFUNCTION_H_

#include "CovariateDistance2NetworkFunction.h"

namespace siena
{

class CovariateDistance2EgoAltSimNetworkFunction: public
	CovariateDistance2NetworkFunction
{
public:
	CovariateDistance2EgoAltSimNetworkFunction(string networkName,
		string covariateName, bool excludeMissing, bool incoming);
	virtual double value(int alter);

private:
	bool lexcludeMissing;
	bool lincoming;
};

}

#endif /* COVARIATEDISTANCE2INSIMNETWORKFUNCTION_H_ */
