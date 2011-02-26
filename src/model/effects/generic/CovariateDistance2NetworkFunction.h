/******************************************************************************
 * SIENA: Simulation Investigation for Empirical Network Analysis
 *
 * Web: http://www.stats.ox.ac.uk/~snijders/siena/
 *
 * File: CovariateDistance2NetworkFunction.h
 *
 * Description: This file contains the definition of the
 * CovariateDistance2NetworkFunction class.
 *****************************************************************************/

#ifndef COVARIATEDISTANCE2NETWORKFUNCTION_H_
#define COVARIATEDISTANCE2NETWORKFUNCTION_H_

#include "CovariateNetworkAlterFunction.h"

using namespace std;

namespace siena
{

class CovariateDistance2NetworkFunction: public CovariateNetworkAlterFunction
{
public:
	CovariateDistance2NetworkFunction(string networkName, string covariateName);
	virtual ~CovariateDistance2NetworkFunction();
	virtual void initialize(const Data * pData,
		State * pState,
		int period,
		Cache * pCache);
	virtual void preprocessEgo(int ego);


protected:
	bool missingDummy(int i) const;
	double averageAlterValue(int i) const;
	double similarityNetwork(int i, int j) const;

private:
	double * laverageAlterValues;
	bool * laverageAlterMissing;

};



}

#endif /* COVARIATEDISTANCE2NETWORKFUNCTION_H_ */
