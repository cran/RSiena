/******************************************************************************
 * SIENA: Simulation Investigation for Empirical Network Analysis
 *
 * Web: http://www.stats.ox.ac.uk/~snijders/siena/
 *
 * File: CovariateMixedNetworkAlterFunction.h
 *
 * Description: This file contains the definition of the
 * CovariateMixedNetworkAlterFunction class.
 *****************************************************************************/

// combines MixedNetworkAlterFunction and CovariateNetworkAlterFunction

#ifndef COVARIATEMIXEDNETWORKALTERFUNCTION_H_
#define COVARIATEMIXEDNETWORKALTERFUNCTION_H_

#include <string>
#include "AlterFunction.h"
#include "utils/NamedObject.h"
#include "MixedNetworkAlterFunction.h"

namespace siena
{

// ----------------------------------------------------------------------------
// Section: Forward declarations
// ----------------------------------------------------------------------------

class Network;
class TwoNetworkCache;
class NetworkCache;
class ConstantCovariate;
class ChangingCovariate;
class BehaviorVariable;
class BehaviorLongitudinalData;

class CovariateMixedNetworkAlterFunction: public MixedNetworkAlterFunction
{
public:
	CovariateMixedNetworkAlterFunction(string firstNetworkName,
		string secondNetworkName, string covariateName);
	virtual ~CovariateMixedNetworkAlterFunction();

	virtual void initialize(const Data * pData,
		State * pState, int period, Cache * pCache);

protected:
	double value(int i) const;
	bool missing(int i) const;
	double similarity(int i, int j) const;
	ConstantCovariate * pConstantCovariate() const;
	ChangingCovariate * pChangingCovariate() const;
	BehaviorLongitudinalData * pBehaviorData() const;

private:
	string lcovariateName;
	int lperiod;
	ConstantCovariate * lpConstantCovariate;
	ChangingCovariate * lpChangingCovariate;
	BehaviorLongitudinalData * lpBehaviorData;
	const Network * lpFirstNetwork;
	const Network * lpSecondNetwork;
	TwoNetworkCache * lpTwoNetworkCache;

	// The current value of a behavior variable per each actor.
	// This array is 0 for covariate-based effects.
	const int * lvalues;
};


}

#endif /* COVARIATEMIXEDNETWORKALTERFUNCTION_H_ */
