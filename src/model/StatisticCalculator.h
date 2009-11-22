/******************************************************************************
 * SIENA: Simulation Investigation for Empirical Network Analysis
 *
 * Web: http://www.stats.ox.ac.uk/~snijders/siena/
 *
 * File: StatisticCalculator.h
 *
 * Description: This file contains the definition of the
 * StatisticCalculator class.
 *****************************************************************************/

#ifndef STATISTICCALCULATOR_H_
#define STATISTICCALCULATOR_H_

#include <map>

using namespace std;

namespace siena
{

// ----------------------------------------------------------------------------
// Section: Forward declarations
// ----------------------------------------------------------------------------

class Data;
class Model;
class State;
class EffectInfo;
class LongitudinalData;
class NetworkLongitudinalData;
class BehaviorLongitudinalData;
class Network;
class EffectFactory;


// ----------------------------------------------------------------------------
// Section: StatisticCalculator class
// ----------------------------------------------------------------------------

/**
 * Provides means for calculating statistics corresponding to effects of a
 * certain model. In order to calculate the statistics for effects of a
 * model, one has to provide the underlying observed data, the period under
 * consideration, and the current state of all dependent variables. This
 * information should be given at the construction of the statistic calculator,
 * after which individual statistics can be queried. Example:
 *
 * StatisticCalculator calculator(pData, pModel, pState, period);
 * double statistic1 = calculator.statistic(pEffectInfo1);
 * double statistic2 = calculator.statistic(pEffectInfo2);
 */
class StatisticCalculator
{
public:
	StatisticCalculator(const Data * pData,
		const Model * pModel,
		State * pState,
		int period);
	virtual ~StatisticCalculator();

	double statistic(EffectInfo * pEffectInfo) const;
	int distance(LongitudinalData * pData, int period) const;

private:
	void calculateStatistics();
	void calculateNetworkRateStatistics(
		NetworkLongitudinalData * pNetworkData);
	void calculateNetworkEvaluationStatistics(
		NetworkLongitudinalData * pNetworkData);
	void calculateNetworkEndowmentStatistics(
		NetworkLongitudinalData * pNetworkData);
	void calculateBehaviorStatistics(BehaviorLongitudinalData * pBehaviorData);
	void calculateBehaviorRateStatistics(BehaviorLongitudinalData * pBehaviorData);
	void subtractNetwork(Network * pNetwork,
		const Network * pMissingTieNetwork) const;

	void replaceNetwork(Network * pNetwork,
		const Network * pValueNetwork,
		const Network * pStructuralTieNetwork) const;

	// The data to be used for calculating the statistics
	const Data * lpData;

	// The model containing the effects whose statistics are to be calculated
	const Model * lpModel;

	// The current state of the dependent variables
	State * lpState;

	// The period of the evolution we are interested in
	int lperiod;

	// The resulting map of statistic values
	map<EffectInfo *, double> lstatistics;

	// Array of simulated distances per variable
	map<LongitudinalData *, int *> ldistances;

	State * lpPredictorState;
};

}

#endif /*STATISTICCALCULATOR_H_*/
