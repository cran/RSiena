/******************************************************************************
 * SIENA: Simulation Investigation for Empirical Network Analysis
 *
 * Web: http://www.stats.ox.ac.uk/~snijders/siena/
 *
 * File: Chain.h
 *
 * Description: This file contains the definition of the Chain class.
 *****************************************************************************/

#ifndef CHAIN_H_
#define CHAIN_H_

#include <vector>
#include <map>
#include "model/ml/Option.h"

using namespace std;

namespace siena
{

// ----------------------------------------------------------------------------
// Section: Forward declarations
// ----------------------------------------------------------------------------

class MiniStep;
class Data;
class State;


// ----------------------------------------------------------------------------
// Section: Class definition
// ----------------------------------------------------------------------------


/**
 * Defines a sequence of ministeps, which is the basic structure for
 * Maximum-Likelihood calculations.
 */
class Chain
{
public:
	Chain(Data * pData);
	virtual ~Chain();

	// Initialization

	void setupInitialState(bool copyValues);

	// Chain modifications

	void clear();
	void insertBefore(MiniStep * pNewMiniStep, MiniStep * pExistingMiniStep);
	void remove(MiniStep * pMiniStep);
	void connect(int period);
	void onReciprocalRateChange(const MiniStep * pMiniStep, double newValue);
	void changeInitialState(const MiniStep * pMiniStep);

	// Accessors

	void period(int period);
	int period() const;
	MiniStep * pFirst() const;
	MiniStep * pLast() const;
	const State * pInitialState() const;
	int ministepCount() const;
	int diagonalMinistepCount() const;
	int consecutiveCancelingPairCount() const;
	int missingNetworkMiniStepCount() const;
	int missingBehaviorMiniStepCount() const;
	double mu() const;
	double sigma2() const;
	void printConsecutiveCancelingPairs() const;

	// Intervals

	int intervalLength(const MiniStep * pFirstMiniStep,
		const MiniStep * pLastMiniStep) const;

	// Same option related

	MiniStep * firstMiniStepForOption(const Option & rOption) const;
	MiniStep * nextMiniStepForOption(const Option & rOption,
		const MiniStep * pFirstMiniStep) const;

	// Random draws

	MiniStep * randomMiniStep() const;
	MiniStep * randomDiagonalMiniStep() const;
	MiniStep * randomMiniStep(MiniStep * pFirstMiniStep,
		MiniStep * pLastMiniStep) const;
	MiniStep * randomConsecutiveCancelingPair() const;
	MiniStep * randomMissingNetworkMiniStep() const;
	MiniStep * randomMissingBehaviorMiniStep() const;

	// copy
	Chain * copyChain() const;
	void dumpChain() const;

private:
	void resetOrderingKeys();
	void updateSameOptionPointersOnInsert(MiniStep * pMiniStep);
	void updateCCPs(MiniStep * pMiniStep);

	// A dummy first ministep in the chain
	MiniStep * lpFirst;

	// A dummy last ministep in the chain
	MiniStep * lpLast;

	// The underlying observed data
	Data * lpData;

	// The period of changes represented by this chain
	int lperiod;

	// The initial state of the variables (denoted y_init in the specification)
	State * lpInitialState;

	// Stores the ministeps in no particular order.
	// The first (dummy) ministep is not stored in this vector.

	vector<MiniStep *> lminiSteps;

	// Stores the diagonal ministeps in no particular order.
	vector<MiniStep *> ldiagonalMiniSteps;

	// Stores the first ministep of each CCP in no particular order.
	vector<MiniStep *> lccpMiniSteps;

	// Stores the networks ministeps of missing options
	vector<MiniStep *> lmissingNetworkMiniSteps;

	// Storesthe behavior ministeps of missing options
	vector<MiniStep *> lmissingBehaviorMiniSteps;

	// Sum of reciprocal rates over all non-dummy ministeps.
	double lmu;

	// Sum of squared reciprocal rates over all non-dummy ministeps.
	double lsigma2;

	// Maps each option to its first ministep in the chain (if any)
	map<const Option, MiniStep *> lfirstMiniStepPerOption;
};

}

#endif /* CHAIN_H_ */
