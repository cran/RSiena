#ifndef MLSIMULATION_H_
#define MLSIMULATION_H_

#include "model/EpochSimulation.h"

namespace siena
{

// ----------------------------------------------------------------------------
// Section: Enums
// ----------------------------------------------------------------------------

enum Aspect {NETWORK, BEHAVIOR};


// ----------------------------------------------------------------------------
// Section: Forward declarations
// ----------------------------------------------------------------------------

class Chain;
class MiniStep;
class Option;


// ----------------------------------------------------------------------------
// Section: Class definition
// ----------------------------------------------------------------------------

class MLSimulation: public EpochSimulation
{
public:
	MLSimulation(Data * pData, Model * pModel);
	virtual ~MLSimulation();

    void initialize(int period);
	void initializeInitialState(int period);
	void connect(int period);
	void preburnin();
	void runEpoch(int period);
	void MLStep();
	void setUpProbabilityArray();
    void updateProbabilities(Chain * pChain,
    	MiniStep * pFirstMiniStep,
    	MiniStep * pLastMiniStep);
    void executeMiniSteps(MiniStep * pFirstMiniStep, MiniStep * pLastMiniStep);

	int acceptances(int stepType) const;
	int rejections(int stepType) const;
	int aborted(int stepType) const;
	Chain * pChain() const;
	void pChain(Chain * pChain);
	void pChainProbabilities(Chain * pChain, int period);

    // Metropolis-Hastings steps

	bool insertDiagonalMiniStep();
	bool cancelDiagonalMiniStep();
	bool permute(int c0);
	bool insertPermute(int c0);
	bool deletePermute(int c0);
	bool insertMissing();
	bool deleteMissing();
	double proposalProbability() const;
	bool missingData() const;
	Aspect aspect() const;

	void simpleRates(bool flag);
	bool simpleRates() const;

	void missingNetworkProbability(double probability);
	double missingNetworkProbability() const;

	void missingBehaviorProbability(double probability);
	double missingBehaviorProbability() const;

	void currentPermutationLength(int period, double value);
	double currentPermutationLength(int period) const;

	void updateCurrentPermutationLength(bool accept);

	void createEndStateDifferences();

private:
	void setStateBefore(MiniStep * pMiniStep);
	void resetVariables();
	bool validInsertMissingStep(const Option * pOption,
		int d0,
		const MiniStep * pMiniStepA);
	bool validDeleteMissingStep(MiniStep * pMiniStepA, bool applyTwice);
	MiniStep * createMiniStep(const Option * pOption,
		int difference, bool value) const;

	Chain * lpChain;
	bool lsimpleRates;
	double lproposalProbability;
	bool lmissingData;
	Aspect laspect;
	double lprobabilityArray[7];
	int lacceptances[7];
	int lrejections[7];
	int laborted[7];
	double lmissingNetworkProbability;
	double lmissingBehaviorProbability;
	// current length of permuted interval
	double lcurrentPermutationLength;
	int lthisPermutationLength;


	// A vector of options with missing values in the initial observation
	vector<const Option *> linitialMissingOptions;
};

}

#endif /* MLSIMULATION_H_ */
