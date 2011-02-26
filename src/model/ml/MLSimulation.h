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
	void connect(int period);
    void updateProbabilities(const Chain * pChain,
    	MiniStep * pFirstMiniStep,
    	MiniStep * pLastMiniStep);
    void executeMiniSteps(MiniStep * pFirstMiniStep, MiniStep * pLastMiniStep);
	void preburnin();
	void runEpoch(int period);
	void MLStep();
	void setUpProbabilityArray();

	int acceptances(int stepType) const;
	int rejections(int stepType) const;
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

	// Bayesian routines
	void initializeMCMCcycle();
	void MHPstep();
	int BayesAcceptances(unsigned iteration) const;
	double candidates(const EffectInfo * pEffect, unsigned iteration) const;
	void candidates(const EffectInfo * pEffect, double value);

private:
	void setStateBefore(MiniStep * pMiniStep);
	void resetVariables();
	bool validInsertMissingStep(const Option * pOption,
		int d0,
		const MiniStep * pMiniStepA);
	bool validDeleteMissingStep(MiniStep * pMiniStepA, bool applyTwice);
	MiniStep * createMiniStep(const Option * pOption,
		int difference = 0) const;

	Chain * lpChain;
	bool lsimpleRates;
	double lproposalProbability;
	bool lmissingData;
	Aspect laspect;
	double lprobabilityArray[7];
	int lacceptances[7];
	int lrejections[7];
	vector<int> lBayesAcceptances;
	map<const EffectInfo *, vector<double> > lcandidates;
	double lmissingNetworkProbability;
	double lmissingBehaviorProbability;
	// current length of permuted interval
	double lcurrentPermutationLength;
	int lthisPermutationLength;
	bool lphase1;

	// A vector of options with missing values in the initial observation
	vector<const Option *> linitialMissingOptions;
};

}

#endif /* MLSIMULATION_H_ */
