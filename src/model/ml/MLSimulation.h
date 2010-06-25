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


// ----------------------------------------------------------------------------
// Section: Class definition
// ----------------------------------------------------------------------------

class MLSimulation: public EpochSimulation
{
public:
	MLSimulation(Data * pData, Model * pModel);
	virtual ~MLSimulation();

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
	const Chain * pChain() const;
	void pChain(Chain * pChain);
	void pChainProbabilities(Chain * pChain, int period);

    // Metropolis-Hastings steps

	bool insertDiagonalMiniStep();
	bool cancelDiagonalMiniStep();
	bool permute(int c0);
	bool insertPermute(int c0);
	bool deletePermute(int c0);
	double proposalProbability() const;
	bool missingData() const;
	Aspect aspect() const;

	void simpleRates(bool flag);
	bool simpleRates() const;

	void insertDiagonalProbability(double probability);
	double insertDiagonalProbability() const;

	void cancelDiagonalProbability(double probability);
	double cancelDiagonalProbability() const;

	void permuteProbability(double probability);
	double permuteProbability() const;

	void insertPermuteProbability(double probability);
	double insertPermuteProbability() const;

	void deletePermuteProbability(double probability);
	double deletePermuteProbability() const;

	void randomMissingProbability(double probability);
	double randomMissingProbability() const;

	void missingNetworkProbability(double probability);
	double missingNetworkProbability() const;

	void missingBehaviorProbability(double probability);
	double missingBehaviorProbability() const;

	// Bayesian routines
	void initializeMCMCcycle();
	void MHPstep();
	int BayesAcceptances(unsigned iteration) const;
	double sampledBasicRates(unsigned iteration) const;
	void sampledBasicRates(double value);
	int sampledBasicRatesDistributions(unsigned iteration) const;
	void sampledBasicRatesDistributions(int value);
	double candidates(const EffectInfo * pEffect, unsigned iteration) const;
	void candidates(const EffectInfo * pEffect, double value);

private:
	void setStateBefore(MiniStep * pMiniStep);
	void resetVariables();

	Chain * lpChain;
	bool lsimpleRates;
	double linsertDiagonalProbability;
	double lcancelDiagonalProbability;
	double lpermuteProbability;
	double linsertPermuteProbability;
	double ldeletePermuteProbability;
	double lrandomMissingProbability;
	double lmissingNetworkProbability;
	double lmissingBehaviorProbability;
	double lproposalProbability;
	bool lmissingData;
	Aspect laspect;
	double lprobabilityArray[6];
	int lacceptances[6];
	int lrejections[6];
	vector<int> lBayesAcceptances;
	vector<double> lsampledBasicRates;
	vector<int> lsampledBasicRatesDistributions;
	map<const EffectInfo *, vector<double> > lcandidates;

};

}

#endif /* MLSIMULATION_H_ */
