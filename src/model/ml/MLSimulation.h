#ifndef MLSIMULATION_H_
#define MLSIMULATION_H_

#include "model/EpochSimulation.h"

namespace siena
{

// ----------------------------------------------------------------------------
// Section: Enums
// ----------------------------------------------------------------------------

enum Aspect {NETWORK, BEHAVIOR};

/**
 * This enumeration defines the possible types of MH step
 */

enum MHStepType {INSDIAG, CANCDIAG, PERMUTE, INSPERM, DELPERM, INSMISS,
				 DELMISS, INSMISDAT, DELMISDAT, NBRTYPES};

// ----------------------------------------------------------------------------
// Section: Forward declarations
// ----------------------------------------------------------------------------

class MiniStep;
class Option;
class DependentVariable;
class NetworkVariable;
    
/**
 * This class provides the functionality necessary for simulating an ML model
 * between two observations.
 */

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

	void missingNetworkProbability(double probability);
	double missingNetworkProbability() const;

	void missingBehaviorProbability(double probability);
	double missingBehaviorProbability() const;

	void currentPermutationLength(double value);
	double currentPermutationLength() const;

	void updateCurrentPermutationLength(bool accept);

	void createEndStateDifferences();

	void recordOutcome(const MiniStep & miniStep, bool accept,
		int stepType, bool misdat);

    bool neighbourhoodChange(MiniStep * pMiniStep1, MiniStep * pMiniStep2,
                             DependentVariable * pVariable,
                             NetworkVariable * pNetworkVariable,
                             int ego1, int alter1);
							 
	void gotoLastState(); 
    

private:
	void setStateBefore(MiniStep * pMiniStep);
	void resetVariables();
	bool validInsertMissingStep(const Option * pOption,
		int d0,
		const MiniStep * pMiniStepA);
	bool validDeleteMissingStep(MiniStep * pMiniStepA, bool applyTwice);
	MiniStep * createMiniStep(const Option * pOption,
		int difference, bool value) const;

	double lproposalProbability;
	bool lmissingData;
	Aspect laspect;
	double lprobabilityArray[7];
	int lacceptances[NBRTYPES];
	int lrejections[NBRTYPES];
	int laborted[NBRTYPES];
	double lmissingNetworkProbability;
	double lmissingBehaviorProbability;
	// current length of permuted interval
	double lcurrentPermutationLength;
	int lthisPermutationLength;

	// A vector of options with missing values in the initial observation
	std::vector<const Option *> linitialMissingOptions;
};

}

#endif /* MLSIMULATION_H_ */
