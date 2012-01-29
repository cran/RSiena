/******************************************************************************
 * SIENA: Simulation Investigation for Empirical Network Analysis
 *
 * Web: http://www.stats.ox.ac.uk/~snijders/siena/
 *
 * File: Model.h
 *
 * Description: This file contains the definition of the Model class.
 *****************************************************************************/

#ifndef MODEL_H_
#define MODEL_H_

#include <map>
#include <vector>
#include <string>

using namespace std;

namespace siena
{

// ----------------------------------------------------------------------------
// Section: Enumerations
// ----------------------------------------------------------------------------

/**
 * This enumeration defines the possible types of model for symmetric,
 * undirected networks.
 */
	enum ModelType { NOTUSED, NORMAL, AFORCE, AAGREE, BFORCE, BAGREE, BJOINT };
// ----------------------------------------------------------------------------
// Section: Forward declarations
// ----------------------------------------------------------------------------

class Data;
class LongitudinalData;
class Effect;
class StructuralRateEffect;
class Function;
class EffectInfo;
class Chain;


// ----------------------------------------------------------------------------
// Section: Model class
// ----------------------------------------------------------------------------

/**
 * This class defines the actor-based models for longitudinal network and
 * behavioral data.
 */
class Model
{
public:
	Model();
	virtual ~Model();

	// Basic rate effects

	void basicRateParameter(LongitudinalData * pDependentVariableData,
		int period,
		double value);
	double basicRateParameter(LongitudinalData * pDependentVariableData,
		int period) const;

	// Other effects

	EffectInfo * addEffect(string variableName,
		string effectName,
		string effectType,
		double parameter,
		double internalEffectParameter = 0,
		string interactionName1 = "",
		string interactionName2 = "",
		string rateType = "");
	EffectInfo * addInteractionEffect(string variableName,
		string effectName,
		string effectType,
		double parameter,
		const EffectInfo * pEffect1,
		const EffectInfo * pEffect2,
		const EffectInfo * pEffect3 = 0);

	const vector<EffectInfo *> & rRateEffects(string variableName) const;
	const vector<EffectInfo *> & rEvaluationEffects(string variableName) const;
	const vector<EffectInfo *> & rEndowmentEffects(string variableName) const;
	const vector<EffectInfo *> & rCreationEffects(string variableName) const;

	void chainStore(const Chain& chain, int periodFromStart);
	vector <Chain *> & rChainStore( int periodFromStart);
	void clearChainStore(int keep, int periodFromStart);
	void setupChainStore(int numberOfPeriods);
	void deleteLastChainStore(int periodFromStart);
	void numberOfPeriods(int numberOfPeriods);
	int numberOfPeriods();

	// Various flags

	void conditional(bool flag);
	bool conditional() const;

	int targetChange(const Data * pData, int period) const;
	void targetChange(const Data * pData, int period, int change);

	string conditionalDependentVariable() const;
	void conditionalDependentVariable(string variableName);

	void needChain(bool flag);
	bool needChain() const;

	void needScores(bool flag);
	bool needScores() const;

	void needDerivatives(bool flag);
	bool needDerivatives() const;

	void parallelRun(bool flag);
	bool parallelRun() const;

	void modelType(int type);
	ModelType modelType() const;
	bool modelTypeB() const;

	// various stores for ML

	void numberMLSteps(int value);
	int numberMLSteps() const;

	void maximumPermutationLength(double value);
	double maximumPermutationLength() const;

	void minimumPermutationLength(double value);
	double minimumPermutationLength() const;

	void initialPermutationLength(double value);
	double initialPermutationLength() const;

	void initializeCurrentPermutationLength();
	double currentPermutationLength(int period) const;
	void currentPermutationLength(int period, double value);

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

	void insertRandomMissingProbability(double probability);
	double insertRandomMissingProbability() const;

	void deleteRandomMissingProbability(double probability);
	double deleteRandomMissingProbability() const;

	void missingNetworkProbability(double probability);
	double missingNetworkProbability(int periodFromStart) const;

	void missingBehaviorProbability(double probability);
	double missingBehaviorProbability(int periodFromStart) const;

	// simple rates flag for ML
	bool simpleRates() const;
	void simpleRates(bool simpleRates);

private:
	// Indicates if conditional simulation has to be carried out
	bool lconditional;

	// name of conditional dependent variable
	string lconditionalDependentVariable;

	// Targets for conditional dependent variable per each data object
	// and period.

	map<const Data *, int *> ltargetChanges;

	// An array of doubles per each longitudinal data object storing
	// the basic rate parameters for all periods

	map<const LongitudinalData *, double *> lbasicRateParameters;

	// A vector of effects other than the basic rate effects.
	vector<EffectInfo *> leffects;

	// A vector of rate effects (except the basic rate effects) per variable
	map<string, vector<EffectInfo *> > lrateEffects;

	// A vector of pointers to evaluation effects per variable
	map<string, vector<EffectInfo *> > levaluationEffects;

	// A vector of pointers to endowment effects per variable
	map<string, vector<EffectInfo *> > lendowmentEffects;

	// A vector of pointers to creation effects per variable
	map<string, vector<EffectInfo *> > lcreationEffects;

	// A dummy vector of effect infos in case we need a reference to
	// non-existent vectors

	const vector<EffectInfo *> lemptyEffectVector;

	// indicates whether we need to keep a chain of ministeps
	bool lneedChain;

	// indicates whether we need to accumulate scores in this iteration
	bool lneedScores;

	// indicates whether we need to accumulate derivatives for ML in
	// this iteration
	bool lneedDerivatives;

	// indicates whether we need to store the change contributions on
	// the ministeps
	bool lneedChangeContributions;

	// indicates whether we need to match Siena3 in use of random variables
	// and score calculations
	bool lparallelRun;

	// number of steps in a run for ML
	int lnumberMLSteps;

	// maximum length of permuted interval
	double lmaximumPermutationLength;

	// minimum length of permuted interval
	double lminimumPermutationLength;

	// initial length of permuted interval
	double linitialPermutationLength;

	// current length of permuted interval: varies by period
	vector <double> lcurrentPermutationLength;

	// probabilities of the different ML steps
	double linsertDiagonalProbability;
	double lcancelDiagonalProbability;
	double lpermuteProbability;
	double linsertPermuteProbability;
	double ldeletePermuteProbability;
	double linsertRandomMissingProbability;
	double ldeleteRandomMissingProbability;

	bool lsimpleRates;

	vector <double> lmissingNetworkProbability;
	vector <double> lmissingBehaviorProbability;

	// chain storage: vector of chains for each period for each set of samples
	// lchainStore[i] is set of entries for periodFromStart <i>,
	// which incorporates both the group and period.
	vector <vector <Chain *> > lchainStore;

	int lnumberOfPeriods;

	ModelType lmodelType;
};

}

#endif /*MODEL_H_*/
