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

	void chainStore(Chain& chain);
	vector<Chain *> *chainStore();
	
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

	void needChangeContributions(bool flag);
	bool needChangeContributions() const;

	void parallelRun(bool flag);
	bool parallelRun() const;

	void numberMLSteps(int value);
	int numberMLSteps() const;

	void maximumPermutationLength(double value);
	double maximumPermutationLength() const;

	void minimumPermutationLength(double value);
	double minimumPermutationLength() const;

	void initialPermutationLength(double value);
	double initialPermutationLength() const;

	void numberMHBatches(int value);
	int numberMHBatches() const;

	void BayesianScaleFactor(double value);
	double BayesianScaleFactor() const;

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

	// A vector of evaluation effects per variable
	map<string, vector<EffectInfo *> > levaluationEffects;

	// A vector of endowment effects per variable
	map<string, vector<EffectInfo *> > lendowmentEffects;

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

	// number of steps in a MCMC run for ML
	int lnumberMHBatches;

	// Bayesian scale factor for normal random variates
	double lBayesianScaleFactor;

	// chain storage
	vector <Chain *> lchainStore;

};

}

#endif /*MODEL_H_*/
