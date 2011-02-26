/******************************************************************************
 * SIENA: Simulation Investigation for Empirical Network Analysis
 *
 * Web: http://www.stats.ox.ac.uk/~snijders/siena/
 *
 * File: Model.cpp
 *
 * Description: This file contains the implementation of the Model class.
 *****************************************************************************/
#include <R_ext/Print.h>
#include <vector>
#include "Model.h"
#include "utils/Utils.h"
#include "data/Data.h"
#include "data/LongitudinalData.h"
#include "model/EffectInfo.h"
#include "model/variables/DependentVariable.h"
#include "model/effects/AllEffects.h"
#include "model/ml/Chain.h"

namespace siena
{

// ----------------------------------------------------------------------------
// Section: Constructors and initializers
// ----------------------------------------------------------------------------

/**
 * Creates a new actor-based model.
 */
Model::Model()
{
	this->lconditional = false;
	this->lneedChain = false;
	this->lneedScores = false;
	this->lneedDerivatives = false;
	this->lneedChangeContributions = false;
	this->lparallelRun = false;
	this->linsertDiagonalProbability = 0;
	this->lcancelDiagonalProbability = 0;
	this->lpermuteProbability = 0;
	this->linsertPermuteProbability = 0;
	this->ldeletePermuteProbability = 0;
	this->linsertRandomMissingProbability = 0;
	this->ldeleteRandomMissingProbability = 0;
	this->lsimpleRates = 0;
	this->lmodelType = NORMAL;
}


/**
 * Deallocates this model.
 */
Model::~Model()
{
	// Delete the arrays of basic rate parameters

	while (!this->lbasicRateParameters.empty())
	{
		double * array = this->lbasicRateParameters.begin()->second;
		this->lbasicRateParameters.erase(this->lbasicRateParameters.begin());
		delete[] array;
	}

	deallocateVector(this->leffects);

	// Delete the arrays of target changes

	while (!this->ltargetChanges.empty())
	{
		int * array = this->ltargetChanges.begin()->second;
		this->ltargetChanges.erase(this->ltargetChanges.begin());
		delete[] array;
	}

	for (unsigned i = 0; i < this->lchainStore.size(); i++)
	{
		deallocateVector(this->lchainStore[i]);
	}
}


// ----------------------------------------------------------------------------
// Section: Accessors
// ----------------------------------------------------------------------------

/**
 * Stores if conditional simulation has to be carried out.
 */
void Model::conditional(bool flag)
{
	this->lconditional = flag;
}


/**
 * Returns if conditional simulation has to be carried out.
 */
bool Model::conditional() const
{
	return this->lconditional;
}


/**
 * Stores the conditional variable name for conditional simulation
 */
void Model::conditionalDependentVariable(string variableName)
{
	this->lconditionalDependentVariable = variableName;
}


/**
 * Returns the conditional variable name for conditional simulation
 */
string Model::conditionalDependentVariable() const
{
	return this->lconditionalDependentVariable;
}
/**
 * Stores if a chain is to be built in the current simulation
 */
void Model::needChain(bool flag)
{
	this->lneedChain = flag;
}


/**
 * Returns if a chain is to be built in the current simulation
 */
bool Model::needChain() const
{
	return this->lneedChain;
}

/**
 * Stores if scores are to accumulated in the current simulation
 */
void Model::needScores(bool flag)
{
	this->lneedScores = flag;
}


/**
 * Returns if scores are to accumulated in the current simulation
 */
bool Model::needScores() const
{
	return this->lneedScores;
}
/**
 * Stores if derivatives are to accumulated in the current simulation
 */
void Model::needDerivatives(bool flag)
{
	this->lneedDerivatives = flag;
}


/**
 * Returns if derivatives are to accumulated in the current simulation
 */
bool Model::needDerivatives() const
{
	return this->lneedDerivatives;
}
/**
 * Stores if change contributions are to be stored on ministeps
 */
void Model::needChangeContributions(bool flag)
{
	this->lneedChangeContributions = flag;
}


/**
 * Returns if change contributions are to be stroed on ministeps
 */
bool Model::needChangeContributions() const
{
	return this->lneedChangeContributions;
}

/**
 * Stores if this is a parallel run
 */
void Model::parallelRun(bool flag)
{
	this->lparallelRun = flag;
}


/**
 * Returns if this is a parallel run
 */
bool Model::parallelRun() const
{
	return this->lparallelRun;
}

/**
 * Stores the number of ML steps
 */
void Model::numberMLSteps(int value)
{
	this->lnumberMLSteps = value;
}


/**
 * Returns the number of ML steps
 */
int Model::numberMLSteps() const
{
	return this->lnumberMLSteps;
}

/**
 * Stores the maximum permutation length
 */
void Model::maximumPermutationLength(double value)
{
	this->lmaximumPermutationLength = value;
}


/**
 * Returns the maximum permutation length
 */
double Model::maximumPermutationLength() const
{
	return this->lmaximumPermutationLength;
}

/**
 * Stores the minimum permutation length
 */
void Model::minimumPermutationLength(double value)
{
	this->lminimumPermutationLength = value;
}


/**
 * Returns the minimum permutation length
 */
double Model::minimumPermutationLength() const
{
	return this->lminimumPermutationLength;
}

/**
 * Stores the initial permutation length
 */
void Model::initialPermutationLength(double value)
{
	this->linitialPermutationLength = value;
}


/**
 * Returns the initial permutation length
 */
double Model::initialPermutationLength() const
{
	return this->linitialPermutationLength;
}

/**
 * Stores the number of MH batches
 */
void Model::numberMHBatches(int value)
{
	this->lnumberMHBatches = value;
}

/**
 * Returns the number of MH batches
 */
int Model::numberMHBatches() const
{
	return this->lnumberMHBatches;
}
/**
 * Stores the Bayesian scale factor
 */
void Model::BayesianScaleFactor(double value)
{
	this->lBayesianScaleFactor = value;
}


/**
 * Returns the Bayesian scale factor
 */
double Model::BayesianScaleFactor() const
{
	return this->lBayesianScaleFactor;
}
// ----------------------------------------------------------------------------
// Section: Effect management
// ----------------------------------------------------------------------------

/**
 * Stores the basic rate parameter for the given dependent variable at the
 * given period.
 */
void Model::basicRateParameter(LongitudinalData * pDependentVariableData,
	int period,
	double value)
{
	if (!this->lbasicRateParameters[pDependentVariableData])
	{
		double * array =
			//	new double[pDependentVariableData->observationCount() - 1];
			new double[pDependentVariableData->observationCount() ];

		// The default basic rate is 1.

		for (int i = 0;
			i < pDependentVariableData->observationCount() - 1;
			i++)
		{
			array[i] = 1;
		}

		this->lbasicRateParameters[pDependentVariableData] = array;
	}
	this->lbasicRateParameters[pDependentVariableData][period] = value;
}


/**
 * Returns the basic rate parameter for the given dependent variable at the
 * given period.
 */
double Model::basicRateParameter(LongitudinalData * pDependentVariableData,
	int period) const
{
	std::map<const LongitudinalData *, double *>::const_iterator iter =
		this->lbasicRateParameters.find(pDependentVariableData);
	double value = 1;

	if (iter != this->lbasicRateParameters.end())
	{
		value = iter->second[period];
	}

	return value;
}


/**
 * Adds a new effect to this model and returns the parameters wrapped into
 * an EffectInfo object for future reference.
 * @param[in] variableName the name of the variable this effect is associated
 * with
 * @param[in] effectName the name of the effect
 * @param[in] effectType the type of the effect ("rate", "eval", or "endow")
 * @param[in] parameter the multiplicative weight of the effect
 * @param[in] internalEffectParameter the internal effect parameter
 * (if applicable)
 * @param[in] interactionName1 the name of another variable or covariate
 * (if any) the effect interacts with. If the effect interacts with two
 * variables or covariates, this parameter specifies one of them.
 * @param[in] interactionName2 the name of the other interacting variable or
 * covariate, if the effect has two such interactions. We should make sure
 * that the order of these interaction variables is the same both in R and C++.
 * @param[in] rateType distinguishes between structural rate effects and
 * covariate rate effects
 */
EffectInfo * Model::addEffect(string variableName,
	string effectName,
	string effectType,
	double parameter,
	double internalEffectParameter,
	string interactionName1,
	string interactionName2,
	string rateType)
{
	EffectInfo * pInfo = new EffectInfo(variableName,
		effectName,
		effectType,
		parameter,
		internalEffectParameter,
		interactionName1,
		interactionName2,
		rateType);

	this->leffects.push_back(pInfo);

	if (effectType == "rate")
	{
		this->lrateEffects[variableName].push_back(pInfo);
	}
	else if (effectType == "eval")
	{
		this->levaluationEffects[variableName].push_back(pInfo);
	}
	else if (effectType == "endow")
	{
		this->lendowmentEffects[variableName].push_back(pInfo);
	}
	else
	{
		throw invalid_argument("Unexpected effect type '" + effectType + "'.");
	}

	return pInfo;
}


/**
 * Adds a new interaction effect between two or three other effects to this model.
 * The parameters are wrapped into an EffectInfo object and returned for future
 * reference.
 * @param[in] variableName the name of the variable this effect is associated
 * with
 * @param[in] effectName the name of the effect
 * @param[in] effectType the type of the effect ("rate", "eval", or "endow")
 * @param[in] parameter the multiplicative weight of the effect
 * @param[in] pEffect1 the first of the interacting effects
 * @param[in] pEffect2 the second of the interacting effects
 * @param[in] pEffect3 the third of the interacting effects (may be undefined)
 */
EffectInfo * Model::addInteractionEffect(string variableName,
	string effectName,
	string effectType,
	double parameter,
	const EffectInfo * pEffect1,
	const EffectInfo * pEffect2,
	const EffectInfo * pEffect3)
{
	EffectInfo * pInfo = new EffectInfo(variableName,
		effectName,
		effectType,
		parameter,
		pEffect1,
		pEffect2,
		pEffect3);

	this->leffects.push_back(pInfo);

	if (effectType == "rate")
	{
		this->lrateEffects[variableName].push_back(pInfo);
	}
	else if (effectType == "eval")
	{
		this->levaluationEffects[variableName].push_back(pInfo);
	}
	else if (effectType == "endow")
	{
		this->lendowmentEffects[variableName].push_back(pInfo);
	}
	else
	{
		throw invalid_argument("Unexpected effect type '" + effectType + "'.");
	}

	return pInfo;
}


/**
 * Returns the rate effects for the given dependent variable (excluding the
 * basic rate effects, which are treated differently).
 */
const vector<EffectInfo *> & Model::rRateEffects(string variableName) const
{
	map<string, vector<EffectInfo *> >::const_iterator iter =
		this->lrateEffects.find(variableName);

	if (iter == this->lrateEffects.end())
	{
		return this->lemptyEffectVector;
	}

	return iter->second;
}


/**
 * Returns the evaluation effects for the given dependent variable.
 */
const vector<EffectInfo *> & Model::rEvaluationEffects(string variableName)
	const
{
	map<string, vector<EffectInfo *> >::const_iterator iter =
		this->levaluationEffects.find(variableName);

	if (iter == this->levaluationEffects.end())
	{
		return this->lemptyEffectVector;
	}

	return iter->second;
}


/**
 * Returns the endowment effectvs for the given dependent variable.
 */
const vector<EffectInfo *> & Model::rEndowmentEffects(string variableName)
	const
{
	map<string, vector<EffectInfo *> >::const_iterator iter =
		this->lendowmentEffects.find(variableName);

	if (iter == this->lendowmentEffects.end())
	{
		return this->lemptyEffectVector;
	}

	return iter->second;
}


// ----------------------------------------------------------------------------
// Section: Target change
// ----------------------------------------------------------------------------

/**
 * Stores the target change for the given period to be used in conditional
 * simulation. Since the change is specific to a certain data object, the data
 * object should be provided as a parameter too.
 */
void Model::targetChange(const Data * pData, int period, int change)
{
	if (!this->ltargetChanges[pData])
	{
		int * array =
			new int[pData->observationCount() - 1];

		// Initialize the array

		for (int i = 0; i < pData->observationCount() - 1; i++)
		{
			array[i] = 0;
		}

		this->ltargetChanges[pData] = array;
	}

	this->ltargetChanges[pData][period] = change;
}


/**
 * Returns the target change for the given period to be used in conditional
 * simulation. Since the change is specific to a certain data object, the data
 * object should be provided as a parameter too.
 */
int Model::targetChange(const Data * pData, int period) const
{
	map<const Data *, int *>::const_iterator iter =
		this->ltargetChanges.find(pData);
	int value = 1;

	if (iter != this->ltargetChanges.end())
	{
		value = iter->second[period];
	}

	return value;
}

// ----------------------------------------------------------------------------
// Section: Chain storage
// ----------------------------------------------------------------------------

/**
 * Stores a copy of the chain in the vector of stored chains for this period.
 */
void Model::chainStore(const Chain& chain, int periodFromStart)
{
//	Rprintf(" %d %d\n", this->lchainStore.size(), periodFromStart);
	if (this->lchainStore.size() == 0)
	{
//	Rprintf(" %d\n", this->lnumberOfPeriods);
		this->setupChainStore(this->lnumberOfPeriods);
	}
	// make a copy of the chain
	Chain * pChain = chain.copyChain();

    this->lchainStore[periodFromStart].push_back(pChain);
}


/**
 * Returns reference to the vector of stored chains for this period.
 */
vector<Chain *> & Model::rChainStore(int periodFromStart)
{
	return this->lchainStore[periodFromStart];
}

void Model::partClearChainStore()
{
	for (unsigned i = 0; i < this->lchainStore.size(); i++)
	{
		// need to keep the final one to continue from
		for (unsigned j = 0; j < this->lchainStore[i].size() - 1; j++)
		{
			//delete the chain pointed to
			delete this->lchainStore[i][j];
		}
		vector<Chain *>::iterator iter = this->lchainStore[i].begin();
		for (unsigned j = 0; j < this->lchainStore[i].size() - 1; j++)
		{
			//erase the entry
			this->lchainStore[i].erase(iter);
			iter++;
		}

	}
}
void Model::clearChainStore()
{
	for (unsigned i = 0; i < this->lchainStore.size(); i++)
	{
		deallocateVector(this->lchainStore[i]);
	}

}

void Model::clearChainStore(int periodFromStart)
{
	deallocateVector(this->lchainStore[periodFromStart]);
}

void Model::setupChainStore(int numberPeriods)
{
	this->lchainStore.resize(numberPeriods);
}

void Model::deleteLastChainStore(int periodFromStart)
{
	this->lchainStore[periodFromStart].pop_back();
}

void Model::numberOfPeriods(int numberOfPeriods)
{
	this->lnumberOfPeriods = numberOfPeriods;
}
int Model::numberOfPeriods()
{
	return this->lnumberOfPeriods;
}

void Model::modelType(int type)
{
	this->lmodelType = ModelType(type);
}

ModelType Model::modelType() const
{
	return this->lmodelType;
}
bool Model::modelTypeB() const
{
	return this->lmodelType == BFORCE || 
		this->lmodelType == BAGREE || this->lmodelType == BJOINT;
}

// ----------------------------------------------------------------------------
// Section: Probabilities
//-----------------------------------------------------------------------------
/**
 * Stores the probability associated with the insertDiagonalMiniStep
 * operation.
 */
void Model::insertDiagonalProbability(double probability)
{
	this->linsertDiagonalProbability = probability;
}


/**
 * Returns the probability associated with the insertDiagonalMiniStep
 * operation.
 */
double Model::insertDiagonalProbability() const
{
	return this->linsertDiagonalProbability;
}


/**
 * Stores the probability associated with the cancelDiagonalMiniStep
 * operation.
 */
void Model::cancelDiagonalProbability(double probability)
{
	this->lcancelDiagonalProbability = probability;
}


/**
 * Returns the probability associated with the cancelDiagonalMiniStep
 * operation.
 */
double Model::cancelDiagonalProbability() const
{
	return this->lcancelDiagonalProbability;
}


/**
 * Stores the probability associated with the permute
 * operation.
 */
void Model::permuteProbability(double probability)
{
	this->lpermuteProbability = probability;
}


/**
 * Returns the probability associated with the permute
 * operation.
 */
double Model::permuteProbability() const
{
	return this->lpermuteProbability;
}

/**
 * Stores the probability associated with the insertPermute
 * operation.
 */
void Model::insertPermuteProbability(double probability)
{
	this->linsertPermuteProbability = probability;
}


/**
 * Returns the probability associated with the insertPermute
 * operation.
 */
double Model::insertPermuteProbability() const
{
	return this->linsertPermuteProbability;
}


/**
 * Stores the probability associated with the deletePermute
 * operation.
 */
void Model::deletePermuteProbability(double probability)
{
	this->ldeletePermuteProbability = probability;
}


/**
 * Returns the probability associated with the deletePermute
 * operation.
 */
double Model::deletePermuteProbability() const
{
	return this->ldeletePermuteProbability;
}

/**
 * Stores the probability associated with the insertRandomMissing
 * operation.
 */
void Model::insertRandomMissingProbability(double probability)
{
	this->linsertRandomMissingProbability = probability;
}


/**
 * Returns the probability associated with the insertRandomMissing
 * operation.
 */
double Model::insertRandomMissingProbability() const
{
	return this->linsertRandomMissingProbability;
}

/**
 * Stores the probability associated with the deleteRandomMissing
 * operation.
 */
void Model::deleteRandomMissingProbability(double probability)
{
	this->ldeleteRandomMissingProbability = probability;
}


/**
 * Returns the probability associated with the deleteRandomMissing
 * operation.
 */
double Model::deleteRandomMissingProbability() const
{
	return this->ldeleteRandomMissingProbability;
}

/**
 * Adds the missing network probability (prmin) to the vector for
 * the specified group
 */
void Model::missingNetworkProbability(double probability)
{
	this->lmissingNetworkProbability.push_back(probability);
}


/**
 * Returns the missing network probability (prmin) for the specified group
 * and period
 */
double Model::missingNetworkProbability( int periodFromStart) const
{
	return this->lmissingNetworkProbability[periodFromStart];
}

/**
 * Adds the missing behavior probability (prmib) to the vector for
 * the specified group.
 */
void Model::missingBehaviorProbability(double probability)
{
	this->lmissingBehaviorProbability.push_back(probability);
}


/**
 * Returns the missing behavior probability (prmib) for the specified
 * period from the start
 */
double Model::missingBehaviorProbability(int periodFromStart) const
{
	return this->lmissingBehaviorProbability[periodFromStart];
}
/**
 * Stores the simpleRates flag for the ML model.
 */
void Model::simpleRates(bool simpleRates)
{
	this->lsimpleRates = simpleRates;
}


/**
 * Returns the simpleRates flag for the ML model.
 */
bool Model::simpleRates() const
{
	return this->lsimpleRates;
}
}
