/******************************************************************************
 * SIENA: Simulation Investigation for Empirical Network Analysis
 *
 * Web: http://www.stats.ox.ac.uk/~snijders/siena/
 *
 * File: Model.cpp
 *
 * Description: This file contains the implementation of the Model class.
 *****************************************************************************/

#include "Model.h"
#include "utils/Utils.h"
#include "data/Data.h"
#include "data/LongitudinalData.h"
#include "model/EffectInfo.h"
#include "model/variables/DependentVariable.h"
#include "model/effects/AllEffects.h"

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
	this->lneedScores = false;
	this->lparallelRun = false;
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

}
