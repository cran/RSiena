/******************************************************************************
 * SIENA: Simulation Investigation for Empirical Network Analysis
 *
 * Web: http://www.stats.ox.ac.uk/~snijders/siena/
 *
 * File: DependentVariable.cpp
 *
 * Description: This file contains the implementation of the
 * DependentVariable class.
 *****************************************************************************/

#include <cmath>
#include <stdexcept>
#include <R.h>

#include "BehaviorVariable.h"
#include "DependentVariable.h"
#include "utils/Utils.h"
#include "data/ActorSet.h"
#include "data/ConstantCovariate.h"
#include "data/ChangingCovariate.h"
#include "data/Network.h"
#include "data/OneModeNetwork.h"
#include "data/LongitudinalData.h"
#include "data/NetworkLongitudinalData.h"
#include "model/Model.h"
#include "model/EffectInfo.h"
#include "model/EpochSimulation.h"
#include "model/effects/AllEffects.h"
#include "model/effects/EffectFactory.h"
#include "model/variables/NetworkVariable.h"
#include "model/variables/EffectValueTable.h"

namespace siena
{

// ----------------------------------------------------------------------------
// Section: Helpers
// ----------------------------------------------------------------------------

/**
 * An iterator type definition for convenient use in several places.
 */
typedef std::map<const NetworkVariable *,
	EffectValueTable *>::const_iterator Iterator;


/**
 * An identity function simply returning its argument.
 */
inline double identity(int x)
{
	return x;
}


/**
 * Returns the inverse of (<i>x</i> + 1).
 */
inline double invertor(int x)
{
	return 1.0 / (x + 1);
}


// ----------------------------------------------------------------------------
// Section: Construction area
// ----------------------------------------------------------------------------

/**
 * Creates a dependent variable for the given set of actors and the given
 * number of observations.
 * @param pSimulation the model simulation, which becomes the owner of
 * this variable
 */
DependentVariable::DependentVariable(string name,
	const ActorSet * pActorSet,
	int observationCount,
	EpochSimulation * pSimulation) : NamedObject(name)
{
	this->lpActorSet = pActorSet;
	this->lpSimulation = pSimulation;
	this->ltotalRate = 0;
	this->lrate = new double[this->n()];
	this->lcovariateRates = new double[this->n()];
	this->lpEvaluationFunction = new Function();
	this->lpEndowmentFunction = new Function();
}


/**
 * Reads the parameters of rate effects from the model and stores them
 * internally. This method must be called right after the creation of
 * all dependent variables.
 */
void DependentVariable::initializeRateFunction()
{
	const Data * pData = this->lpSimulation->pData();
	const Model * pModel = this->lpSimulation->pModel();

	const vector<EffectInfo *> & rRateEffects =
		pModel->rRateEffects(this->name());

	for (unsigned i = 0; i < rRateEffects.size(); i++)
	{
		EffectInfo * pEffectInfo = rRateEffects[i];
		double parameter = pEffectInfo->parameter();
		string effectName = pEffectInfo->effectName();
		string interactionName = pEffectInfo->interactionName1();
		string rateType = pEffectInfo->rateType();

		if (rateType == "covariate")
		{
			//	Rprintf("covariate\n");
			// Covariate-dependent rate effect

			//	if (parameter != 0)
			//	{
				ConstantCovariate * pConstantCovariate =
					pData->pConstantCovariate(interactionName);
				ChangingCovariate * pChangingCovariate =
					pData->pChangingCovariate(interactionName);
				const BehaviorVariable * pBehaviorVariable =
					(const BehaviorVariable *)
						this->lpSimulation->pVariable(interactionName);

				if (pConstantCovariate)
				{
					if (this->lpActorSet != pConstantCovariate->pActorSet())
					{
						throw domain_error("Mismatch of actor sets");
					}

					this->lconstantCovariateParameters[pConstantCovariate] =
						parameter;
					this->lconstantCovariateScores[pConstantCovariate] = 0;
				}
				else if (pChangingCovariate)
				{
					if (this->lpActorSet != pChangingCovariate->pActorSet())
					{
						throw domain_error("Mismatch of actor sets");
					}

					this->lchangingCovariateParameters[pChangingCovariate] =
						parameter;
					this->lchangingCovariateScores[pChangingCovariate] = 0;
				}
				else if (pBehaviorVariable)
				{
					if (this->lpActorSet != pBehaviorVariable->pActorSet())
					{
						throw domain_error("Mismatch of actor sets");
					}

					this->lbehaviorVariableParameters[pBehaviorVariable] =
						parameter;
					this->lbehaviorVariableScores[pBehaviorVariable] = 0;
				}
				else
					throw logic_error(
						"No individual covariate named '" +
						interactionName +
						"'.");
				//	}
		}
		else
		{
			// We expect a structural (network-dependent) rate effect here.

			const NetworkVariable * pVariable;

			if (interactionName == "")
			{
				 pVariable = dynamic_cast<const NetworkVariable *>(this);
			}
			else
			{
				 pVariable = dynamic_cast<const NetworkVariable *>(
					this->lpSimulation->pVariable(interactionName));
			}

			if (!pVariable)
			{
				throw logic_error("The structural rate effect " +
					effectName +
					" for dependent variable " +
					this->name() +
					" refers to a non-existing network variable " +
					interactionName);
			}

			if (effectName == "outRate")
			{
				this->outDegreeRateParameter(pVariable, parameter);
				this->loutDegreeScores[pVariable] = 0;
			}
			else if (effectName == "inRate")
			{
				this->inDegreeRateParameter(pVariable, parameter);
				this->linDegreeScores[pVariable] = 0;
			}
			else if (effectName == "recipRate")
			{
				this->reciprocalDegreeRateParameter(pVariable, parameter);
				this->lreciprocalDegreeScores[pVariable] = 0;
			}
			else if (effectName == "outRateInv")
			{
				this->inverseOutDegreeRateParameter(pVariable, parameter);
				this->linverseOutDegreeScores[pVariable] = 0;
			}
			else
			{
				throw domain_error("Unexpected rate effect " + effectName);
			}
		}
	}

	// If there are no rate effects depending on changing covariates,
    // or behavior variables then the covariate based rates can be calculated
	// just once.

	if (this->lchangingCovariateParameters.empty() &&
		this->lbehaviorVariableParameters.empty())
	{
		this->updateCovariateRates();
	}
}


void DependentVariable::initializeEvaluationFunction()
{
	this->initializeFunction(this->lpEvaluationFunction,
		this->lpSimulation->pModel()->rEvaluationEffects(this->name()));
}


void DependentVariable::initializeEndowmentFunction()
{
	this->initializeFunction(this->lpEndowmentFunction,
		this->lpSimulation->pModel()->rEndowmentEffects(this->name()));
}


void DependentVariable::initializeFunction(Function * pFunction,
	const vector<EffectInfo *> & rEffects) const
{
	EffectFactory factory;

	for (unsigned i = 0; i < rEffects.size(); i++)
	{
		EffectInfo * pEffectInfo = rEffects[i];
		Effect * pEffect = factory.createEffect(pEffectInfo);
		pFunction->addEffect(pEffect);
		pEffect->initialize(this->lpSimulation);
	}
}


/**
 * Deallocates this dependent variable.
 */
DependentVariable::~DependentVariable()
{
	delete this->lpEvaluationFunction;
	delete this->lpEndowmentFunction;
	delete[] this->lrate;
	delete[] this->lcovariateRates;

	// Delete the structural rate effect tables.

	clearMap(this->loutDegreeRateEffects, false, true);
	clearMap(this->linDegreeRateEffects, false, true);
	clearMap(this->lreciprocalDegreeRateEffects, false, true);
	clearMap(this->linverseOutDegreeRateEffects, false, true);

	// Nullify the fields

	this->lpSimulation = 0;
	this->lrate = 0;
	this->lcovariateRates = 0;
	this->lpEvaluationFunction = 0;
	this->lpEndowmentFunction = 0;
}


// ----------------------------------------------------------------------------
// Section: Accessors
// ----------------------------------------------------------------------------

/**
 * Returns the set of actors underlying this dependent variable.
 */
const ActorSet * DependentVariable::pActorSet() const
{
	return this->lpActorSet;
}


/**
 * Returns the number of actors for this variable.
 */
int DependentVariable::n() const
{
	return this->lpActorSet->n();
}


/**
 * Stores the distance of this variable to the observed data at the
 * beginning of the current period.
 */
void DependentVariable::simulatedDistance(int distance)
{
	this->lsimulatedDistance = distance;
}


/**
 * Returns the distance of this variable to the observed data at the
 * beginning of the current period.
 */
int DependentVariable::simulatedDistance() const
{
	return this->lsimulatedDistance;
}

/**
 * sets the changing covariate period.
 */
void DependentVariable::changingCovariatePeriod(int changingCovariatePeriod)
{
	this->lchangingCovariatePeriod = changingCovariatePeriod;
}

/**
 * Returns the changing covariate period
 */
int DependentVariable::changingCovariatePeriod() const
{
	return this->lchangingCovariatePeriod;
}

// ----------------------------------------------------------------------------
// Section: Initialization
// ----------------------------------------------------------------------------

/**
 * Initializes this dependent variable as of the beginning of the given period.
 */
void DependentVariable::initialize(int period)
{
	this->lperiod = period;
	this->lsimulatedDistance = 0;
	this->lbasicRateScore = 0;
	this->lbasicRate =
		this->lpSimulation->pModel()->basicRateParameter(this->pData(),
			period);

	if (!this->lchangingCovariateParameters.empty() ||
		!this->lbehaviorVariableParameters.empty())
	{
		// The changing covariates may have different values in different
		// periods, hence the covariate-based rates have to be recalculated.

		// TODO: Behavior variables change all the time, so I'm not sure
		// it is enough to update the covariate rates at the start of the
		// period only.

		this->updateCovariateRates();
	}
}


// ----------------------------------------------------------------------------
// Section: Rate calculation
// ----------------------------------------------------------------------------

/**
 * Calculates the rate of change or each actor and the total rate.
 */
void DependentVariable::calculateRates()
{
	this->ltotalRate = 0;

	for (int i = 0; i < this->n(); i++)
	{
		// If an actor cannot make a change with respect to this variable,
		// then its rate is 0.

		if (this->canMakeChange(i))
		{
			this->lrate[i] = this->calculateRate(i);
		}
		else
		{
			this->lrate[i] = 0;
		}

		this->ltotalRate += this->lrate[i];
	}
}


/**
 * Returns if the given actor can change the current state of this variable.
 */
bool DependentVariable::canMakeChange(int actor) const
{
	return this->lpSimulation->active(this->pActorSet(), actor);
}


/**
 * Calculates the rate of the given actor.
 */
double DependentVariable::calculateRate(int i)
{
	// The rate is the product of the basic rate parameter for the current
	// period, exponentials of some covariate-based effects, and exponentials
	// of some effects depending on the structure of certain networks. The
	// later two components are precomputed for efficiency.

	return this->basicRate() *
		this->lcovariateRates[i] *
		this->structuralRate(i);
}


/**
 * Returns the total rate of change over all actors.
 */
double DependentVariable::totalRate() const
{
	return this->ltotalRate;
}


/**
 * Returns the rate of change for the given actor. It is assumed that the
 * rates have been calculated already by calling the method calculateRates.
 */
double DependentVariable::rate(int actor) const
{
	return this->lrate[actor];
}


/**
 * Recalculates the covariate-based components of the rate functions using
 * the current values of parameters and changing covariates.
 */
void DependentVariable::updateCovariateRates()
{
	// Nullify the array.

	for (int i = 0; i < this->n(); i++)
	{
		this->lcovariateRates[i] = 0;
	}

	// Add the contributions of each constant covariate with non-zero parameter
	// for the rate functions. The contribution of a constant covariate v to
	// the rate function of an actor i is exp(\alpha v[i]), where \alpha is
	// the corresponding parameter. The calculation of the exponentials is
	// postponed, though.

	for (std::map<const ConstantCovariate *, double>::iterator iter =
			this->lconstantCovariateParameters.begin();
		iter != this->lconstantCovariateParameters.end();
		iter++)
	{
		const ConstantCovariate * pCovariate = iter->first;
		double parameter = iter->second;

		for (int i = 0; i < this->n(); i++)
		{
			this->lcovariateRates[i] += parameter * pCovariate->value(i);
		}
	}

	// Add the contributions of each changing covariate with non-zero parameter
	// for the rate functions. The contribution of a changing covariate v to
	// the rate function of an actor i is exp(\alpha v[i][h]), where \alpha is
	// the corresponding parameter, and h is the current period.
	// Again, the calculation of the exponentials is postponed.

	for (std::map<const ChangingCovariate *, double>::iterator iter =
			this->lchangingCovariateParameters.begin();
		iter != this->lchangingCovariateParameters.end();
		iter++)
	{
		const ChangingCovariate * pCovariate = iter->first;
		double parameter = iter->second;
		for (int i = 0; i < this->n(); i++)
		{
			this->lcovariateRates[i] +=
				parameter * pCovariate->value(i, this->period());
		}
	}

	// Add the contributions of each behavior variable with non-zero parameter
	// for the rate functions. The contribution of a behavior variable v to
	// the rate function of an actor i is exp(\alpha v[i]), where \alpha is
	// the corresponding parameter.
	// Again, the calculation of the exponentials is postponed.

	for (std::map<const BehaviorVariable *, double>::iterator iter =
			this->lbehaviorVariableParameters.begin();
		iter != this->lbehaviorVariableParameters.end();
		iter++)
	{
		const BehaviorVariable * pBehavior = iter->first;
		double parameter = iter->second;

		for (int i = 0; i < this->n(); i++)
		{
			this->lcovariateRates[i] += parameter * pBehavior->value(i);
		}
	}

	// Okay, now we take the exponentials of the sums of contributions over
	// all covariates. This is valid, because exp(x_1 + ... + x_k) =
	// exp(x_1) * ... * exp(x_k).

	for (int i = 0; i < this->n(); i++)
	{
		this->lcovariateRates[i] = exp(this->lcovariateRates[i]);
	}
}


/**
 * Stores the parameter value for the structural rate effect depending on
 * the out-degrees in the given one-mode network variable.
 */
void DependentVariable::outDegreeRateParameter(
	const NetworkVariable * pVariable,
	double value)
{
	if (this->lpActorSet != pVariable->pSenders())
	{
		throw std::invalid_argument("Mismatch of actor sets");
	}

	EffectValueTable * pTable = this->loutDegreeRateEffects[pVariable];

	if (!pTable)
	{
		pTable = new EffectValueTable(pVariable->m(), identity);
		this->loutDegreeRateEffects[pVariable] = pTable;
	}

	pTable->parameter(value);
}


/**
 * Stores the parameter value for the structural rate effect depending on
 * the in-degrees in the given one-mode network variable.
 */
void DependentVariable::inDegreeRateParameter(
	const NetworkVariable * pVariable,
	double value)
{
	if (this->lpActorSet != pVariable->pReceivers())
	{
		throw std::invalid_argument("Mismatch of actor sets");
	}

	EffectValueTable * pTable = this->linDegreeRateEffects[pVariable];

	if (!pTable)
	{
		pTable = new EffectValueTable(pVariable->n(), identity);
		this->linDegreeRateEffects[pVariable] = pTable;
	}

	pTable->parameter(value);
}


/**
 * Stores the parameter value for the structural rate effect depending on
 * the reciprocal degrees in the given one-mode network variable.
 */
void DependentVariable::reciprocalDegreeRateParameter(
	const NetworkVariable * pVariable,
	double value)
{
	if (!pVariable->oneModeNetwork())
	{
		throw std::invalid_argument("One-mode network variable expected");
	}

	if (this->lpActorSet != pVariable->pSenders())
	{
		throw std::invalid_argument("Mismatch of actor sets");
	}

	EffectValueTable * pTable = this->lreciprocalDegreeRateEffects[pVariable];

	if (!pTable)
	{
		pTable = new EffectValueTable(pVariable->n(), identity);
		this->lreciprocalDegreeRateEffects[pVariable] = pTable;
	}

	pTable->parameter(value);
}


/**
 * Stores the parameter value for the structural rate effect depending on
 * the inverse out-degrees in the given one-mode network variable.
 */
void DependentVariable::inverseOutDegreeRateParameter(
	const NetworkVariable * pVariable,
	double value)
{
	if (this->lpActorSet != pVariable->pSenders())
	{
		throw std::invalid_argument("Mismatch of actor sets");
	}

	EffectValueTable * pTable = this->linverseOutDegreeRateEffects[pVariable];

	if (!pTable)
	{
		pTable = new EffectValueTable(pVariable->m(), invertor);
		this->linverseOutDegreeRateEffects[pVariable] = pTable;
	}

	pTable->parameter(value);
}


/**
 * Returns the component of the rate function of actor <i>i</i> depending
 * on structural effects.
 */
double DependentVariable::structuralRate(int i) const
{
	double rate = 1;

	// Calculate the rate for outdegree effects.

	for (Iterator iter = this->loutDegreeRateEffects.begin();
		iter != this->loutDegreeRateEffects.end();
		iter++)
	{
		// The network variable
		const NetworkVariable * pVariable = iter->first;

		// The current state of the network variable
		Network * pNetwork = pVariable->pNetwork();

		// The table of precalculated contributions for each outdegree
		EffectValueTable * pTable = iter->second;

		// Update the rate
		rate *= pTable->value(pNetwork->outDegree(i));
	}

	// Other structural rate effects are treated similarly.

	for (Iterator iter = this->linDegreeRateEffects.begin();
		iter != this->linDegreeRateEffects.end();
		iter++)
	{
		const NetworkVariable * pVariable = iter->first;
		Network * pNetwork = pVariable->pNetwork();
		EffectValueTable * pTable = iter->second;

		rate *= pTable->value(pNetwork->inDegree(i));
	}

	for (Iterator iter = this->lreciprocalDegreeRateEffects.begin();
		iter != this->lreciprocalDegreeRateEffects.end();
		iter++)
	{
		const NetworkVariable * pVariable = iter->first;
		OneModeNetwork * pNetwork = (OneModeNetwork *) pVariable->pNetwork();
		EffectValueTable * pTable = iter->second;

		rate *= pTable->value(pNetwork->reciprocalDegree(i));
	}

	for (Iterator iter = this->linverseOutDegreeRateEffects.begin();
		iter != this->linverseOutDegreeRateEffects.end();
		iter++)
	{
		const NetworkVariable * pVariable = iter->first;
		Network * pNetwork = pVariable->pNetwork();
		EffectValueTable * pTable = iter->second;

		rate *= pTable->value(pNetwork->outDegree(i));
	}

	return rate;
}


/**
 * Updates the rate score functions for this event for this variable.
 * @param[in] tau the time increment in the current step of the simulation
 * @param[in] pSelectedVariable the variable, which has been selected in
 * the current step of the simulation (0, if none is selected)
 * @param[in] selectedActor the actor, which has been selected to change
 * the selected variable, if any. If no variable has been selected, this
 * parameter is ignored.
 */
void DependentVariable::accumulateRateScores(double tau,
	const DependentVariable * pSelectedVariable,
	int selectedActor)
{
	// Update the score for the basic rate parameter

	if (this == pSelectedVariable)
	{
		this->lbasicRateScore += 1.0 / this->basicRate();
	}

	this->lbasicRateScore -= this->totalRate() * tau / this->basicRate();

	// Update scores for covariate dependent rate parameters

	for (std::map<const ConstantCovariate *, double>::iterator iter =
			this->lconstantCovariateScores.begin();
		iter != this->lconstantCovariateScores.end();
		iter++)
	{
		const ConstantCovariate * pCovariate = iter->first;

		if (this == pSelectedVariable)
		{
			iter->second += pCovariate->value(selectedActor);
		}

		for (int i = 0; i < this->n(); i++)
		{
			iter->second -= this->lrate[i] * pCovariate->value(i) * tau;
		}
	}

	for (std::map<const ChangingCovariate *, double>::iterator iter =
			this->lchangingCovariateScores.begin();
		iter != this->lchangingCovariateScores.end();
		iter++)
	{
		const ChangingCovariate * pCovariate = iter->first;
		if (this == pSelectedVariable)
		{
			iter->second += pCovariate->value(selectedActor, this->period());
		}

		for (int i = 0; i < this->n(); i++)
		{
			iter->second -=
				this->lrate[i] * pCovariate->value(i, this->period()) * tau;
		}
	}

	for (std::map<const BehaviorVariable *, double>::iterator iter =
			this->lbehaviorVariableScores.begin();
		iter != this->lbehaviorVariableScores.end();
		iter++)
	{
		const BehaviorVariable * pBehavior = iter->first;

		if (this == pSelectedVariable)
		{
			iter->second += pBehavior->value(selectedActor);
		}

		for (int i = 0; i < this->n(); i++)
		{
			iter->second -= this->lrate[i] * pBehavior->value(i) * tau;
		}
	}

	// Update scores for structural rate parameters

	for (std::map<const NetworkVariable *, double>::iterator iter =
			this->loutDegreeScores.begin();
		iter != this->loutDegreeScores.end();
		iter++)
	{
		const Network * pNetwork = iter->first->pNetwork();

		if (this == pSelectedVariable)
		{
			iter->second += pNetwork->outDegree(selectedActor);
		}

		for (int i = 0; i < this->n(); i++)
		{
			iter->second -= this->rate(i) *	pNetwork->outDegree(i) * tau;
			//		Rprintf("%d %f %d %f\n",i,this->rate(i), pNetwork->outDegree(i), tau);
		}
		//	Rprintf("%f 1\n", iter->second);
	}

	for (std::map<const NetworkVariable *, double>::iterator iter =
			this->linDegreeScores.begin();
		iter != this->linDegreeScores.end();
		iter++)
	{
		const Network * pNetwork = iter->first->pNetwork();

		if (this == pSelectedVariable)
		{
			iter->second += pNetwork->inDegree(selectedActor);
		}

		for (int i = 0; i < this->n(); i++)
		{
			iter->second -= this->rate(i) *	pNetwork->inDegree(i) * tau;
		}
		//	Rprintf("%f 1\n", iter->second);
	}

	for (std::map<const NetworkVariable *, double>::iterator iter =
			this->lreciprocalDegreeScores.begin();
		iter != this->lreciprocalDegreeScores.end();
		iter++)
	{
		const OneModeNetwork * pNetwork =
			(const OneModeNetwork *) iter->first->pNetwork();

		if (this == pSelectedVariable)
		{
			iter->second += pNetwork->reciprocalDegree(selectedActor);
		}

		for (int i = 0; i < this->n(); i++)
		{
			iter->second -=
				this->rate(i) *	pNetwork->reciprocalDegree(i) * tau;
		}
	}

	for (std::map<const NetworkVariable *, double>::iterator iter =
			this->linverseOutDegreeScores.begin();
		iter != this->linverseOutDegreeScores.end();
		iter++)
	{
		const Network * pNetwork = iter->first->pNetwork();

		if (this == pSelectedVariable)
		{
			iter->second += invertor(pNetwork->outDegree(selectedActor));
		}

		for (int i = 0; i < this->n(); i++)
		{
			iter->second -=
				this->rate(i) * invertor(pNetwork->outDegree(i)) * tau;
		}
	}
}


double DependentVariable::basicRateScore() const
{
	return this->lbasicRateScore;
}


double DependentVariable::constantCovariateScore(
	const ConstantCovariate * pCovariate) const
{
	map<const ConstantCovariate *, double>::const_iterator iter =
		this->lconstantCovariateScores.find(pCovariate);

	if (iter == this->lconstantCovariateScores.end())
	{
		throw invalid_argument(
			string("Unknown covariate: The given covariate rate ") +
			string("effect is not part of the model."));
	}

	return iter->second;
}


double DependentVariable::changingCovariateScore(
	const ChangingCovariate * pCovariate) const
{
	map<const ChangingCovariate *, double>::const_iterator iter =
		this->lchangingCovariateScores.find(pCovariate);
	if (iter == this->lchangingCovariateScores.end())
	{
		throw invalid_argument(
			string("Unknown covariate: The given covariate rate ") +
			string("effect is not part of the model."));
	}

	return iter->second;
}


double DependentVariable::behaviorVariableScore(
	const BehaviorVariable * pBehavior) const
{
	map<const BehaviorVariable *, double>::const_iterator iter =
		this->lbehaviorVariableScores.find(pBehavior);

	if (iter == this->lbehaviorVariableScores.end())
	{
		throw invalid_argument(
			string("Unknown behavior variable: ") +
			"The given covariate rate effect is not part of the model.");
	}

	return iter->second;

}


double DependentVariable::outDegreeScore(
	const NetworkVariable * pNetworkData) const
{
	map<const NetworkVariable *, double>::const_iterator iter =
		this->loutDegreeScores.find(pNetworkData);

	if (iter == this->loutDegreeScores.end())
	{
		throw invalid_argument(
			string("Unknown network: ") +
			"The given outdegree rate effect is not part of the model.");
	}

	return iter->second;

}


double DependentVariable::inDegreeScore(
	const NetworkVariable * pNetworkData) const
{
	map<const NetworkVariable *, double>::const_iterator iter =
		this->linDegreeScores.find(pNetworkData);

	if (iter == this->linDegreeScores.end())
	{
		throw invalid_argument(
			string("Unknown network: ") +
			"The given indegree rate effect is not part of the model.");
	}

	return iter->second;
}


double DependentVariable::reciprocalDegreeScore(
	const NetworkVariable * pNetworkData) const
{
	map<const NetworkVariable *, double>::const_iterator iter =
		this->lreciprocalDegreeScores.find(pNetworkData);

	if (iter == this->lreciprocalDegreeScores.end())
	{
		throw invalid_argument(
			string("Unknown network: ") +
			"The given reciprocal degree rate effect is not " +
			"part of the model.");
	}

	return iter->second;
}


double DependentVariable::inverseOutDegreeScore(
	const NetworkVariable * pNetworkData) const
{
	map<const NetworkVariable *, double>::const_iterator iter =
		this->linverseOutDegreeScores.find(pNetworkData);

	if (iter == this->linverseOutDegreeScores.end())
	{
		throw invalid_argument(
			string("Unknown network: ") +
			"The given inverse outdegree rate effect is not " +
			"part of the model.");
	}

	return iter->second;
}

}
