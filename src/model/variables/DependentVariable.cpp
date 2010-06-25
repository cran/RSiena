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

#include "BehaviorVariable.h"
#include "DependentVariable.h"
#include "utils/Utils.h"
#include "utils/Random.h"
#include "data/ActorSet.h"
#include "data/ConstantCovariate.h"
#include "data/ChangingCovariate.h"
#include "network/Network.h"
#include "network/OneModeNetwork.h"
#include "data/LongitudinalData.h"
#include "data/NetworkLongitudinalData.h"
#include "model/Model.h"
#include "model/EffectInfo.h"
#include "model/EpochSimulation.h"
#include "model/ml/MLSimulation.h"
#include "model/SimulationActorSet.h"
#include "model/effects/AllEffects.h"
#include "model/effects/EffectFactory.h"
#include "model/effects/StructuralRateEffect.h"
#include "model/variables/NetworkVariable.h"
#include "model/variables/EffectValueTable.h"

namespace siena
{

// ----------------------------------------------------------------------------
// Section: Construction area
// ----------------------------------------------------------------------------

/**
 * Creates a dependent variable for the given set of actors.
 * @param pSimulation the model simulation, which becomes the owner of
 * this variable
 */
DependentVariable::DependentVariable(string name,
	const ActorSet * pActorSet,
	EpochSimulation * pSimulation) : NamedObject(name)
{
	this->lpActorSet = pSimulation->pSimulationActorSet(pActorSet);
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
					if (this->lpActorSet->pOriginalActorSet() !=
						pConstantCovariate->pActorSet())
					{
						throw domain_error("Mismatch of actor sets");
					}

					this->lconstantCovariateParameters[pConstantCovariate] =
						parameter;
					this->lconstantCovariateScores[pConstantCovariate] = 0;
				}
				else if (pChangingCovariate)
				{
					if (this->lpActorSet->pOriginalActorSet() !=
						pChangingCovariate->pActorSet())
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
				if (this->lpActorSet != pVariable->pSenders())
				{
					throw std::invalid_argument("Mismatch of actor sets");
				}

				this->lstructuralRateEffects.push_back(
					new StructuralRateEffect(pVariable,
						OUT_DEGREE_RATE,
						parameter));
				this->loutDegreeScores[pVariable] = 0;
			}
			else if (effectName == "inRate")
			{
				if (this->lpActorSet != pVariable->pReceivers())
				{
					throw std::invalid_argument("Mismatch of actor sets");
				}

				this->lstructuralRateEffects.push_back(
					new StructuralRateEffect(pVariable,
						IN_DEGREE_RATE,
						parameter));
				this->linDegreeScores[pVariable] = 0;
			}
			else if (effectName == "recipRate")
			{
				if (!pVariable->oneModeNetwork())
				{
					throw std::invalid_argument(
						"One-mode network variable expected");
				}

				if (this->lpActorSet != pVariable->pSenders())
				{
					throw std::invalid_argument("Mismatch of actor sets");
				}

				this->lstructuralRateEffects.push_back(
					new StructuralRateEffect(pVariable,
						RECIPROCAL_DEGREE_RATE,
						parameter));
				this->lreciprocalDegreeScores[pVariable] = 0;
			}
			else if (effectName == "outRateInv")
			{
				if (this->lpActorSet != pVariable->pSenders())
				{
					throw std::invalid_argument("Mismatch of actor sets");
				}

				this->lstructuralRateEffects.push_back(
					new StructuralRateEffect(pVariable,
						INVERSE_OUT_DEGREE_RATE,
						parameter));
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
	EffectFactory factory(this->pSimulation()->pData());

	for (unsigned i = 0; i < rEffects.size(); i++)
	{
		EffectInfo * pEffectInfo = rEffects[i];
		Effect * pEffect = factory.createEffect(pEffectInfo);
		pFunction->addEffect(pEffect);
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

	// Delete the structural rate effects.
	deallocateVector(this->lstructuralRateEffects);

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
 * Returns the number of actors for this variable.
 */
int DependentVariable::n() const
{
	return this->lpActorSet->n();
}


/**
 * Returns the ID of this dependent variable, which is the same as the ID
 * of the underlying observed data object.
 */
int DependentVariable::id() const
{
	return this->pData()->id();
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
	this->lbasicRateDerivative = 0;
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

	// Make sure the rates will be calculated whenever calculateRates()
	// is called next time.

	this->lvalidRates = false;
}


// ----------------------------------------------------------------------------
// Section: Rate calculation
// ----------------------------------------------------------------------------

/**
 * Calculates the rate of change or each actor and the total rate.
 */
void DependentVariable::calculateRates()
{
	if (!this->constantRates() || !this->lvalidRates)
	{
		this->ltotalRate = 0;
		int n = this->n();

		for (int i = 0; i < n; i++)
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

		this->lvalidRates = true;
	}
}


/**
 * Ensures that the rates will be recalculated when the method
 * calculateRates() is called next time.
 */
void DependentVariable::invalidateRates()
{
	this->lvalidRates = false;
}


/**
 * Returns if the given actor can change the current state of this variable.
 */
bool DependentVariable::canMakeChange(int actor) const
{
	return this->lpActorSet->active(actor);
}


/**
 * Returns if the rates are constant over the current period, and can be
 * thus calculated just once.
 */
bool DependentVariable::constantRates() const
{
	return this->lstructuralRateEffects.empty() &&
		this->lbehaviorVariableParameters.empty();
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
 * Returns the component of the rate function of actor <i>i</i> depending
 * on structural effects.
 */
double DependentVariable::structuralRate(int i) const
{
	double rate = 1;
	int effectCount = this->lstructuralRateEffects.size();

	for (int effectIndex = 0; effectIndex < effectCount; effectIndex++)
	{
		rate *= this->lstructuralRateEffects[effectIndex]->value(i);
	}

	return rate;
}


// ----------------------------------------------------------------------------
// Section: Scores
// ----------------------------------------------------------------------------

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

/**
 * Calculates the rate score functions for this chain for this variable.
 * @param[in] activeMiniStepCount the number of non-structurally determined
 * links in the current chain for this variable.
 */

void DependentVariable::calculateMaximumLikelihoodRateScores(int
	activeMiniStepCount)
{
//	Rprintf("%d %d %f\n", this->n(), activeMiniStepCount, this->basicRate());
	this->lbasicRateScore =
		- this->n() + activeMiniStepCount / this->basicRate();
}

/**
 * Calculates the rate derivative functions for this chain for this variable.
 * @param[in] activeMiniStepCount the number of non-structurally determined
 * links in the current chain for this variable.
 */

void DependentVariable::calculateMaximumLikelihoodRateDerivatives(int
	activeMiniStepCount)
{
//	Rprintf("%d %d %f\n", this->n(), activeMiniStepCount, this->basicRate());
	this->lbasicRateDerivative =
		- activeMiniStepCount / this->basicRate()/ this->basicRate();
}

double DependentVariable::basicRateScore() const
{
	return this->lbasicRateScore;
}

double DependentVariable::basicRateDerivative() const
{
	return this->lbasicRateDerivative;
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


// ----------------------------------------------------------------------------
// Section: Composition change
// ----------------------------------------------------------------------------

/**
 * Updates this variable when an actor becomes active.
 */
void DependentVariable::actOnJoiner(const SimulationActorSet * pActorSet,
	int actor)
{
	if (pActorSet == this->lpActorSet)
	{
		this->invalidateRates();
	}
}


/**
 * Updates this variable when an actor becomes inactive.
 */
void DependentVariable::actOnLeaver(const SimulationActorSet * pActorSet,
	int actor)
{
	if (pActorSet == this->lpActorSet)
	{
		this->invalidateRates();
	}
}


// ----------------------------------------------------------------------------
// Section: Maximum likelihood related methods
// ----------------------------------------------------------------------------

/**
 * Returns whether applying the given ministep on the current state of this
 * variable would be valid with respect to all constraints.
 */
bool DependentVariable::validMiniStep(const MiniStep * pMiniStep) const
{
	return true;
}


/**
 * Updates effect parameters
 */

void DependentVariable::updateEffectParameters()
{
	// find the Evaluation effectInfos
	const vector<EffectInfo *>  rEffects=
		this->lpSimulation->pModel()->rEvaluationEffects(this->name());

	const Function * pFunction = this->pEvaluationFunction();

	for (unsigned i = 0; i < pFunction->rEffects().size(); i++)
	{
		Effect * pEffect = pFunction->rEffects()[i];
		pEffect->parameter(rEffects[i]->parameter());
	}
	// find the Endowment effectInfos
	const vector<EffectInfo *>  rEffects2=
	 this->lpSimulation->pModel()->rEndowmentEffects(this->name());

	pFunction = this->pEndowmentFunction();

	for (unsigned i = 0; i < pFunction->rEffects().size(); i++)
	{
		Effect * pEffect = pFunction->rEffects()[i];
		pEffect->parameter(rEffects2[i]->parameter());
	}
}
/**
 * Updates effect info parameters from the effects (used if accepted)
 */

void DependentVariable::updateEffectInfoParameters()
{
	// find the Evaluation effectInfos
	const vector<EffectInfo *>  rEffects=
		this->lpSimulation->pModel()->rEvaluationEffects(this->name());

	const Function * pFunction = this->pEvaluationFunction();

	for (unsigned i = 0; i < pFunction->rEffects().size(); i++)
	{
		Effect * pEffect = pFunction->rEffects()[i];
		rEffects[i]->parameter(pEffect->parameter());
	}
	// find the Endowment effectInfos
	const vector<EffectInfo *>  rEffects2=
	 this->lpSimulation->pModel()->rEndowmentEffects(this->name());

	pFunction = this->pEndowmentFunction();

	for (unsigned i = 0; i < pFunction->rEffects().size(); i++)
	{
		Effect * pEffect =	pFunction->rEffects()[i];
		rEffects2[i]->parameter(pEffect->parameter());
	}
}

/**
 * Samples from distribution for basic rate parameters
 *
 */
void DependentVariable::sampleBasicRate(int miniStepCount)
{
	this->lbasicRate = nextGamma(miniStepCount + 1, 1.0 / this->n());
	MLSimulation * pMLSimulation =
		dynamic_cast<MLSimulation * >(this->pSimulation());
	pMLSimulation->sampledBasicRates(this->lbasicRate);
	pMLSimulation->sampledBasicRatesDistributions(miniStepCount + 1);
	this->lvalidRates = false;
}
/**
 * Samples from distribution for non basic rate parameters
 *
 */
double DependentVariable::sampleParameters(double scaleFactor)
{
	double priorSD = 10000;
	double priornew = 0;
	double priorold = 0;

	// need to access this later
	MLSimulation * pMLSimulation =
		dynamic_cast<MLSimulation * >(this->pSimulation());

	// create candidate set of parameters and copy them to the effects

	const Function * pFunction = this->pEvaluationFunction();

	for (unsigned i = 0; i < pFunction->rEffects().size(); i++)
	{
		Effect * pEffect = pFunction->rEffects()[i];
		double candidate = pEffect->parameter() + nextNormal(0, scaleFactor);
		priornew += normalDensity(candidate, 0, priorSD, true);

		pMLSimulation->candidates(pEffect->pEffectInfo(), candidate);
		priorold += normalDensity(pEffect->parameter(), 0, priorSD, true);
		pEffect->parameter(candidate);

	}

	pFunction = this->pEndowmentFunction();

	for (unsigned i = 0; i < pFunction->rEffects().size(); i++)
	{
		Effect * pEffect = pFunction->rEffects()[i];
		double candidate = pEffect->parameter() + nextNormal(0, scaleFactor);
		priornew += normalDensity(candidate, 0, priorSD, true);
		pMLSimulation->candidates(pEffect->pEffectInfo(), candidate);
		priorold += normalDensity(pEffect->parameter(), 0, priorSD,	true);
		pEffect->parameter(candidate);
	}
	return priornew - priorold;
}

// ----------------------------------------------------------------------------
// Section: Properties
// ----------------------------------------------------------------------------

/**
 * Returns if this is a network variable.
 */
bool DependentVariable::networkVariable() const
{
	// This method is overriden in NetworkVariable. Here we return false.
	return false;
}


/**
 * Returns if this is a behavior variable.
 */
bool DependentVariable::behaviorVariable() const
{
	// This method is overriden in BehaviorVariable. Here we return false.
	return false;
}

}
