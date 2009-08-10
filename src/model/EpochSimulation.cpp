/******************************************************************************
 * SIENA: Simulation Investigation for Empirical Network Analysis
 *
 * Web: http://www.stats.ox.ac.uk/~snijders/siena/
 *
 * File: EpochSimulation.cpp
 *
 * Description: This file contains the implementation of the
 * EpochSimulation class.
 *****************************************************************************/

#include <algorithm>
#include <R.h>

#include "EpochSimulation.h"
#include "utils/Random.h"
#include "utils/Utils.h"
#include "data/ActorSet.h"
#include "data/ExogenousEvent.h"
#include "data/LongitudinalData.h"
#include "model/variables/DependentVariable.h"
#include "model/variables/NetworkVariable.h"
#include "model/variables/BehaviorVariable.h"
#include "model/Model.h"
#include "model/SimulationActorSet.h"
#include "model/effects/Effect.h"

namespace siena
{

// ----------------------------------------------------------------------------
// Section: Constructors and destructors
// ----------------------------------------------------------------------------

/**
 * Creates an epoch simulation object for the given observed data and an
 * actor-based model for that data.
 */
EpochSimulation::EpochSimulation(Data * pData, Model * pModel)
{
    this->lpData = pData;
    this->lpModel = pModel;
    this->lpConditioningVariable = 0;

	// Create a wrapper for each actor set for simulation purposes,
	// and find the maximum number of actors in any actor set.

	int maxN = 0;

	for (unsigned i = 0; i < pData->rActorSets().size(); i++)
	{
		const ActorSet * pActorSet = pData->rActorSets()[i];
		SimulationActorSet * pSimulationActorSet =
			new SimulationActorSet(pActorSet);

		this->lsimulationActorSets.push_back(pSimulationActorSet);
		this->lactorSetMap[pActorSet] = pSimulationActorSet;
		maxN = std::max(maxN, pActorSet->n());
	}

    // Create the dependent variables from the observed data

    for (unsigned i = 0; i < pData->rDependentVariableData().size(); i++)
    {
    	DependentVariable * pVariable =
    		pData->rDependentVariableData()[i]->createVariable(this);
        this->lvariables.push_back(pVariable);

        if (pModel->conditional() &&
        	pModel->conditionalDependentVariable() == pVariable->name())
        {
        	this->lpConditioningVariable = pVariable;
        }
    }

    // Initialize the rate, evaluation, and endowment
    // functions of all variables.

    for (unsigned i = 0; i < this->lvariables.size(); i++)
    {
    	this->lvariables[i]->initializeRateFunction();
    	this->lvariables[i]->initializeEvaluationFunction();
    	this->lvariables[i]->initializeEndowmentFunction();
    }

    // Allocate a helper array

    this->lcummulativeRates =
        new double[std::max(maxN, (int) this->lvariables.size())];
    this->ltargetChange = 0;

}


/**
 * Deallocates this simulation object.
 */
EpochSimulation::~EpochSimulation()
{
    delete[] this->lcummulativeRates;

    this->lcummulativeRates = 0;

    deallocateVector(this->lvariables);
    deallocateVector(this->lsimulationActorSets);
}


// ----------------------------------------------------------------------------
// Section: Model simulations
// ----------------------------------------------------------------------------

/**
 * Initializes the dependent variables as of the beginning of the specified
 * period.
 */
void EpochSimulation::initialize(int period)
{
	this->lperiod = period;

	// Initialize the active actor indicators

    for (unsigned i = 0; i < this->lsimulationActorSets.size(); i++)
    {
        SimulationActorSet * pActorSet = this->lsimulationActorSets[i];

        for (int i = 0; i < pActorSet->n(); i++)
        {
        	pActorSet->active(i,
        		this->lpData->active(pActorSet->pOriginalActorSet(),
        			i,
        			period));
        }
    }

    // Initialize each dependent variable

    for (unsigned i = 0; i < this->lvariables.size(); i++)
    {
      this->lvariables[i]->initialize(period);
    }

    // Initialize the effects for the upcoming simulation

    for (unsigned i = 0; i < this->lvariables.size(); i++)
    {
    	const Function * pFunction =
    		this->lvariables[i]->pEvaluationFunction();

    	for (unsigned j = 0; j < pFunction->rEffects().size(); j++)
    	{
    		pFunction->rEffects()[j]->initializeBeforeSimulation(period);
    	}

    	pFunction = this->lvariables[i]->pEndowmentFunction();

    	for (unsigned j = 0; j < pFunction->rEffects().size(); j++)
    	{
    		pFunction->rEffects()[j]->initializeBeforeSimulation(period);
    	}
    }

    // Reset the time
    this->ltime = 0;

    // Exogenous events
    this->lpEvents = this->lpData->pEventSet(period);
    this->lnextEvent = this->lpEvents->begin();

    // targets for conditional simulation
    if (this->lpModel->conditional())
    {
        this->ltargetChange = this->lpModel->rTargetChange(period);
    }
    else
    {
        this->ltargetChange = 0;
    }

    // Reset scores
    this->lscores.clear();
}


/**
 * Simulates one complete period for the model and data
 */
void EpochSimulation::runEpoch(int period)
{
    this->initialize(period);

    for(int nIter = 0; ; nIter++)
    {
        this->runStep();

        if (this->lpModel->conditional())
        {
            if (this->lpConditioningVariable->simulatedDistance() >=
            	this->ltargetChange)
            {
                break;
            }
        }
        else
        {
            if (this->ltime >= 1)
            {
                break;
            }
            else if (nIter > 1000000)
            {
#ifdef STANDALONE
				exit(1);
#endif
#ifndef STANDALONE
				error("Unlikely to terminate this epoch:",
					" more than 1000000 steps");
#endif
            }
        }
   }
	if (this->lpEvents->size())
	{
		this->setLeaversBack();
	}
}


/**
 * Simulates a single step of the actor-oriented model.
 */
void EpochSimulation::runStep()
{
    this->calculateRates();
    this->drawTimeIncrement();
    double nextTime = this->ltime + this->ltau;

    DependentVariable * pSelectedVariable = 0;
    int selectedActor = 0;

	if (this->lpModel->conditional() || nextTime < 1)
	{
		if (this->reachedCompositionChange())
		{
			this->makeNextCompositionChange();
			if (this->pModel()->needScores())
			{
				// commented out for parallel testing: bug in Siena3
				//	this->accumulateRateScores(this->ltau);
			}
		}
		else
		{
			this->ltime = nextTime;

			pSelectedVariable = this->chooseVariable();
			selectedActor = this->chooseActor(pSelectedVariable);
			// Update the scores for rate parameters
			if (this->pModel()->needScores())
			{
				this->accumulateRateScores(this->ltau,
					pSelectedVariable,
					selectedActor);
			}

			pSelectedVariable->makeChange(selectedActor);
		}
	}
	else
	{
		// Make sure we stop at 1.0 precisely.

		this->ltau = 1 - this->ltime;
		this->ltime = 1;

		// Update rate scores
		if (this->pModel()->needScores())
		{
			this->accumulateRateScores(this->ltau);
		}
	}
}


/**
 * Calculates the rates of chagne of each actor for each dependent variable and
 * the total rates of change for each variable summed over all actors.
 */
void EpochSimulation::calculateRates()
{
    for (unsigned i = 0; i < this->lvariables.size(); i++)
    {
        this->lvariables[i]->calculateRates();
    }
}


/**
 * Generates an exponential variate tau with the sum ot total rates over all
 * dependent variables as the distribution parameter. It is used later to
 * increment the current time of the simulation.
 */
void EpochSimulation::drawTimeIncrement()
{
    double totalRate = 0;

    for (unsigned i = 0; i < this->lvariables.size(); i++)
    {
        totalRate += this->lvariables[i]->totalRate();
    }

    // revert to non QAD when we have finished parallel runs
    double tau = nextExponentialQAD(totalRate);

	this->ltau = tau;
}


/**
 * Returns if the simulation has reached the time point of the next
 * exogenous event of composition change.
 */
bool EpochSimulation::reachedCompositionChange() const
{
    return this->lnextEvent != this->lpEvents->end() &&
        (*this->lnextEvent)->time() <= this->ltime + this->ltau;
}


/**
 * Makes the current composition change and resets the time of this simulation
 * to the time of the composition change.
 */
void EpochSimulation::makeNextCompositionChange()
{
	ExogenousEvent * pEvent = *this->lnextEvent;
	this->lnextEvent++;

	SimulationActorSet * pActorSet =
		this->lactorSetMap[pEvent->pActorSet()];

	if (pEvent->type() == JOINING)
	{
		pActorSet->active(pEvent->actor(), true);

		for (unsigned i = 0; i < this->lvariables.size(); i++)
		{
			this->lvariables[i]->actOnJoiner(pActorSet, pEvent->actor());
		}
	}
	else if (pEvent->type() == LEAVING)
	{
		pActorSet->active(pEvent->actor(), false);

		for (unsigned i = 0; i < this->lvariables.size(); i++)
		{
			this->lvariables[i]->actOnLeaver(pActorSet, pEvent->actor());
		}
	}

	this->ltau = pEvent->time() - this->ltime;
	this->ltime = pEvent->time();
}


/**
 * Resets the values for any actors who left the system during the current
 * period to their value at the start of the period. It will then not affect
 * the calculation of statistics. In fact resets values for all non active
 * actors.
 */
void EpochSimulation::setLeaversBack()
{
	for (unsigned i = 0; i < this->lvariables.size(); i++)
	{
// 		for (EventSet::iterator  iter = this->lpEvents->begin();
// 			 iter!=this->lpEvents->end();
// 			 iter++)
// 		{
// 			ExogenousEvent * pEvent = *iter;

// 			if (pEvent->type() == LEAVING)
// 			{
// 				this->lvariables[i]->setLeaverBack(pEvent->pActorSet(),
// 					pEvent->actor());
// 			}
// 		}

		DependentVariable *pVariable = this->lvariables[i];
		const SimulationActorSet *pActorSet = pVariable->pActorSet();

		for (int j = 0; j < pVariable->n(); j++)
		{
			if (!pActorSet->active(j))
			{
				pVariable->setLeaverBack(pActorSet, j);
			}
		}
	}
}

/**
 * Chooses one of the dependent varaibles randomly with probabilities
 * proportional to the total rate of each variable.
 */
DependentVariable * EpochSimulation::chooseVariable() const
{
	int index = 0;

	if (this->lvariables.size() > 1)
	{
		for (unsigned i = 0; i < this->lvariables.size(); i++)
		{
			this->lcummulativeRates[i] = this->lvariables[i]->totalRate();

			if (i > 0)
			{
				this->lcummulativeRates[i] += this->lcummulativeRates[i - 1];
			}
		}

		index =
			nextIntWithCumulativeProbabilities(this->lvariables.size(),
				this->lcummulativeRates);
		//	Rprintf(" %d %f %f %f\n", index, this->lcummulativeRates[0],
		//this->lcummulativeRates[1],
		//  this->lcummulativeRates[2]);
	}

	return this->lvariables[index];
}


/**
 * Chooses a random actor with probabilities proportional to the rate of change
 * for the given variable.
 */
int EpochSimulation::chooseActor(const DependentVariable * pVariable) const
{
    for (int i = 0; i < pVariable->n(); i++)
    {
        this->lcummulativeRates[i] = pVariable->rate(i);

		if (i > 0)
        {
            this->lcummulativeRates[i] += this->lcummulativeRates[i - 1];
        }
    }

    return nextIntWithCumulativeProbabilities(pVariable->n(),
        this->lcummulativeRates);
}


/**
 * Accumulates the scores for the rate parameters.
 */
void EpochSimulation::accumulateRateScores(double tau,
	const DependentVariable * pSelectedVariable,
	int selectedActor)
{
	for (unsigned i = 0; i < this->lvariables.size(); i++)
	{
		this->lvariables[i]->accumulateRateScores(tau,
			pSelectedVariable,
			selectedActor);
	}
}


// ----------------------------------------------------------------------------
// Section: Accessors
// ----------------------------------------------------------------------------

/**
 * Returns the data object underlying this simulation.
 */
const Data * EpochSimulation::pData() const
{
	return this->lpData;
}


/**
 * Returns the actor-based model simulated by this simulation object.
 */
const Model * EpochSimulation::pModel() const
{
    return this->lpModel;
}


/**
 * Returns the dependent variable with the given name if it exists;
 * otherwise 0 is returned.
 */
const DependentVariable * EpochSimulation::pVariable(string name) const
{
	DependentVariable * pVariable = 0;

	for (unsigned i = 0; i < this->lvariables.size() && !pVariable; i++)
	{
		if (this->lvariables[i]->name() == name)
		{
			pVariable = this->lvariables[i];
		}
	}

	return pVariable;
}


/**
 * Returns a reference to the vector of dependent variables.
 */
const vector<DependentVariable *> & EpochSimulation::rVariables() const
{
	return this->lvariables;
}


/**
 * Returns the wrapper actor set corresponding to the given original actor set.
 */
const SimulationActorSet * EpochSimulation::pSimulationActorSet(
	const ActorSet * pOriginalActorSet) const
{
	map<const ActorSet *, SimulationActorSet *>::const_iterator iter =
		this->lactorSetMap.find(pOriginalActorSet);

	const SimulationActorSet * pSimulationActorSet = 0;

	if (iter != this->lactorSetMap.end())
	{
		pSimulationActorSet = iter->second;
	}

	return pSimulationActorSet;
}


/**
 * Returns the currently simulated period.
 */
int EpochSimulation::period() const
{
	return this->lperiod;
}


/**
 * Returns the time taken in the simulation.
 */
double EpochSimulation::time() const
{
    return this->ltime;
}


/**
 * Returns the current score for the given effect. The scores are updated
 * in each ministep of the simulation.
 */
double EpochSimulation::score(const EffectInfo * pEffect) const
{
	map<const EffectInfo *, double>::const_iterator iter =
		this->lscores.find(pEffect);
	double score = 0;

	if (iter != this->lscores.end())
	{
		score = iter->second;
	}

	return score;
}


/**
 * Sets the score for the given effect to the given value.
 */
void EpochSimulation::score(const EffectInfo * pEffect, double value)
{
	this->lscores[pEffect] = value;
}

}
