/******************************************************************************
 * SIENA: Simulation Investigation for Empirical Network Analysis
 *
 * Web: http://www.stats.ox.ac.uk/~snijders/siena/
 *
 * File: BehaviorVariable.cpp
 *
 * Description: This file contains the implementation of the
 * BehaviorVariable class.
 *****************************************************************************/

#include <cmath>
#include <string>
#include "data/ActorSet.h"
#include "utils/Random.h"
#include "BehaviorVariable.h"
#include "data/BehaviorLongitudinalData.h"
#include "model/EpochSimulation.h"
#include "model/variables/DependentVariable.h"
#include "model/Model.h"
#include "model/effects/BehaviorEffect.h"
#include "model/EffectInfo.h"
#include "model/SimulationActorSet.h"

namespace siena
{

/**
 * Creates a new behavior variable for the given observed data.
 * @param pModel the owner model of this variable
 */
BehaviorVariable::BehaviorVariable(BehaviorLongitudinalData * pData,
	EpochSimulation * pSimulation) :
		DependentVariable(pData->name(),
			pData->pActorSet(),
			pData->observationCount(),
			pSimulation)
{
	this->lpData = pData;
	this->lvalues = new int[this->n()];
	this->lpredictorValues = new int[this->n()];
	this->levaluationEffectContribution = new double * [3];
	this->lendowmentEffectContribution = new double * [3];
	this->lprobabilities = new double[3];

	for (int i = 0; i < 3; i++)
	{
		this->levaluationEffectContribution[i] =
			new double[pSimulation->pModel()->rEvaluationEffects(pData->name()).size()];
		this->lendowmentEffectContribution[i] =
			new double[pSimulation->pModel()->rEvaluationEffects(pData->name()).size()];
	}
}


/**
 * Deallocates this variable object.
 */
BehaviorVariable::~BehaviorVariable()
{
	delete[] this->lvalues;
	delete[] this->lpredictorValues;

	this->lpData = 0;
	this->lvalues = 0;
	this->lpredictorValues = 0;
	delete[] this->lprobabilities;
	// Delete arrays of contributions

	for (int i = 0; i < 3; i++)
	{
		delete[] this->levaluationEffectContribution[i];
		delete[] this->lendowmentEffectContribution[i];
	}

	delete[] this->levaluationEffectContribution;
	delete[] this->lendowmentEffectContribution;

	this->levaluationEffectContribution = 0;
	this->lendowmentEffectContribution = 0;
	this->lprobabilities = 0;
}


// ----------------------------------------------------------------------------
// Section: Accessors
// ----------------------------------------------------------------------------

/**
 * Returns the second dimension of this variable, namely, how many values
 * correspond to each actor. This number is 1 for behavior variables.
 */
int BehaviorVariable::m() const
{
	return 1;
}


/**
 * Returns the longitudinal data object this variable is based on.
 */
LongitudinalData * BehaviorVariable::pData() const
{
	return this->lpData;
}


/**
 * Returns if the observed start value of the given actor is missing
 * for the current period.
 */
bool BehaviorVariable::missingStartValue(int actor) const
{
	return this->lpData->missing(this->period(), actor);
}


/**
 * Returns the current value on this behavior for the given actor.
 */
int BehaviorVariable::value(int actor) const
{
	return this->lvalues[actor];
}


/**
 * Returns the current value on this behavior for the given actor, which is
 * centered around the overall mean of the observed values.
 */
double BehaviorVariable::centeredValue(int actor) const
{
	return this->lvalues[actor] - this->lpData->overallMean();
}


/**
 * Returns the centered similarity of the given actors.
 */
double BehaviorVariable::similarity(int i, int j) const
{
	return this->lpData->similarity(this->lvalues[i], this->lvalues[j]);
}


/**
 * Returns the predictor value on this behavior for the given actor.
 */
int BehaviorVariable::predictorValue(int actor) const
{
	return this->lpredictorValues[actor];
}


/**
 * Returns the predictor value on this behavior for the given actor, which is
 * centered around the overall mean of the observed values.
 */
double BehaviorVariable::centeredPredictorValue(int actor) const
{
	return this->lpredictorValues[actor] - this->lpData->overallMean();
}


/**
 * Returns the array of current values for this behavior variable.
 */
const int * BehaviorVariable::values() const
{
	return this->lvalues;
}


/**
 * Returns the range of observed values for this behavior variable.
 */
int BehaviorVariable::range() const
{
	return this->lpData->range();
}

/**
 * Returns the similarity mean for this behavior variable.
 */
double BehaviorVariable::similarityMean() const
{
	return this->lpData->similarityMean();
}

// ----------------------------------------------------------------------------
// Section: Initialization at the beginning of a period
// ----------------------------------------------------------------------------

/**
 * Initializes this variable as of the beginning of the given period.
 */
void BehaviorVariable::initialize(int period)
{
	DependentVariable::initialize(period);

	// Copy the values from the corresponding observation.
	// Calculate the mean over active actors.

	this->lmean = 0;

	for (int i = 0; i < this->n(); i++)
	{
		this->lvalues[i] = this->lpData->value(period, i);

		if (this->pActorSet()->active(i))
		{
			this->lmean += this->lvalues[i];
		}
	}

	this->lmean /= this->pActorSet()->activeActorCount();
}

/**
 * Set the value of the variable for use as a predictor rather than a
 * dependent variable .
 */

void BehaviorVariable::predictorValue(int actor, int value)
{
	this->lpredictorValues[actor] = value;
}

// ----------------------------------------------------------------------------
// Section: Composition change
// ----------------------------------------------------------------------------

/**
 * Updates the state of this variable when an actor becomes active.
 */
void BehaviorVariable::actOnJoiner(const SimulationActorSet * pActorSet,
	int actor)
{
	if (pActorSet == this->pActorSet())
	{
		// Update the mean value of this variable over all active actors.

		// First calculate the sum over previously active actors (all active
		// actors except for the joiner itself).

		this->lmean *= pActorSet->activeActorCount() - 1;

		// Add the value of the joiner to the sum.
		this->lmean += this->lvalues[actor];

		// Divide by the number of active actors, so we get the average.
		this->lmean /= pActorSet->activeActorCount();
	}
}


/**
 * Updates the state of this variable when an actor becomes inactive.
 */
void BehaviorVariable::actOnLeaver(const SimulationActorSet * pActorSet,
	int actor)
{
	if (pActorSet == this->pActorSet())
	{
		// Update the mean value of this variable over all active actors.

		// First calculate the sum over previously active actors (all active
		// actors plus the leaver, who is inactive now).

		this->lmean *= pActorSet->activeActorCount() + 1;

		// Subtract the value of the leaver from the sum.
		this->lmean -= this->lvalues[actor];

		// Divide by the number of active actors, so we get the average.
		this->lmean /= pActorSet->activeActorCount();
	}
}


/**
 * Updates the current network and other variables when an actor becomes
 * inactive.
 */
void BehaviorVariable::setLeaverBack(const SimulationActorSet * pActorSet,
	int actor)
{
	if (pActorSet == this->pActorSet())
	{
		// Reset ties from the given actor to values at start

		for (int i = 0; i < this->n(); i++)
		{
			this->lvalues[actor] =	this->lpData->value(this->period(), actor);
		}
	}
}


// ----------------------------------------------------------------------------
// Section: Changing the behavior variable
// ----------------------------------------------------------------------------

/**
 * Simulates a change of the behavior according to the choice of the given
 * actor.
 */
void BehaviorVariable::makeChange(int actor)
{
	// The array of probabilities for downward change, no change, and upward
	// change, respectivelly

	double probabilities[3];

	// Calculate the probability for downward change

	if (this->lvalues[actor] > this->lpData->min() &&
		!this->lpData->upOnly(this->period()))
	{
		probabilities[0] =
			exp(this->totalEvaluationContribution(actor,
					-1) +
					 +
				this->totalEndowmentContribution(actor,
					-1));
	}
	else
	{
		probabilities[0] = 0;
	}
	// No change means zero contribution, but exp(0) = 1
	probabilities[1] = 1;

	// Calculate the probability for upward change

	if (this->lvalues[actor] < this->lpData->max() &&
		!this->lpData->downOnly(this->period()))
	{
		probabilities[2] =
			exp(this->totalEvaluationContribution(actor,
				1));
	}
	else
	{
		probabilities[2] = 0;
	}

	if (this->pSimulation()->pModel()->needScores())
	{
		double sum = 0;
		sum = probabilities[0] + 1 + probabilities[2];
		this->lprobabilities[0] =  probabilities[0] / sum;
		this->lprobabilities[1] =  probabilities[1] / sum;
		this->lprobabilities[2] =  probabilities[2] / sum;
	}

	// Transform to cumulative probabilities

	probabilities[1] += probabilities[0];
	probabilities[2] += probabilities[1];

	// Choose the change
	int difference = nextIntWithCumulativeProbabilities(3, probabilities) - 1;

	if (this->pSimulation()->pModel()->needScores())
	{
		this->accumulateScores(difference + 1);
	}

	// Make the change

	if (difference != 0)
	{
		int oldValue = this->lvalues[actor];

		// Make the change
		this->lvalues[actor] += difference;

		// Update the mean over active actors

		this->lmean +=
			((double) difference) /
				this->pActorSet()->activeActorCount();

		// Update the distance from the observed data at the beginning of the
		// period. Actors with missing values at any of the endpoints of the
		// period don't contribute to the distance

		if (!this->lpData->missing(this->period(), actor) &&
			!this->lpData->missing(this->period() + 1, actor))
		{
			int observedValue = this->lpData->value(this->period(), actor);
			this->simulatedDistance(this->simulatedDistance() +
				abs(this->lvalues[actor] - observedValue) -
				abs(oldValue - observedValue));
		}
	}
}


/**
 * Returns the total contribution of all effects in the given function if
 * the behavior of the given actor is changed by the given amount.
 */
double BehaviorVariable::totalEvaluationContribution(int actor,
	int difference) const
{
	double contribution = 0;
	const Function * pFunction = this->pEvaluationFunction();

	for (unsigned i = 0; i < pFunction->rEffects().size(); i++)
	{
		BehaviorEffect * pEffect =
			(BehaviorEffect *) pFunction->rEffects()[i];
		double thisContribution =
			pEffect->calculateChangeContribution(actor, difference);
		this->levaluationEffectContribution[difference+1][i] =
			thisContribution;
		contribution += pEffect->weight() * thisContribution;
	}
	return contribution;
}

double BehaviorVariable::totalEndowmentContribution(int actor,
	int difference) const
{
	double contribution = 0;
	const Function * pFunction = this->pEndowmentFunction();

	for (unsigned i = 0; i < pFunction->rEffects().size(); i++)
	{
		BehaviorEffect * pEffect =
			(BehaviorEffect *) pFunction->rEffects()[i];
		double thisContribution =
			pEffect->calculateChangeContribution(actor, difference);
		this->lendowmentEffectContribution[difference+1][i] =
			thisContribution;
		contribution += pEffect->weight() * thisContribution;
	}

	return contribution;
}
/**
 * Updates the scores for evaluation and endowment function effects according
 * to the current step in the simulation.
 */
void BehaviorVariable::accumulateScores(int difference) const
{
	for (unsigned i = 0;
		i < this->pEvaluationFunction()->rEffects().size();
		i++)
	{
		if (difference == 1) // no change, but not initialised
		{
			this->levaluationEffectContribution[difference][i] = 0;
		}
		Effect * pEffect = this->pEvaluationFunction()->rEffects()[i];
		double score = this->levaluationEffectContribution[difference][i];

		for (int j = 0; j < 3; j+=2)
		{
			score -=
				this->levaluationEffectContribution[j][i] *
					this->lprobabilities[j];
		}

		this->pSimulation()->score(pEffect->pEffectInfo(),
			this->pSimulation()->score(pEffect->pEffectInfo()) + score);
	}

	for (unsigned i = 0;
		i < this->pEndowmentFunction()->rEffects().size();
		i++)
	{
		if (difference == 1) // no change, but not initialised
		{
			this->lendowmentEffectContribution[difference][i] = 0;
		}
		if (difference == 2) // up has no effect on endowment function
		{
			this->lendowmentEffectContribution[difference][i] = 0;
		}
		Effect * pEffect = this->pEndowmentFunction()->rEffects()[i];
		double score = this->lendowmentEffectContribution[difference][i];

		for (int j = 0; j < 2; j+=2)
		{
			score -=
				this->lendowmentEffectContribution[j][i] *
					this->lprobabilities[j];

		}

		this->pSimulation()->score(pEffect->pEffectInfo(),
			this->pSimulation()->score(pEffect->pEffectInfo()) + score);
	}
}


}
