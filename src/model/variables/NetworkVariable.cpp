/******************************************************************************
 * SIENA: Simulation Investigation for Empirical Network Analysis
 *
 * Web: http://www.stats.ox.ac.uk/~snijders/siena/
 *
 * File: NetworkVariable.cpp
 *
 * Description: This file contains the implementation of the
 * NetworkVariable class.
 *****************************************************************************/
#include <algorithm>
#include <cmath>
#include <R.h>
#include "NetworkVariable.h"
#include "utils/Utils.h"
#include "utils/Random.h"
#include "data/ActorSet.h"
#include "data/Network.h"
#include "data/OneModeNetwork.h"
#include "data/IncidentTieIterator.h"
#include "data/NetworkLongitudinalData.h"
#include "data/OneModeNetworkLongitudinalData.h"
#include "model/EpochSimulation.h"
#include "model/Model.h"
#include "model/EffectInfo.h"
#include "model/effects/NetworkEffect.h"
#include "model/tables/TwoPathTable.h"
#include "model/tables/CriticalInStarTable.h"

namespace siena
{

// ----------------------------------------------------------------------------
// Section: Construction
// ----------------------------------------------------------------------------

/**
 * Creates a new network variable for the given observed data.
 * @param pSimulation the owner of this variable
 */
NetworkVariable::NetworkVariable(NetworkLongitudinalData * pData,
	EpochSimulation * pSimulation) :
		DependentVariable(pData->name(),
			pData->pActorSet(),
			pData->observationCount(),
			pSimulation)
{
	this->lpData = pData;
	this->lpNetwork = 0;
	this->lpPredictorNetwork = 0;
	this->lactiveStructuralTieCount = new int[this->n()];

	this->lpHasOutTie = new bool[this->m()];
	this->lpermitted = new bool[this->m()];
	this->lprobabilities = new double[this->m()];
	this->lpHasInTie = 0;

	if (this->oneModeNetwork())
	{
		this->lpHasInTie = new bool[this->n()];
	}

	this->levaluationEffectContribution = new double * [this->m()];
	this->lendowmentEffectContribution = new double * [this->m()];

	for (int i = 0; i < this->m(); i++)
	{
		this->levaluationEffectContribution[i] =
			new double[pSimulation->pModel()->rEvaluationEffects(pData->name()).size()];
		this->lendowmentEffectContribution[i] =
			new double[pSimulation->pModel()->rEvaluationEffects(pData->name()).size()];
	}

	this->initConfigurationTables();
}


/**
 * Deallocates this variable object.
 */
NetworkVariable::~NetworkVariable()
{
	this->deleteConfigurationTables();
	delete this->lpNetwork;
	delete this->lpPredictorNetwork;

	delete[] this->lactiveStructuralTieCount;
	delete[] this->lpHasOutTie;
	delete[] this->lpHasInTie;
	delete[] this->lpermitted;
	delete[] this->lprobabilities;


	// Delete arrays of contributions

	for (int i = 0; i < this->m(); i++)
	{
		delete[] this->levaluationEffectContribution[i];
		delete[] this->lendowmentEffectContribution[i];
	}

	delete[] this->levaluationEffectContribution;
	delete[] this->lendowmentEffectContribution;

	this->levaluationEffectContribution = 0;
	this->lendowmentEffectContribution = 0;
	this->lpData = 0;
	this->lpNetwork = 0;
	this->lpPredictorNetwork = 0;
	this->lactiveStructuralTieCount = 0;
	this->lpHasOutTie = 0;
	this->lpHasInTie = 0;
	this->lpermitted = 0;
	this->lprobabilities = 0;
}


/**
 * Creates various configuration tables, one of each of the available types.
 * These configurations are stored explicitly and also in a list.
 */
void NetworkVariable::initConfigurationTables()
{
	// Two-paths, reverse two-paths, in-stars, and out-stars can be seen
	// as special cases of generalized two-paths, where one can specify
	// the directions of the first and the second step. For instance, the
	// number of in-stars between i and j equals the number of possible ways
	// of reaching j, if we have to traverse one outgoing tie of actor i,
	// say (i,h), followed by on incoming tie of actor h.

	if (this->oneModeNetwork())
	{
		this->lpTwoPathTable = new TwoPathTable(this, FORWARD, FORWARD);
		this->lpReverseTwoPathTable =
			new TwoPathTable(this, BACKWARD, BACKWARD);
		this->lpOutStarTable = new TwoPathTable(this, BACKWARD, FORWARD);
		this->lpCriticalInStarTable = new CriticalInStarTable(this);

		this->lconfigurationTables.push_back(this->pTwoPathTable());
		this->lconfigurationTables.push_back(this->pReverseTwoPathTable());
		this->lconfigurationTables.push_back(this->pOutStarTable());
		this->lconfigurationTables.push_back(this->pCriticalInStarTable());
	}

	this->lpInStarTable = new TwoPathTable(this, FORWARD, BACKWARD);
	this->lconfigurationTables.push_back(this->pInStarTable());
}


/**
 * Deletes all configuration tables of this function.
 */
void NetworkVariable::deleteConfigurationTables()
{
	deallocateVector(this->lconfigurationTables);

	// It is a good idea to nullify the pointers.

	this->lpTwoPathTable = 0;
	this->lpReverseTwoPathTable = 0;
	this->lpInStarTable = 0;
	this->lpOutStarTable = 0;
	this->lpCriticalInStarTable = 0;
}


// ----------------------------------------------------------------------------
// Section: Accessors
// ----------------------------------------------------------------------------

/**
 * Returns the set of actors acting as tie senders.
 */
const ActorSet * NetworkVariable::pSenders() const
{
	return this->lpData->pSenders();
}


/**
 * Returns the set of actors acting as tie receivers.
 */
const ActorSet * NetworkVariable::pReceivers() const
{
	return this->lpData->pReceivers();
}


/**
 * Returns the second dimension of this variable, namely, how many values
 * correspond to each actor. This number equals the number of receivers for
 * network variables.
 */
int NetworkVariable::m() const
{
	return this->lpData->pReceivers()->n();
}


/**
 * Indicates if this is a one-mode network, namely, if the senders and
 * receivers are the same set of actors.
 */
bool NetworkVariable::oneModeNetwork() const
{
	return this->pSenders() == this->pReceivers();
}


/**
 * Returns the longitudinal data object this variable is based on.
 */
LongitudinalData * NetworkVariable::pData() const
{
	return this->lpData;
}


// ----------------------------------------------------------------------------
// Section: Initialization
// ----------------------------------------------------------------------------

/**
 * Initializes this variable as of the beginning of the given period.
 */
void NetworkVariable::initialize(int period)
{
	DependentVariable::initialize(period);

	// Copy the respective observation

	delete this->lpNetwork;
	delete this->lpPredictorNetwork;

	if (this->oneModeNetwork())
	{
		OneModeNetwork * pOneModeNetwork =
			(OneModeNetwork *) this->lpData->pNetwork(period);
		this->lpNetwork = new OneModeNetwork(*pOneModeNetwork);
	}
	else
	{
		this->lpNetwork = new Network(*this->lpData->pNetwork(period));
	}

	// Initialize the counters with all structural ties, including those to
	// inactive actors.

	for (int i = 0; i < this->n(); i++)
	{
		this->lactiveStructuralTieCount[i] =
			this->lpData->structuralTieCount(i, period);
	}

	// Now subtract the structural ties to initially inactive actors.

	for (int i = 0; i < this->m(); i++)
	{
		if (!this->pSimulation()->active(this->pReceivers(), i))
		{
			for (IncidentTieIterator iter =
					this->lpData->pStructuralTieNetwork(period)->inTies(i);
				iter.valid();
				iter.next())
			{
				this->lactiveStructuralTieCount[iter.actor()]--;
			}
		}
	}

}


/**
 * Returns the current state of the network.
 */
Network * NetworkVariable::pNetwork() const
{
	return this->lpNetwork;
}

/**
 * sets the predictor network.
 */
void NetworkVariable::pPredictorNetwork(Network * pPredictorNetwork)
{
	this->lpPredictorNetwork = pPredictorNetwork;
}

/**
 * Returns the state of this variable at the
 * beginning of the current period, less missings.
 */
Network * NetworkVariable::pPredictorNetwork() const
{
	return this->lpPredictorNetwork;
}

// ----------------------------------------------------------------------------
// Section: Composition change
// ----------------------------------------------------------------------------

/**
 * Updates the current network and other variables when an actor becomes
 * active.
 */
void NetworkVariable::actOnJoiner(const ActorSet * pActorSet, int actor)
{
	const Network * pStartNetwork = this->lpData->pNetwork(this->period());

	if (pActorSet == this->pSenders())
	{
		// Activate the ties to active receivers according to the
		// initial observation of the period.

		for (IncidentTieIterator iter = pStartNetwork->outTies(actor);
			iter.valid();
			iter.next())
		{
			if (this->pSimulation()->active(this->pReceivers(), iter.actor()))
			{
				this->lpNetwork->setTieValue(actor,
					iter.actor(),
					iter.value());
			}
		}
	}

	if (pActorSet == this->pReceivers())
	{
		// Activate the ties from active senders according to the
		// initial observation of the period.

		for (IncidentTieIterator iter = pStartNetwork->inTies(actor);
			iter.valid();
			iter.next())
		{
			if (this->pSimulation()->active(this->pSenders(), iter.actor()))
			{
				this->lpNetwork->setTieValue(iter.actor(),
					actor,
					iter.value());
			}
		}

		// Update the numbers of structural ties to active actors, as one actor
		// becomes active now.

		for (IncidentTieIterator iter =
				this->lpData->pStructuralTieNetwork(this->period())->inTies(
					actor);
			iter.valid();
			iter.next())
		{
			this->lactiveStructuralTieCount[iter.actor()]++;
		}
	}
}


/**
 * Updates the current network and other variables when an actor becomes
 * inactive.
 */
void NetworkVariable::actOnLeaver(const ActorSet * pActorSet, int actor)
{
	if (pActorSet == this->pSenders())
	{
		// Remove the ties from the given actor.
		this->lpNetwork->clearOutTies(actor);
	}

	if (pActorSet == this->pReceivers())
	{
		// Remove the ties to the given actor.
		this->lpNetwork->clearInTies(actor);

		// Update the numbers of structural ties to active actors, as one actor
		// becomes inactive now.

		for (IncidentTieIterator iter =
				this->lpData->pStructuralTieNetwork(this->period())->inTies(
					actor);
			iter.valid();
			iter.next())
		{
			this->lactiveStructuralTieCount[iter.actor()]--;
		}
	}
}

/**
 * Sets leavers values back to the value at the start of the simulation.
 *
 */
void NetworkVariable::setLeaverBack(const ActorSet * pActorSet, int actor)
{
	if (pActorSet == this->pSenders())
	{
		// Reset ties from the given actor to values at start
		for (int i = 0; i < this->m(); i++)
		{
			if (i != actor)
			{
				this->lpNetwork->setTieValue(actor, i,
				this->lpData->pNetwork(this->period())->tieValue(actor, i));
			}
		}
	}
	if (pActorSet == this->pReceivers())
	{
		for (int i = 0; i < this->n(); i++)
		{
			if (i != actor)
			{
			this->lpNetwork->setTieValue(i, actor,
				this->lpData->pNetwork(this->period())->tieValue(i, actor));
			}
		}
	}
}

// ----------------------------------------------------------------------------
// Section: Changing the network
// ----------------------------------------------------------------------------

/**
 * Returns if the given actor can change the current state of this variable,
 * namely, it is active and at least one of the tie variables to active
 * receivers is not structurally determined.
 */
bool NetworkVariable::canMakeChange(int actor) const
{
	int activeAlterCount =
		this->pSimulation()->activeActorCount(this->pReceivers());

	if (this->oneModeNetwork())
	{
		// No loops are possible in one mode networks
		activeAlterCount--;
	}

	return this->pSimulation()->active(this->pSenders(), actor) &&
		this->lactiveStructuralTieCount[actor] < activeAlterCount;
}


/**
 * Simulates a change of the network according to the choice of the given
 * actor. First, the actor chooses the alter based on the evaluation and
 * endowment functions, and then the tie to the selected alter is flipped.
 */
void NetworkVariable::makeChange(int actor)
{
	this->lego = actor;
	this->preprocessEgo();
	this->calculatePermissibleChanges();
	this->calculateTieFlipContributions();
	this->calculateTieFlipProbabilities();
	int alter = nextIntWithProbabilities(this->m(), this->lprobabilities);

	if (this->pSimulation()->pModel()->needScores())
	{
		this->accumulateScores(alter);
	}
	// Make a change if we have a real alter (other than the ego)

	if (!this->oneModeNetwork() || this->lego != alter)
	{
		int currentValue = this->lpNetwork->tieValue(this->lego, alter);

		// Update the distance from the observed data at the beginning of the
		// period. Ties missing at any of the endpoints of the period
		// don't contribute to the distance

		if (!this->lpData->missing(this->lego, alter, this->period()) &&
			!this->lpData->missing(this->lego, alter, this->period() + 1))
		{
			if (this->lpData->tieValue(this->lego, alter, this->period()) ==
				currentValue)
			{
				this->simulatedDistance(this->simulatedDistance() + 1);
			}
			else
			{
				this->simulatedDistance(this->simulatedDistance() - 1);
			}
		}

		this->lpNetwork->setTieValue(this->lego, alter, 1 - currentValue);
	}
}


/**
 * This method does some preprocessing to speed up subsequent queries regarding
 * the current ego.
 */
void NetworkVariable::preprocessEgo()
{
	// Initialize the arrays of incoming and outgoing tie indicators.

	for (int i = 0; i < this->lpNetwork->m(); i++)
	{
		this->lpHasOutTie[i] = false;

		if (this->oneModeNetwork())
		{
			this->lpHasInTie[i] = false;
		}
	}

	// Mark the actors receiving a tie from the ego.

	for (IncidentTieIterator iter = this->lpNetwork->outTies(this->lego);
		iter.valid();
		iter.next())
	{
		this->lpHasOutTie[iter.actor()] = true;
	}

	if (this->oneModeNetwork())
	{
		// Mark the actors sending a tie to the ego.

		for (IncidentTieIterator iter = this->lpNetwork->inTies(this->lego);
			iter.valid();
			iter.next())
		{
			this->lpHasInTie[iter.actor()] = true;
		}
	}

	// Invalidate all configuration tables such that they are recalculated
	// the first time they are queried.

	std::for_each(this->lconfigurationTables.begin(),
		this->lconfigurationTables.end(),
		std::mem_fun(&ConfigurationTable::invalidate));

	// Calculate the configuration tables that are required for calculating
	// the selected effects.

	for (unsigned i = 0; i < this->lconfigurationTables.size(); i++)
	{
		ConfigurationTable * pTable = this->lconfigurationTables[i];

		if (this->required(pTable))
		{
			pTable->calculate();
		}
	}

	// Let the effects do their preprocessing.

	const Function * pFunction = this->pEvaluationFunction();

	for (unsigned i = 0; i < pFunction->rEffects().size(); i++)
	{
		NetworkEffect * pEffect =
			(NetworkEffect *) pFunction->rEffects()[i];
		pEffect->preprocessEgo();
	}

	pFunction = this->pEndowmentFunction();

	for (unsigned i = 0; i < pFunction->rEffects().size(); i++)
	{
		NetworkEffect * pEffect =
			(NetworkEffect *) pFunction->rEffects()[i];
		pEffect->preprocessEgo();
	}
}


/**
 * Tests if the given configuration table is required for the calculation
 * of at least one of the selected effects.
 */
bool NetworkVariable::required(const ConfigurationTable * pTable) const
{
	bool rc = false;
	const Function * pFunction = this->pEvaluationFunction();

	for (unsigned i = 0; i < pFunction->rEffects().size() && !rc; i++)
	{
		NetworkEffect * pEffect =
			(NetworkEffect *) pFunction->rEffects()[i];

		if (pEffect->usesTable(pTable))
		{
			rc = true;
		}
	}

	pFunction = this->pEndowmentFunction();

	for (unsigned i = 0; i < pFunction->rEffects().size() && !rc; i++)
	{
		NetworkEffect * pEffect =
			(NetworkEffect *) pFunction->rEffects()[i];

		if (pEffect->usesTable(pTable))
		{
			rc = true;
		}
	}

	return rc;
}


/**
 * Determines the set of actors that are allowed to act as alters in the next
 * tie flip, and stores this information in <code>lpermitted</code>.
 */
void NetworkVariable::calculatePermissibleChanges()
{
	NetworkLongitudinalData * pData =
		(NetworkLongitudinalData *) this->pData();

	// Test each alter if a tie flip to that alter is permitted according
	// to upOnly/downOnly flags and the max degree.

	for (int i = 0; i < this->m(); i++)
	{
		if (this->lpHasOutTie[i])
		{
			this->lpermitted[i] = !pData->upOnly(this->period());
		}
		else if (i != this->lego || !this->oneModeNetwork())
		{
			this->lpermitted[i] =
				// for comparability with siena3 comment this out
				//	this->pSimulation()->active(this->pReceivers(), i) &&
				!pData->downOnly(this->period()) &&
				this->lpNetwork->outDegree(this->lego) < pData->maxDegree();
		}
		else
		{
			// It is okay to not make any change at all
			this->lpermitted[i] = true;
		}
	}

	// Prohibit the change of structural ties

	for (IncidentTieIterator iter =
			pData->pStructuralTieNetwork(this->period())->outTies(this->lego);
		iter.valid();
		iter.next())
	{
		this->lpermitted[iter.actor()] = false;
	}
}


/**
 * For each alter, this method calculates the contribution of each evaluation
 * and endowment effect if a tie from the ego to this alter was flipped.
 * These contributions are stored in arrays
 * <code>levaluationEffectContribution</code> and
 * <code>lendowmentEffectContribution</code>.
 */
void NetworkVariable::calculateTieFlipContributions()
{
	int evaluationEffectCount = this->pEvaluationFunction()->rEffects().size();
	int endowmentEffectCount = this->pEndowmentFunction()->rEffects().size();

	for (int alter = 0; alter < this->m(); alter++)
	{
		// alter = ego for one-mode networks means no change.
		// No change, no contribution.

		if (this->lpermitted[alter] &&
			(!this->oneModeNetwork() || alter != this->lego))
		{
			for (int i = 0; i < evaluationEffectCount; i++)
			{
				NetworkEffect * pEffect =
					(NetworkEffect *)
					this->pEvaluationFunction()->rEffects()[i];
				this->levaluationEffectContribution[alter][i] =
					pEffect->calculateTieFlipContribution(alter);
			}
		}
		else
		{
			for (int i = 0; i < evaluationEffectCount; i++)
			{
				this->levaluationEffectContribution[alter][i] = 0;
			}
		}

		// The endowment effects have non-zero contributions on tie
		// withdrawals only

		if (this->lpermitted[alter] &&
			(!this->oneModeNetwork() || alter != this->lego) &&
			this->lpHasOutTie[alter])
		{
			for (int i = 0; i < endowmentEffectCount; i++)
			{
				NetworkEffect * pEffect =
					(NetworkEffect *)
						this->pEndowmentFunction()->rEffects()[i];
				this->lendowmentEffectContribution[alter][i] =
					pEffect->calculateTieFlipContribution(alter);
			}
		}
		else
		{
			for (int i = 0; i < endowmentEffectCount; i++)
			{
				this->lendowmentEffectContribution[alter][i] = 0;
			}
		}
	}
}


/**
 * Calculates the probability of each actor of being chosen for the next
 * tie flip.
 */
void NetworkVariable::calculateTieFlipProbabilities()
{
	int evaluationEffectCount = this->pEvaluationFunction()->rEffects().size();
	int endowmentEffectCount = this->pEndowmentFunction()->rEffects().size();

	double total = 0;

	for (int alter = 0; alter < this->m(); alter++)
	{
		if (this->lpermitted[alter])
		{
			// Calculate the total contribution of all effects

			double contribution = 0;

			for (int i = 0; i < evaluationEffectCount; i++)
			{
				Effect * pEffect = this->pEvaluationFunction()->rEffects()[i];
				contribution +=
					pEffect->weight() *
						this->levaluationEffectContribution[alter][i];
			}

			if (this->lpHasOutTie[alter])
			{
				for (int i = 0; i < endowmentEffectCount; i++)
				{
					Effect * pEffect =
						this->pEndowmentFunction()->rEffects()[i];
					contribution +=
						pEffect->weight() *
							this->lendowmentEffectContribution[alter][i];
				}
			}

			// The selection probability is the exponential of the total
			// contribution.

			this->lprobabilities[alter] = exp(contribution);
		}
		else
		{
			this->lprobabilities[alter] = 0;
		}

		total += this->lprobabilities[alter];
	}

	// Normalize

	if (total != 0)
	{
		for (int alter = 0; alter < this->m(); alter++)
		{
			this->lprobabilities[alter] /= total;
		}
	}
}


/**
 * Updates the scores for evaluation and endowment function effects according
 * to the current step in the simulation.
 */
void NetworkVariable::accumulateScores(int alter) const
{
	for (unsigned i = 0;
		i < this->pEvaluationFunction()->rEffects().size();
		i++)
	{
		Effect * pEffect = this->pEvaluationFunction()->rEffects()[i];
		double score = this->levaluationEffectContribution[alter][i];
		//Rprintf("score 1 %f\n", score);

		for (int j = 0; j < this->m(); j++)
		{
			score -=
				this->levaluationEffectContribution[j][i] *
					this->lprobabilities[j];
			//		Rprintf("score 2 %d %f %f %f\n", j, score,
			//			this->levaluationEffectContribution[j][i],
			//		this->lprobabilities[j]);
		}

		this->pSimulation()->score(pEffect->pEffectInfo(),
			this->pSimulation()->score(pEffect->pEffectInfo()) + score);
//	Rprintf("score %f\n", score);
	}
	for (unsigned i = 0;
		i < this->pEndowmentFunction()->rEffects().size();
		i++)
	{
		Effect * pEffect = this->pEndowmentFunction()->rEffects()[i];
		double score = this->lendowmentEffectContribution[alter][i];

		for (int j = 0; j < this->m(); j++)
		{
			score -=
				this->lendowmentEffectContribution[j][i] *
					this->lprobabilities[j];
		}

		this->pSimulation()->score(pEffect->pEffectInfo(),
			this->pSimulation()->score(pEffect->pEffectInfo()) + score);
//	Rprintf(" endow score %f\n", score);
	}
}

}
