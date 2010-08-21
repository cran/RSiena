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
#include <R_ext/Print.h>
#include <R_ext/Arith.h>
#include <Rinternals.h>
#include <algorithm>
#include <cmath>
#include "NetworkVariable.h"
#include "utils/Utils.h"
#include "utils/Random.h"
#include "data/ActorSet.h"
#include "network/Network.h"
#include "network/OneModeNetwork.h"
#include "network/IncidentTieIterator.h"
#include "data/NetworkLongitudinalData.h"
#include "data/OneModeNetworkLongitudinalData.h"
#include "model/EpochSimulation.h"
#include "model/SimulationActorSet.h"
#include "model/Model.h"
#include "model/EffectInfo.h"
#include "model/effects/NetworkEffect.h"
#include "model/tables/Cache.h"
#include "model/tables/NetworkCache.h"
#include "model/ml/NetworkChange.h"
#include "model/filters/PermittedChangeFilter.h"
#include "model/ml/Chain.h"

namespace siena
{
SEXP getMiniStepDF(const MiniStep& miniStep);

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
			pSimulation)
{
	this->lpData = pData;
	this->lpSenders = pSimulation->pSimulationActorSet(pData->pSenders());
	this->lpReceivers = pSimulation->pSimulationActorSet(pData->pReceivers());
	this->lpNetwork = 0;
	this->lactiveStructuralTieCount = new int[this->n()];

	this->lpermitted = new bool[this->m()];

	int numberOfAlters;
	if (this->oneModeNetwork())
	{
		this->lpNetwork = new OneModeNetwork(this->n(), false);
		numberOfAlters = this->m();
	}
	else
	{
		this->lpNetwork = new Network(this->n(), this->m());
		numberOfAlters = this->m() + 1;
	}

	this->lprobabilities = new double[numberOfAlters];
	this->levaluationEffectContribution = new double * [numberOfAlters];
	this->lendowmentEffectContribution = new double * [numberOfAlters];

	for (int i = 0; i < numberOfAlters; i++)
	{
		this->levaluationEffectContribution[i] =
			new double[pSimulation->pModel()->
				rEvaluationEffects(pData->name()).size()];
		this->lendowmentEffectContribution[i] =
			new double[pSimulation->pModel()->
				rEvaluationEffects(pData->name()).size()];
	}

	this->lpNetworkCache =
		pSimulation->pCache()->pNetworkCache(this->lpNetwork);
}


/**
 * Deallocates this variable object.
 */
NetworkVariable::~NetworkVariable()
{
	delete this->lpNetwork;

	delete[] this->lactiveStructuralTieCount;
	delete[] this->lpermitted;
	delete[] this->lprobabilities;


	// Delete arrays of contributions
	int numberOfAlters;
	if (this->oneModeNetwork())
	{
		numberOfAlters = this->m();
	}
	else
	{
		numberOfAlters = this->m() + 1;
	}
	for (int i = 0; i < numberOfAlters; i++)
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
	this->lactiveStructuralTieCount = 0;
	this->lpermitted = 0;
	this->lprobabilities = 0;

	deallocateVector(this->lpermittedChangeFilters);
}


/**
 * Adds the given filter of permissible changes to the filters of this
 * variable. This variable becomes the owner of the filter, which means that
 * the filter will be deleted as soon as this variable is deleted.
 */
void NetworkVariable::addPermittedChangeFilter(PermittedChangeFilter * pFilter)
{
	this->lpermittedChangeFilters.push_back(pFilter);
}



// ----------------------------------------------------------------------------
// Section: Accessors
// ----------------------------------------------------------------------------

/**
 * Returns the set of actors acting as tie senders.
 */
const SimulationActorSet * NetworkVariable::pSenders() const
{
	return this->lpSenders;
}


/**
 * Returns the set of actors acting as tie receivers.
 */
const SimulationActorSet * NetworkVariable::pReceivers() const
{
	return this->lpReceivers;
}


/**
 * Returns the second dimension of this variable, namely, how many values
 * correspond to each actor. This number equals the number of receivers for
 * network variables.
 */
int NetworkVariable::m() const
{
	return this->lpReceivers->n();
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

	if (this->oneModeNetwork())
	{
		OneModeNetwork * pNetwork =
			(OneModeNetwork *) this->lpNetwork;
		const OneModeNetwork * pObservedNetwork =
			(const OneModeNetwork *) this->lpData->pNetwork(period);

		// Use the copy assignment operator
		(*pNetwork) = (*pObservedNetwork);
	}
	else
	{
		// Use the copy assignment operator
		(*this->lpNetwork) = (*this->lpData->pNetwork(period));
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
		if (!this->pReceivers()->active(i))
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


// ----------------------------------------------------------------------------
// Section: Composition change
// ----------------------------------------------------------------------------

/**
 * Updates the current network and other variables when an actor becomes
 * active.
 */
void NetworkVariable::actOnJoiner(const SimulationActorSet * pActorSet,
	int actor)
{
	DependentVariable::actOnJoiner(pActorSet, actor);

	const Network * pStartNetwork = this->lpData->pNetwork(this->period());

	if (pActorSet == this->pSenders())
	{
		// Activate the ties to active receivers according to the
		// initial observation of the period.

		for (IncidentTieIterator iter = pStartNetwork->outTies(actor);
			iter.valid();
			iter.next())
		{
			if (this->pReceivers()->active(iter.actor()))
			{
				this->lpNetwork->setTieValue(actor,
					iter.actor(),
					iter.value());
			}
		}

		// The rates need to be recalculated.
		this->invalidateRates();
	}

	if (pActorSet == this->pReceivers())
	{
		// Activate the ties from active senders according to the
		// initial observation of the period.

		for (IncidentTieIterator iter = pStartNetwork->inTies(actor);
			iter.valid();
			iter.next())
		{
			if (this->pSenders()->active(iter.actor()))
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

		// The rates need to be recalculated.
		this->invalidateRates();
	}
}


/**
 * Updates the current network and other variables when an actor becomes
 * inactive.
 */
void NetworkVariable::actOnLeaver(const SimulationActorSet * pActorSet,
	int actor)
{
	DependentVariable::actOnLeaver(pActorSet, actor);

	if (pActorSet == this->pSenders())
	{
		// Remove the ties from the given actor.
		this->lpNetwork->clearOutTies(actor);

		// The rates need to be recalculated.
		this->invalidateRates();
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

		// The rates need to be recalculated.
		this->invalidateRates();
	}
}


/**
 * Sets leavers values back to the value at the start of the simulation.
 *
 */
void NetworkVariable::setLeaverBack(const SimulationActorSet * pActorSet,
	int actor)
{
	if (pActorSet == this->pSenders())
	{
		// Reset ties from the given actor to values at start

		for (int i = 0; i < this->m(); i++)
		{
			if (i != actor)
			{
				this->lpNetwork->setTieValue(actor,
					i,
					this->lpData->tieValue(actor, i, this->period()));
			}
		}
	}

	if (pActorSet == this->pReceivers())
	{
		for (int i = 0; i < this->n(); i++)
		{
			if (i != actor)
			{
				this->lpNetwork->setTieValue(i,
					actor,
					this->lpData->tieValue(i, actor, this->period()));
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
	bool rc = DependentVariable::canMakeChange(actor);

	if (rc)
	{
		int activeAlterCount =
			this->lpReceivers->activeActorCount();

		if (this->oneModeNetwork())
		{
			// No loops are possible in one mode networks
			activeAlterCount--;
		}

		rc &= this->lpSenders->active(actor) &&
			this->lactiveStructuralTieCount[actor] < activeAlterCount;
	}

	return rc;
}


/**
 * Simulates a change of the network according to the choice of the given
 * actor. First, the actor chooses the alter based on the evaluation and
 * endowment functions, and then the tie to the selected alter is flipped.
 */
void NetworkVariable::makeChange(int actor)
{
	int m;
	this->lego = actor;
	this->calculateTieFlipProbabilities();
	if (this->oneModeNetwork())
	{
		m = this->m();
	}
	else
	{
		m = this->m() + 1;
	}

	int alter = nextIntWithProbabilities(m, this->lprobabilities);

	if (this->pSimulation()->pModel()->needScores())
	{
		this->accumulateScores(alter);
	}
	if (this->pSimulation()->pModel()->needChain())
	{
		// add ministep to chain
		MiniStep * pMiniStep =
			new NetworkChange(this->lpData, actor, alter);
		this->pSimulation()->pChain()->insertBefore(pMiniStep,
			this->pSimulation()->pChain()->pLast());
		pMiniStep->logChoiceProbability(log(this->lprobabilities[alter]));
		if (this->pSimulation()->pModel()->needChangeContributions())
		{
			this->copyChangeContributions(pMiniStep);
		}
	}
	// Make a change if we have a real alter (other than the ego or
	// the dummy for bipartite networks)

	if ((!this->oneModeNetwork() && alter < this->m()) ||
		(this->oneModeNetwork() && this->lego != alter))
	{
		int currentValue = this->lpNetwork->tieValue(this->lego, alter);

		// Update the distance from the observed data at the beginning of the
		// period. Ties missing at any of the endpoints of the period
		// don't contribute to the distance

		// If network is symmetric, changes are in steps of two, otherwise 1
		int change = 1;
		if (this->oneModeNetwork())
		{
			if ((dynamic_cast<const OneModeNetworkLongitudinalData *>
					(this->lpData))->symmetric())
			{
				change = 2;
			}
		}
		if (!this->lpData->missing(this->lego, alter, this->period()) &&
			!this->lpData->missing(this->lego, alter, this->period() + 1))
		{
			if (this->lpData->tieValue(this->lego, alter, this->period()) ==
				currentValue)
			{
				this->simulatedDistance(this->simulatedDistance() + change);
			}
			else
			{
				this->simulatedDistance(this->simulatedDistance() - change);
			}
		}

		this->lpNetwork->setTieValue(this->lego, alter, 1 - currentValue);

		if (this->oneModeNetwork())
		{
			const OneModeNetworkLongitudinalData * pData =
				dynamic_cast<const OneModeNetworkLongitudinalData *>(
					this->pData());

			if (pData->symmetric())
			{
				this->lpNetwork->setTieValue(alter,
					this->lego,
					1 - currentValue);
			}
		}
	}
}


/**
 * This method does some preprocessing to speed up subsequent queries regarding
 * the current ego.
 */
void NetworkVariable::preprocessEgo()
{
	// Let the effects do their preprocessing.

	const Function * pFunction = this->pEvaluationFunction();

	for (unsigned i = 0; i < pFunction->rEffects().size(); i++)
	{
		NetworkEffect * pEffect =
			(NetworkEffect *) pFunction->rEffects()[i];
		pEffect->preprocessEgo(this->lego);
	}

	pFunction = this->pEndowmentFunction();

	for (unsigned i = 0; i < pFunction->rEffects().size(); i++)
	{
		NetworkEffect * pEffect =
			(NetworkEffect *) pFunction->rEffects()[i];
		pEffect->preprocessEgo(this->lego);
	}
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

	int m = this->m();

	for (int i = 0; i < m; i++)
	{
		if (this->lpNetworkCache->outTieExists(i))
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

	// Run the filters that may forbid some more changes

	for (unsigned i = 0; i < this->lpermittedChangeFilters.size(); i++)
	{
		PermittedChangeFilter * pFilter = this->lpermittedChangeFilters[i];
		pFilter->filterPermittedChanges(this->lego, this->lpermitted);
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
	const vector<Effect *> & rEvaluationEffects =
		this->pEvaluationFunction()->rEffects();
	const vector<Effect *> & rEndowmentEffects =
		this->pEndowmentFunction()->rEffects();
	bool twoModeNetwork = !this->oneModeNetwork();
	int m = this->m();

	for (int alter = 0; alter < m; alter++)
	{
		// alter = ego for one-mode networks means no change.
		// No change, no contribution.

		if (this->lpermitted[alter] &&
			(twoModeNetwork || alter != this->lego))
		{
			for (int i = 0; i < evaluationEffectCount; i++)
			{
				NetworkEffect * pEffect =
					(NetworkEffect *) rEvaluationEffects[i];
				double contribution =
					pEffect->calculateContribution(alter);

				// Tie withdrawals contribute in the opposite way

				if (this->lpNetworkCache->outTieExists(alter))
				{
					contribution = -contribution;
				}

				this->levaluationEffectContribution[alter][i] = contribution;
			}
		}
		else
		{
			for (int i = 0; i < evaluationEffectCount; i++)
			{
				if (!this->lpermitted[alter])
				{
					this->levaluationEffectContribution[alter][i] = R_NaN;
				}
				else
				{
				this->levaluationEffectContribution[alter][i] = 0;
				}
			}
		}

		// The endowment effects have non-zero contributions on tie
		// withdrawals only

		if (this->lpNetworkCache->outTieExists(alter) &&
			this->lpermitted[alter])
		{
			for (int i = 0; i < endowmentEffectCount; i++)
			{
				NetworkEffect * pEffect =
					(NetworkEffect *) rEndowmentEffects[i];
				this->lendowmentEffectContribution[alter][i] =
					-pEffect->calculateContribution(alter);
			}
		}
		else
		{
			for (int i = 0; i < endowmentEffectCount; i++)
			{
				this->lendowmentEffectContribution[alter][i] = R_NaN;
			}
		}
	}

	// need to initialise the no-effect option of alter = this->m() for
	// twomode networks

	if (twoModeNetwork)
	{
		for (int i = 0; i < evaluationEffectCount; i++)
		{
			this->levaluationEffectContribution[m][i] = 0;
		}
		for (int i = 0; i < evaluationEffectCount; i++)
		{
			this->lendowmentEffectContribution[m][i] = 0;
		}
	}
}


/**
 * Calculates the probability of each actor of being chosen for the next
 * tie flip.
 */
void NetworkVariable::calculateTieFlipProbabilities()
{
	this->preprocessEgo();
	this->calculatePermissibleChanges();
	this->calculateTieFlipContributions();

	int evaluationEffectCount = this->pEvaluationFunction()->rEffects().size();
	int endowmentEffectCount = this->pEndowmentFunction()->rEffects().size();

	double total = 0;
	int m = this->m();

	for (int alter = 0; alter < m; alter++)
	{
		if (this->lpermitted[alter])
		{
			// Calculate the total contribution of all effects

			double contribution = 0;

			for (int i = 0; i < evaluationEffectCount; i++)
			{
				Effect * pEffect = this->pEvaluationFunction()->rEffects()[i];
				contribution +=
					pEffect->parameter() *
						this->levaluationEffectContribution[alter][i];
			}

			if (this->lpNetworkCache->outTieExists(alter))
			{
				for (int i = 0; i < endowmentEffectCount; i++)
				{
					Effect * pEffect =
						this->pEndowmentFunction()->rEffects()[i];
					contribution +=
						pEffect->parameter() *
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
	if (!this->oneModeNetwork())
	{
		total += 1;
		this->lprobabilities[m] = 1;
		m++;
	}

	// Normalize

	if (total != 0)
	{
		for (int alter = 0; alter < m; alter++)
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

	int m = this->m();

	for (unsigned i = 0;
		i < this->pEvaluationFunction()->rEffects().size();
		i++)
	{
		Effect * pEffect = this->pEvaluationFunction()->rEffects()[i];
		double score = this->levaluationEffectContribution[alter][i];
		//	Rprintf("score 1 %f\n", score);

		for (int j = 0; j < m; j++)
		{
			if (this->lpermitted[j])
			{
				score -=
					this->levaluationEffectContribution[j][i] *
					this->lprobabilities[j];
			//Rprintf("score 2 %d %f %f %f\n", j, score,
			//		this->levaluationEffectContribution[j][i],
			//	this->lprobabilities[j]);
			}
		}

		this->pSimulation()->score(pEffect->pEffectInfo(),
			this->pSimulation()->score(pEffect->pEffectInfo()) + score);
		//Rprintf("score %f\n", score);
	}
	for (unsigned i = 0;
		 i < this->pEndowmentFunction()->rEffects().size();
		 i++)
	{
		Effect * pEffect = this->pEndowmentFunction()->rEffects()[i];

		double score = 0;
		if (this->lpNetworkCache->outTieExists(alter))
		{
			score += this->lendowmentEffectContribution[alter][i];
		}

		for (int j = 0; j < m; j++)
		{
			if (this->lpNetworkCache->outTieExists(j) &&
				this->lpermitted[j])
			{
				score -=
					this->lendowmentEffectContribution[j][i] *
					this->lprobabilities[j];
			}
		}

		this->pSimulation()->score(pEffect->pEffectInfo(),
			this->pSimulation()->score(pEffect->pEffectInfo()) + score);
	}
}


// ----------------------------------------------------------------------------
// Section: Maximum likelihood related methods
// ----------------------------------------------------------------------------

/**
 * Calculates the probability of the given ministep assuming that the
 * ego of the ministep will change this variable.
 */
double NetworkVariable::probability(MiniStep * pMiniStep)
{
	// Initialize the cache object for the current ego
	this->pSimulation()->pCache()->initialize(pMiniStep->ego());

	NetworkChange * pNetworkChange =
		dynamic_cast<NetworkChange *>(pMiniStep);
	this->lego = pNetworkChange->ego();
	this->calculateTieFlipProbabilities();
	if (this->pSimulation()->pModel()->needScores())
	{
		this->accumulateScores(pNetworkChange->alter());
	}
	if (this->pSimulation()->pModel()->needDerivatives())
	{
		this->accumulateDerivatives();
	}
	if (this->pSimulation()->pModel()->needChangeContributions())
	{
		this->copyChangeContributions(pMiniStep);
	}
	return this->lprobabilities[pNetworkChange->alter()];
}

/**
 * Updates the derivatives for evaluation and endowment function effects
 * according to the current miniStep in the chain.
 */
void NetworkVariable::accumulateDerivatives() const
{
	int totalEvaluationEffects = this->pEvaluationFunction()->rEffects().size();
	int totalEndowmentEffects = this->pEndowmentFunction()->rEffects().size();
	int totalEffects = totalEvaluationEffects + totalEndowmentEffects;
	Effect * pEffect1;
	Effect * pEffect2;
	double derivative;
	double * product = new double[totalEffects];
	double contribution1 = 0.0;
	double contribution2 = 0.0;

	for (int effect1 = 0; effect1 < totalEffects; effect1++)
	{
		product[effect1] = 0.0;

		if (effect1 < totalEvaluationEffects)
		{
			pEffect1 = this->pEvaluationFunction()->rEffects()[effect1];
		}
		else
		{
			pEffect1 = this->pEndowmentFunction()->rEffects()[effect1];
		}
		for (int alter = 0; alter < this->m(); alter++)
		{
			if (effect1 < totalEvaluationEffects)
			{
				product[effect1] +=
					this->levaluationEffectContribution[alter][effect1] *
					this->lprobabilities[alter];
			}
			else
			{
				product[effect1] +=
					this->lendowmentEffectContribution[alter][effect1] *
					this->lprobabilities[alter];
			}
			//	Rprintf("%d %d %f\n", alter, effect1, product[effect1]);
		}
		for (int effect2 = effect1; effect2 < totalEffects; effect2++)
		{
			derivative = 0.0;
			if (effect2 <= totalEvaluationEffects)
			{
				pEffect2 = this->pEvaluationFunction()->rEffects()[effect2];
			}
			else
			{
				pEffect2 = this->pEndowmentFunction()->rEffects()[effect2];
			}

			for (int alter = 0; alter < this->m(); alter++)
			{
				if (effect1 < totalEvaluationEffects)
				{
					contribution1 =
						this->levaluationEffectContribution[alter][effect1];
				}
				else
				{
					contribution1 =
						this->lendowmentEffectContribution[alter][effect1];
				}
				if (effect2 < totalEvaluationEffects)
				{
					contribution2 =
						this->levaluationEffectContribution[alter][effect2];
				}
				else
				{
					contribution2 =
						this->lendowmentEffectContribution[alter][effect2];
				}

				derivative -=
					contribution1 * contribution2 *	this->lprobabilities[alter];
				//	Rprintf("deriv 2 %d %d %d %f %f %f %f\n", alter, effect1,
				//    effect2,
				//		derivative,
				//		this->levaluationEffectContribution[alter][effect1],
				//		this->levaluationEffectContribution[alter][effect2],
				//		this->lprobabilities[alter]);
			}

			this->pSimulation()->derivative(pEffect1->pEffectInfo(),
				pEffect2->pEffectInfo(),
				this->pSimulation()->derivative(pEffect1->pEffectInfo(),
					pEffect2->pEffectInfo()) +	derivative);
		}
	}

	for (int effect1 = 0; effect1 < totalEffects; effect1++)
	{
		for (int effect2 = effect1; effect2 < totalEffects; effect2++)
		{
			if (effect1 < totalEvaluationEffects)
			{
				pEffect1 = this->pEvaluationFunction()->rEffects()[effect1];
			}
			else
			{
				pEffect1 = this->pEndowmentFunction()->rEffects()[effect1];
			}
			if (effect2 <= totalEvaluationEffects)
			{
				pEffect2 = this->pEvaluationFunction()->rEffects()[effect2];
			}
			else
			{
				pEffect2 = this->pEndowmentFunction()->rEffects()[effect2];
			}
			this->pSimulation()->derivative(pEffect1->pEffectInfo(),
				pEffect2->pEffectInfo(),
				this->pSimulation()->derivative(pEffect1->pEffectInfo(),
					pEffect2->pEffectInfo()) + product[effect1] *
				product[effect2]);
		}
	}
	delete [] product;
}

//	Rprintf("deriv %f\n", derivative;



/**
 * Returns whether applying the given ministep on the current state of this
 * variable would be valid with respect to all constraints.
 */
bool NetworkVariable::validMiniStep(const MiniStep * pMiniStep) const
{
	bool valid = DependentVariable::validMiniStep(pMiniStep);

	if (valid && !pMiniStep->diagonal())
	{
		NetworkLongitudinalData * pData =
			(NetworkLongitudinalData *) this->pData();
		const NetworkChange * pNetworkChange =
			dynamic_cast<const NetworkChange *>(pMiniStep);
		int i = pNetworkChange->ego();
		int j = pNetworkChange->alter();

		if (this->lpNetwork->tieValue(i, j))
		{
			valid = !pData->upOnly(this->period());
		}
		else
		{
			valid = !pData->downOnly(this->period()) &&
				this->lpNetwork->outDegree(i) < pData->maxDegree() &&
				this->pReceivers()->active(j);
		}

		if (valid)
		{
			valid = !pData->structural(i, j, this->period());
		}

		// The filters may add some more conditions.

		for (unsigned i = 0;
			i < this->lpermittedChangeFilters.size() && valid;
			i++)
		{
			PermittedChangeFilter * pFilter = this->lpermittedChangeFilters[i];
			valid = pFilter->validMiniStep(pNetworkChange);
		}
	}

	return valid;
}


/**
 * Generates a random ministep for the given ego.
 */
MiniStep * NetworkVariable::randomMiniStep(int ego)
{
	this->pSimulation()->pCache()->initialize(ego);
	this->lego = ego;
	this->calculateTieFlipProbabilities();
	int alter = nextIntWithProbabilities(this->m(), this->lprobabilities);

	MiniStep * pMiniStep =
		new NetworkChange(this->lpData, ego, alter);
	pMiniStep->logChoiceProbability(log(this->lprobabilities[alter]));

	return pMiniStep;
}


/**
 * Returns if the observed value for the option of the given ministep
 * is missing at either end of the period.
 */
bool NetworkVariable::missing(const MiniStep * pMiniStep) const
{
	const NetworkChange * pNetworkChange =
		dynamic_cast<const NetworkChange *>(pMiniStep);

	return this->lpData->missing(pNetworkChange->ego(),
			pNetworkChange->alter(),
			this->period()) ||
		this->lpData->missing(pNetworkChange->ego(),
			pNetworkChange->alter(),
			this->period() + 1);
}

/**
 * Returns if the given ministep is structurally determined in the period.
 */
bool NetworkVariable::structural(const MiniStep * pMiniStep) const
{
	const NetworkChange * pNetworkChange =
		dynamic_cast<const NetworkChange *>(pMiniStep);
	return this->lpData->structural(pNetworkChange->ego(),
			pNetworkChange->alter(),
			this->period());
}
/**
 * Copies the change contributions for evaluation and endowment function
 * effects according to the current miniStep.
 */
void NetworkVariable::copyChangeContributions(MiniStep * pMiniStep) const
{
	 NetworkChange * pNetworkChange =
		dynamic_cast< NetworkChange *>(pMiniStep);

	 int nEvaluationEffects = this->pEvaluationFunction()->rEffects().size();
	 int nEndowmentEffects = this->pEndowmentFunction()->rEffects().size();
	 pNetworkChange->allocateEffectContributionArrays(nEvaluationEffects,
		 nEndowmentEffects, this->m());

	 for (int alter = 0; alter < this->m(); alter++)
	 {
		for (unsigned i = 0;
			 i < this->pEvaluationFunction()->rEffects().size(); i++)
		{
			pNetworkChange->evaluationEffectContribution(
				this->levaluationEffectContribution[alter][i], alter, i);
		}

		for (unsigned i = 0;
			 i < this->pEndowmentFunction()->rEffects().size(); i++)
		{
			pNetworkChange->endowmentEffectContribution(
				this->lendowmentEffectContribution[alter][i], alter, i);
		}

	}
}
/**
 * Calculates the log probability of the choice of this ministep,
 * using stored change contributions.
 *
 */
double NetworkVariable::calculateChoiceProbability(const MiniStep * pMiniStep)
const
{
	const NetworkChange * pNetworkChange =
		dynamic_cast< const NetworkChange *>(pMiniStep);
	int evaluationEffectCount = this->pEvaluationFunction()->rEffects().size();
	int endowmentEffectCount = this->pEndowmentFunction()->rEffects().size();

	double total = 0;
	double * probabilities = new double[this->m()];
	double value;

	for (int alter = 0; alter < this->m(); alter++)
	{
		double contribution = 0;

		for (int i = 0; i < evaluationEffectCount; i++)
		{
			Effect * pEffect = this->pEvaluationFunction()->rEffects()[i];
			contribution += pEffect->parameter() *
				pNetworkChange->evaluationEffectContribution(alter, i);
			//	Rprintf("%d %d %d %f %f \n",alter, i, pNetworkChange->ego(),
			//	pEffect->parameter(),
			//	pNetworkChange->evaluationEffectContribution(alter, i));
		}

		for (int i = 0; i < endowmentEffectCount; i++)
		{
			Effect * pEffect = this->pEndowmentFunction()->rEffects()[i];
			contribution +=	pEffect->parameter() *
				pNetworkChange->endowmentEffectContribution(alter, i);
		}

		// The selection probability is the exponential of the total
		// contribution.

		probabilities[alter] = exp(contribution);
		if (R_IsNaN(probabilities[alter]))
		{
			PrintValue(getMiniStepDF(*pNetworkChange));
		}
		total += probabilities[alter];

	}

	// Normalize

	if (total != 0)
	{
		for (int alter = 0; alter < this->m(); alter++)
		{
			probabilities[alter] /= total;
		}
	}
	value = log(probabilities[pNetworkChange->alter()]);
	delete[] probabilities;
	return value;
}


// ----------------------------------------------------------------------------
// Section: Properties
// ----------------------------------------------------------------------------

/**
 * Returns if this is a network variable.
 */
bool NetworkVariable::networkVariable() const
{
	return true;
}

}
