/******************************************************************************
 * SIENA: Simulation Investigation for Empirical Network Analysis
 *
 * Web: http://www.stats.ox.ac.uk/~snijders/siena/
 *
 * File: BalanceEffect.cpp
 *
 * Description: This file contains the implementation of the class
 * BalanceEffect.
 *****************************************************************************/

#include <stdexcept>
#include "BalanceEffect.h"
#include "data/OneModeNetworkLongitudinalData.h"
#include "data/Network.h"
#include "data/TieIterator.h"
#include "model/variables/NetworkVariable.h"
#include "model/EffectInfo.h"
#include "model/tables/ConfigurationTable.h"

namespace siena
{

/**
 * Constructor.
 */
BalanceEffect::BalanceEffect(const EffectInfo * pEffectInfo) :
	NetworkEffect(pEffectInfo)
{
	this->lbalanceMean = 0;
}


/**
 * Initializes this effect for the use with the given epoch simulation.
 */
void BalanceEffect::initialize(EpochSimulation * pSimulation)
{
	NetworkEffect::initialize(pSimulation);

	const OneModeNetworkLongitudinalData * pData =
		dynamic_cast<const OneModeNetworkLongitudinalData *>(this->pData());

	if (pData)
	{
		this->lbalanceMean = pData->balanceMean();
	}
	else
	{
		throw logic_error("One-mode network variable '" +
			this->pEffectInfo()->variableName() +
			"' expected.");
	}
}


/**
 * Initializes this effect for calculating the corresponding statistics.
 * @param[in] pData the observed data
 * @param[in] pState the current state of the dependent variables
 * @param[in] period the period of interest
 */
void BalanceEffect::initialize(const Data * pData, State * pState, int period)
{
	NetworkEffect::initialize(pData, pState, period);

	const OneModeNetworkLongitudinalData * pNetworkData =
		dynamic_cast<const OneModeNetworkLongitudinalData *>(this->pData());

	if (pNetworkData)
	{
		this->lbalanceMean = pNetworkData->balanceMean();
	}
	else
	{
		throw logic_error("Data for one-mode network variable '" +
			this->pEffectInfo()->variableName() +
			"' expected.");
	}
}


/**
 * Calculates the contribution of a tie flip to the given actor.
 */
double BalanceEffect::calculateTieFlipContribution(int alter) const
{
	// The formula from SIENA manual:
	// s_i(x) = sum_j x_{ij} sum_{h!=i,j} (b_0 - |x_{ih} - x_{jh}|)
	// Rearrange and obtain s_i(x) = A - B, where
	//    A = x_{i+} (n-2) b_0 and
	//    B = sum_j x_{ij} sum_{h!=i,j} |x_{ih} - x_{jh}|
	// A is easy to compute.
	// B equals the number of non-transitive two-paths starting at i
	// plus the number of out-stars <(i,j),(i,h)> with no tie from
	// j to h.

	// The contribution of a tie flip to A is (n-2) b_0.
	double a = (this->pVariable()->n() - 2) * this->lbalanceMean;

	// These will be used later.

	int twoPathCount = this->pVariable()->pTwoPathTable()->get(alter);
	int inStarCount = this->pVariable()->pInStarTable()->get(alter);

	// First consider how the number of non-transitive two-paths would
	// change after introducing the tie from the ego i to the alter j.

	// x_{j+} - x_{ji} - IS_{ij} non-transitive two-paths would be created.

	double b =
		this->pVariable()->pNetwork()->outDegree(alter) - inStarCount;

	if (this->pVariable()->inTieExists(alter))
	{
		b--;
	}

	// However, TP_{ij} non-transitive two-paths would be lost as the
	// tie (i,j) would make them transitive.

	b -= twoPathCount;

	// The number of ties from the ego to actors other than the alter.

	int outDegree =
		this->pVariable()->pNetwork()->outDegree(this->pVariable()->ego());

	if (this->pVariable()->outTieExists(alter))
	{
		outDegree--;
	}

	// Now, consider the change in the number of out-stars <(i,j),(i,h)> with
	// no tie from j to h.

	// The introduction of the tie from the ego to the alter would create
	// x_{i+} - IS_{ij} new such out-stars, where the alter assumes the
	// role of j, and x_{i+} - TP_{ij} new out-stars, where the alter assumes
	// the role of h.

	b += 2 * outDegree - inStarCount - twoPathCount;

	double change = a - b;

	if (this->pVariable()->outTieExists(alter))
	{
		change = -change;
	}

	return change;
}


/**
 * Returns if the given configuration table is used by this effect
 * during the calculation of tie flip contributions.
 */
bool BalanceEffect::usesTable(const ConfigurationTable * pTable) const
{
	return pTable == this->pVariable()->pTwoPathTable() ||
		pTable == this->pVariable()->pInStarTable();
}


/**
 * Detailed comment in the base class.
 */
double BalanceEffect::statistic(Network * pNetwork,
	Network * pSummationTieNetwork) const
{
	double statistic = 0;
	int n = pNetwork->n();
	const Network * pStartMissingNetwork =
		this->pData()->pMissingTieNetwork(this->period());
	const Network * pEndMissingNetwork =
		this->pData()->pMissingTieNetwork(this->period() + 1);

	// An indicator array for invalid actors.
	// Invariants:
	// A: flag[i] <= round for all actors
	// B: flag[i] == round for invalid actors

	int * flag = new int[n];
	int round = 0;

	for (int i = 0; i < n; i++)
	{
		flag[i] = 0;
	}

	for (TieIterator iter = pSummationTieNetwork->ties();
		iter.valid();
		iter.next())
	{
		// Initially all actors are valid

		round++;
		int validActorCount = n;

		// Mark as invalid the actors h that have a missing tie from the ego or
		// the alter of the current tie at either end of the period.

		this->markInvalidActors(pStartMissingNetwork->outTies(iter.ego()),
			flag,
			round,
			validActorCount);
		this->markInvalidActors(pStartMissingNetwork->outTies(iter.alter()),
			flag,
			round,
			validActorCount);
		this->markInvalidActors(pEndMissingNetwork->outTies(iter.ego()),
			flag,
			round,
			validActorCount);
		this->markInvalidActors(pEndMissingNetwork->outTies(iter.alter()),
			flag,
			round,
			validActorCount);

		// Mark the ego and alter invalid as well.

		if (flag[iter.ego()] < round)
		{
			flag[iter.ego()] = round;
			validActorCount--;
		}

		if (flag[iter.alter()] < round)
		{
			flag[iter.alter()] = round;
			validActorCount--;
		}

		// Now we add the expression
		//   sum_h (b_0 - |x_{ih} - x_{jh}|)
		// to the statistic, where the sum extends over all valid actors h.

		// First add sum_h b_0
		statistic += validActorCount * this->lbalanceMean;

		// Now subtract sum_h |x_{ih} - x_{jh}|

		IncidentTieIterator egoIter = pNetwork->outTies(iter.ego());
		IncidentTieIterator alterIter = pNetwork->outTies(iter.alter());

		while (egoIter.valid() || alterIter.valid())
		{
			if (egoIter.valid() &&
				(!alterIter.valid() || egoIter.actor() < alterIter.actor()))
			{
				if (flag[egoIter.actor()] < round)
				{
					statistic--;
				}

				egoIter.next();
			}
			else if (alterIter.valid() &&
				(!egoIter.valid() || alterIter.actor() < egoIter.actor()))
			{
				if (flag[alterIter.actor()] < round)
				{
					statistic--;
				}

				alterIter.next();
			}
			else
			{
				egoIter.next();
				alterIter.next();
			}
		}
	}

	delete[] flag;
	return statistic;
}


/**
 * This method marks as invalid all actors that are iterated over by the given
 * iterator. The fact that an actor i is invalid is represented by setting
 * flag[i] = markValue. It is assumed that flag[i] <= markValue for all actors,
 * meaning that flag[i] < markValue holds for valid actors. The variable
 * validActorCount keeps track of the still valid actors.
 */
void BalanceEffect::markInvalidActors(IncidentTieIterator iter,
	int * flag,
	int markValue,
	int & validActorCount) const
{
	while (iter.valid())
	{
		if (flag[iter.actor()] < markValue)
		{
			flag[iter.actor()] = markValue;
			validActorCount--;
		}

		iter.next();
	}
}

}
