/******************************************************************************
 * SIENA: Simulation Investigation for Empirical Network Analysis
 *
 * Web: http://www.stats.ox.ac.uk/~snijders/siena/
 *
 * File: BetweennessEffect.cpp
 *
 * Description: This file contains the implementation of the class
 * BetweennessEffect.
 *****************************************************************************/

#include "BetweennessEffect.h"
#include "data/Network.h"
#include "data/IncidentTieIterator.h"
#include "data/NetworkLongitudinalData.h"
#include "model/variables/NetworkVariable.h"
#include "model/tables/ConfigurationTable.h"

namespace siena
{

/**
 * Constructor.
 */
BetweennessEffect::BetweennessEffect(const EffectInfo * pEffectInfo) :
	NetworkEffect(pEffectInfo)
{
}


/**
 * Calculates the contribution of a tie flip to the given actor.
 */
double BetweennessEffect::calculateTieFlipContribution(int alter) const
{
	// Get the number of actors h other than the alter j that have a tie
	// to the ego i.

	int inDegree =
		this->pVariable()->pNetwork()->inDegree(this->pVariable()->ego());

	if (this->pVariable()->inTieExists(alter))
	{
		inDegree--;
	}

	// Now, for each of these actors h, the introduction of the tie (i,j)
	// creates a new non-transitive two-path through i unless there's a tie
	// (h,j), in which case <(h,i),(h,j)> is an out-star.

	int change = inDegree - this->pVariable()->pOutStarTable()->get(alter);

	// If we are withdrawing the tie, the contribution is negative.

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
bool BetweennessEffect::usesTable(const ConfigurationTable * pTable) const
{
	return pTable == this->pVariable()->pOutStarTable();
}


/**
 * Detailed comment in the base class.
 */
double BetweennessEffect::statistic(Network * pNetwork,
	Network * pSummationTieNetwork) const
{
	int statistic = 0;
	int n = pNetwork->n();
	const Network * pStartMissingNetwork =
		this->pData()->pMissingTieNetwork(this->period());
	const Network * pEndMissingNetwork =
		this->pData()->pMissingTieNetwork(this->period() + 1);

	// A helper array of marks

	int * mark = new int[n];
	int currentMark = 0;

	for (int i = 0; i < n; i++)
	{
		mark[i] = 0;
	}

	// Count the number of non-transitive two-paths via each actor i in turn

	for (int i = 0; i < n; i++)
	{
		// Mark the in-neighbors of i by ensuring the following assertion:
		// mark[h] >= currentMark <==> the tie (h,i) exists

		currentMark++;

		for (IncidentTieIterator iter = pNetwork->inTies(i);
			iter.valid();
			iter.next())
		{
			mark[iter.actor()] = currentMark;
		}

		// Remember the current mark such that mark[h] >= baseMark if and only
		// if there's a tie from h to i.

		int baseMark = currentMark;

		// Now go through the ties (i,j) of the summation network and sum up
		// the numbers of non-transitive two-paths through i terminating at j.

		for (IncidentTieIterator iterJ = pSummationTieNetwork->outTies(i);
			iterJ.valid();
			iterJ.next())
		{
			int j = iterJ.actor();

			// For the start, assume that each in-neighbor of i gives one such
			// a two-path.

			statistic += pNetwork->inDegree(i);

			// The actor j itself doesn't count

			if (mark[j] >= baseMark)
			{
				statistic--;
			}

			// There are three cases where an actor h shoudln't contribute:
			// - There is a tie (h,j).
			// - The tie (h,j) is missing at the start of the period.
			// - The tie (h,j) is missing at the end of the period.
			// We should decrement the statistic for each of these actors h,
			// but we shouldn't do that more than once per actor, so we use
			// another mark.

			currentMark++;

			// Ties (h,j)

			for (IncidentTieIterator iterH = pNetwork->inTies(j);
				iterH.valid();
				iterH.next())
			{
				int h = iterH.actor();

				if (mark[h] >= baseMark && mark[h] < currentMark)
				{
					statistic--;
					mark[h] = currentMark;
				}
			}

			// Missing ties (h,j) at the start of the period

			for (IncidentTieIterator iterH = pStartMissingNetwork->inTies(j);
				iterH.valid();
				iterH.next())
			{
				int h = iterH.actor();

				if (mark[h] >= baseMark && mark[h] < currentMark)
				{
					statistic--;
					mark[h] = currentMark;
				}
			}

			// Missing ties (h,j) at the end of the period

			for (IncidentTieIterator iterH = pEndMissingNetwork->inTies(j);
				iterH.valid();
				iterH.next())
			{
				int h = iterH.actor();

				if (mark[h] >= baseMark && mark[h] < currentMark)
				{
					statistic--;
					mark[h] = currentMark;
				}
			}
		}
	}

	delete[] mark;
	return statistic;
}

}
