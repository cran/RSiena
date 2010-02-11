/******************************************************************************
 * SIENA: Simulation Investigation for Empirical Network Analysis
 *
 * Web: http://www.stats.ox.ac.uk/~snijders/siena/
 *
 * File: Chain.cpp
 *
 * Description: This file contains the implementation of the class Chain.
 *****************************************************************************/

#include <vector>
#include "Chain.h"
#include "utils/Random.h"
#include "model/ml/MiniStep.h"
#include "model/ml/BehaviorChange.h"
#include "model/ml/NetworkChange.h"
#include "network/Network.h"
#include "network/IncidentTieIterator.h"
#include "data/Data.h"
#include "data/BehaviorLongitudinalData.h"
#include "data/NetworkLongitudinalData.h"

namespace siena
{

/**
 * Creates an empty chain.
 */
Chain::Chain()
{
	this->lpFirst = new MiniStep(0, 0, 0);
	this->lpLast = new MiniStep(0, 0, 0);

	this->lpFirst->pNext(this->lpLast);
	this->lpLast->pPrevious(this->lpFirst);

	this->lperiod = -1;
}


/**
 * Deallocates this chain.
 */
Chain::~Chain()
{
	this->emptyChain();

	delete this->lpFirst;
	delete this->lpLast;
	this->lpFirst = 0;
	this->lpLast = 0;
}


/**
 * Removes all ministeps from this chain except for the dummy ministeps
 * at the ends of the chain.
 */
void Chain::emptyChain()
{
	MiniStep * pMiniStep = this->lpFirst->pNext();

	while (pMiniStep != this->lpLast)
	{
		MiniStep * pNextMiniStep = pMiniStep->pNext();
		delete pMiniStep;
		pMiniStep = pNextMiniStep;
	}

	this->lpFirst->pNext(this->lpLast);
	this->lpLast->pPrevious(this->lpFirst);
}


/**
 * Inserts a new ministep before a ministep that already belongs to this chain.
 */
void Chain::insertBefore(MiniStep * pNewMiniStep, MiniStep * pExistingMiniStep)
{
	MiniStep * pPreviousMiniStep = pExistingMiniStep->pPrevious();

	pNewMiniStep->pPrevious(pPreviousMiniStep);
	pPreviousMiniStep->pNext(pNewMiniStep);

	pNewMiniStep->pNext(pExistingMiniStep);
	pExistingMiniStep->pPrevious(pNewMiniStep);
}


/**
 * Removes the given ministep from this chain without deleting the ministep.
 */
void Chain::remove(MiniStep * pMiniStep)
{
	pMiniStep->pPrevious()->pNext(pMiniStep->pNext());
	pMiniStep->pNext()->pPrevious(pMiniStep->pPrevious());
	pMiniStep->pNext(0);
	pMiniStep->pPrevious(0);
}


/**
 * Generates a random chain connecting the start and end observations of the
 * given data object for the given period. The chain is simple in the sense
 * that no two ministeps cancel each other out.
 */
void Chain::connect(Data * pData, int period)
{
	this->emptyChain();
	this->lperiod = period;
	vector<MiniStep *> miniSteps;

	// Create the required ministeps

	for (unsigned i = 0; i < pData->rDependentVariableData().size(); i++)
	{
		LongitudinalData * pVariableData = pData->rDependentVariableData()[i];
		NetworkLongitudinalData * pNetworkData =
			dynamic_cast<NetworkLongitudinalData *>(pVariableData);
		BehaviorLongitudinalData * pBehaviorData =
			dynamic_cast<BehaviorLongitudinalData *>(pVariableData);

		if (pNetworkData)
		{
			const Network * pNetwork1 = pNetworkData->pNetwork(period);
			const Network * pNetwork2 = pNetworkData->pNetwork(period + 1);

			for (int i = 0; i < pNetwork1->n(); i++)
			{
				IncidentTieIterator iter1 = pNetwork1->outTies(i);
				IncidentTieIterator iter2 = pNetwork2->outTies(i);

				while (iter1.valid() || iter2.valid())
				{
					if (iter1.valid() &&
						(!iter2.valid() || iter1.actor() < iter2.actor()))
					{
						miniSteps.push_back(
							new NetworkChange(i,
								iter1.actor(),
								pNetworkData->name(),
								-1));
						iter1.next();
					}
					else if (iter2.valid() &&
						(!iter1.valid() || iter2.actor() < iter1.actor()))
					{
						miniSteps.push_back(
							new NetworkChange(i,
								iter2.actor(),
								pNetworkData->name(),
								1));
						iter2.next();
					}
					else
					{
						iter1.next();
						iter2.next();
					}
				}
			}
		}
		else if (pBehaviorData)
		{
			for (int i = 0; i < pBehaviorData->n(); i++)
			{
				int delta = pBehaviorData->value(period + 1, i) -
					pBehaviorData->value(period, i);
				int singleChange = 1;

				if (delta < 0)
				{
					delta = -delta;
					singleChange = -1;
				}

				for (int j = 0; j < delta; j++)
				{
					miniSteps.push_back(
						new BehaviorChange(i,
							pBehaviorData->name(),
							singleChange));
				}
			}
		}
	}

	// Randomize the ministeps

	for (unsigned i = 1; i < miniSteps.size(); i++)
	{
		int j = nextInt(i + 1);
		MiniStep * pTempMiniStep = miniSteps[i];
		miniSteps[i] = miniSteps[j];
		miniSteps[j] = pTempMiniStep;
	}

	// And finally add the ministeps to this chain

	for (unsigned i = 0; i < miniSteps.size(); i++)
	{
		this->insertBefore(miniSteps[i], this->lpLast);
	}
}


/**
 * Returns the period whose end observations are connected by this chain.
 */
int Chain::period() const
{
	return this->lperiod;
}


/**
 * Returns the first (dummy) ministep in this chain.
 */
MiniStep * Chain::pFirst() const
{
	return this->lpFirst;
}


/**
 * Returns the last (dummy) ministep in this chain.
 */
MiniStep * Chain::pLast() const
{
	return this->lpLast;
}

}
