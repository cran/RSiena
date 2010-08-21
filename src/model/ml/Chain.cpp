/******************************************************************************
 * SIENA: Simulation Investigation for Empirical Network Analysis
 *
 * Web: http://www.stats.ox.ac.uk/~snijders/siena/
 *
 * File: Chain.cpp
 *
 * Description: This file contains the implementation of the class Chain.
 *****************************************************************************/
#include <Rinternals.h>

#include <vector>
#include "Chain.h"
#include "utils/Random.h"
#include "model/ml/MiniStep.h"
#include "model/ml/BehaviorChange.h"
#include "model/ml/NetworkChange.h"
#include "model/ml/MLSimulation.h"
#include "network/Network.h"
#include "network/IncidentTieIterator.h"
#include "data/Data.h"
#include "data/BehaviorLongitudinalData.h"
#include "data/NetworkLongitudinalData.h"

namespace siena
{
SEXP getMiniStepDF(const MiniStep& miniStep);

// ----------------------------------------------------------------------------
// Section: Constructors, destructors, initializers
// ----------------------------------------------------------------------------

/**
 * Creates an empty chain.
 */
Chain::Chain(Data * pData)
{
	this->lpFirst = new MiniStep(0, 0);
	this->lpLast = new MiniStep(0, 0);

	this->lpFirst->pNext(this->lpLast);
	this->lpLast->pPrevious(this->lpFirst);

	this->resetOrderingKeys();

	this->lpData = pData;
	this->lperiod = -1;

	this->lminiSteps.push_back(this->lpLast);
	this->lpLast->index(0);

	this->lmu = 0;
	this->lsigma2 = 0;
}


/**
 * Deallocates this chain.
 */
Chain::~Chain()
{
	this->clear();

	delete this->lpFirst;
	delete this->lpLast;
	this->lpFirst = 0;
	this->lpLast = 0;
	this->lpData = 0;

	this->lminiSteps.clear();
}


// ----------------------------------------------------------------------------
// Section: Chain modifications
// ----------------------------------------------------------------------------

/**
 * Removes all ministeps from this chain except for the dummy ministeps
 * at the ends of the chain.
 */
void Chain::clear()
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

	this->lminiSteps.clear();
	this->lminiSteps.push_back(this->lpLast);
	this->lpLast->index(0);

	this->ldiagonalMiniSteps.clear();
	this->lccpMiniSteps.clear();
	this->lmissingNetworkMiniSteps.clear();
	this->lmissingBehaviorMiniSteps.clear();
	this->lfirstMiniStepPerOption.clear();

	this->lmu = 0;
	this->lsigma2 = 0;
}


/**
 * Inserts a new ministep before a ministep that already belongs to this chain.
 */
void Chain::insertBefore(MiniStep * pNewMiniStep, MiniStep * pExistingMiniStep)
{
	MiniStep * pPreviousMiniStep = pExistingMiniStep->pPrevious();

	pNewMiniStep->pChain(this);

	// Update the pointers to next and previous ministeps.

	pNewMiniStep->pPrevious(pPreviousMiniStep);
	pPreviousMiniStep->pNext(pNewMiniStep);

	pNewMiniStep->pNext(pExistingMiniStep);
	pExistingMiniStep->pPrevious(pNewMiniStep);

	// Set the ordering key.

	double leftKey = pPreviousMiniStep->orderingKey();
	double rightKey = pExistingMiniStep->orderingKey();
	double newKey = (leftKey + rightKey) / 2;

	if ((-leftKey + newKey) > 1e-10 && (- newKey +rightKey) > 1e-10)
	{
		pNewMiniStep->orderingKey(newKey);
	}
	else
	{
		// Bad. We have exhausted the double precision.
		// However, using a linear time effort, we can restore a correct order.

		this->resetOrderingKeys();
	}

	// Update the array of diagonal ministeps

	if (pNewMiniStep->diagonal())
	{
		this->ldiagonalMiniSteps.push_back(pNewMiniStep);
		pNewMiniStep->diagonalIndex(this->ldiagonalMiniSteps.size() - 1);
	}

	if (pNewMiniStep->missing(this->lperiod))
	{
		if (pNewMiniStep->networkMiniStep())
		{
			pNewMiniStep->missingIndex(this->lmissingNetworkMiniSteps.size());
			this->lmissingNetworkMiniSteps.push_back(pNewMiniStep);
		}
		else
		{
			pNewMiniStep->missingIndex(this->lmissingBehaviorMiniSteps.size());
			this->lmissingBehaviorMiniSteps.push_back(pNewMiniStep);
		}
	}

	// Update the array of all ministeps

	this->lminiSteps.push_back(pNewMiniStep);
	pNewMiniStep->index(this->lminiSteps.size() - 1);

	this->updateSameOptionPointersOnInsert(pNewMiniStep);

	// Update the vector of CCPs for each ministep whose CCP status
	// might have changed.

	this->updateCCPs(pNewMiniStep->pPrevious());
	this->updateCCPs(pNewMiniStep->pPreviousWithSameOption());
	this->updateCCPs(pNewMiniStep);

	// Update some variables that are constantly maintained

	double rr = pNewMiniStep->reciprocalRate();
	this->lmu += rr;
	this->lsigma2 += rr * rr;
}


/**
 * Removes the given ministep from this chain without deleting the ministep.
 */
void Chain::remove(MiniStep * pMiniStep)
{
	MiniStep * pPrevious = pMiniStep->pPrevious();
	MiniStep * pPreviousWithSameOption = pMiniStep->pPreviousWithSameOption();

	// Updates pointers to the next and previous ministep.

	pPrevious->pNext(pMiniStep->pNext());
	pMiniStep->pNext()->pPrevious(pPrevious);
	pMiniStep->pNext(0);
	pMiniStep->pPrevious(0);

	// Next and previous pointers to ministeps of the same option

	if (this->lfirstMiniStepPerOption[*pMiniStep->pOption()] ==
		pMiniStep)
	{
		this->lfirstMiniStepPerOption[*pMiniStep->pOption()] =
			pMiniStep->pNextWithSameOption();
	}

	if (pMiniStep->pNextWithSameOption())
	{
		pMiniStep->pNextWithSameOption()->pPreviousWithSameOption(
			pPreviousWithSameOption);
	}

	if (pPreviousWithSameOption)
	{
		pPreviousWithSameOption->pNextWithSameOption(
			pMiniStep->pNextWithSameOption());
	}

	pMiniStep->pPreviousWithSameOption(0);
	pMiniStep->pNextWithSameOption(0);

	// Update the vector of CCPs

	this->updateCCPs(pPrevious);
	this->updateCCPs(pPreviousWithSameOption);
	this->updateCCPs(pMiniStep);

	// Update the vector of diagonal ministeps

	if (pMiniStep->diagonal())
	{
		MiniStep * pLastMiniStep =
			this->ldiagonalMiniSteps[this->ldiagonalMiniSteps.size() - 1];
		this->ldiagonalMiniSteps[pMiniStep->diagonalIndex()] = pLastMiniStep;
		pLastMiniStep->diagonalIndex(pMiniStep->diagonalIndex());
		this->ldiagonalMiniSteps.pop_back();
		pMiniStep->diagonalIndex(-1);
	}

	// Update the vectors of ministeps for missing options

	if (pMiniStep->missing(this->lperiod))
	{
		vector<MiniStep *> * pVector = &this->lmissingNetworkMiniSteps;

		if (pMiniStep->behaviorMiniStep())
		{
			pVector = &this->lmissingBehaviorMiniSteps;
		}

		MiniStep * pLast = (*pVector)[pVector->size() - 1];
		(*pVector)[pMiniStep->missingIndex()] = pLast;
		pLast->missingIndex(pMiniStep->missingIndex());
		pVector->pop_back();
		pMiniStep->missingIndex(-1);
	}

	MiniStep * pLastMiniStep = this->lminiSteps[this->lminiSteps.size() - 1];
	this->lminiSteps[pMiniStep->index()] = pLastMiniStep;
	pLastMiniStep->index(pMiniStep->index());
	this->lminiSteps.pop_back();
	pMiniStep->index(-1);

	double rr = pMiniStep->reciprocalRate();
	this->lmu -= rr;
	this->lsigma2 -= rr * rr;

	pMiniStep->pChain(0);
}


/**
 * Generates a random chain connecting the start and end observations of the
 * given data object for the given period. The chain is simple in the sense
 * that no two ministeps cancel each other out.
 */
void Chain::connect(int period)
{
	this->clear();
	this->lperiod = period;
	vector<MiniStep *> miniSteps;

	// Create the required ministeps

	for (unsigned variableIndex = 0;
		variableIndex < this->lpData->rDependentVariableData().size();
		variableIndex++)
	{
		LongitudinalData * pVariableData =
			this->lpData->rDependentVariableData()[variableIndex];
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
						if (!pNetworkData->structural(i, iter1.actor(),
								period)
							//	|| !pNetworkData->structural(i, iter1.actor(),
							//	period + 1)
							)
						{
						miniSteps.push_back(
							new NetworkChange(pNetworkData,
								i,
								iter1.actor()));
						iter1.next();
					}
						else
						{
							// create step in structural subchain?
						}
					}
					else if (iter2.valid() &&
						(!iter1.valid() || iter2.actor() < iter1.actor()))
					{
						if (!pNetworkData->structural(i, iter2.actor(),
								period)
							//	|| !pNetworkData->structural(i, iter2.actor(),
							// period + 1)
							)
						{
						miniSteps.push_back(
							new NetworkChange(pNetworkData,
								i,
								iter2.actor()));
						iter2.next();
					}
					else
					{
							// create step in structural subchain?
						}
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
					if (!pBehaviorData->structural(period, j)
						//|| !pBehaviorData->structural(period, j + 1)
						)

					{
					miniSteps.push_back(
						new BehaviorChange(pBehaviorData,
							i,
							singleChange));
					}
					else
					{
						// create step in structural subchain?
					}
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


// ----------------------------------------------------------------------------
// Section: Various updates
// ----------------------------------------------------------------------------

/**
 * Sets the value of the period of the chain.
 */
void Chain::period(int value)
{
	this->lperiod = value;
}

/**
 * Performs the necessary updates when the reciprocal rate of the
 * given ministep changes to the given value.
 */
void Chain::onReciprocalRateChange(const MiniStep * pMiniStep, double newValue)
{
	double oldValue = pMiniStep->reciprocalRate();

	this->lmu -= oldValue;
	this->lsigma2 -= oldValue * oldValue;

	this->lmu += newValue;
	this->lsigma2 += newValue * newValue;
}


/**
 * Assigns new ordering keys for the ministeps of this chain. The new
 * ordering key of a ministep is its index in the list.
 */
void Chain::resetOrderingKeys()
{
	int key = 0;
	MiniStep * pMiniStep = this->lpFirst;

	while (pMiniStep)
	{
		pMiniStep->orderingKey(key);
		pMiniStep = pMiniStep->pNext();
		key++;
	}
}


/**
 * Updates the "same option" related structures after the insertion of
 * the given ministep.
 */
void Chain::updateSameOptionPointersOnInsert(MiniStep * pNewMiniStep)
{
	MiniStep * pFirstMiniStep =
		this->lfirstMiniStepPerOption[*pNewMiniStep->pOption()];

	if (!pFirstMiniStep ||
		pFirstMiniStep->orderingKey() > pNewMiniStep->orderingKey())
	{
		// The new ministep is the earliest for its option

		pNewMiniStep->pNextWithSameOption(pFirstMiniStep);

		if (pFirstMiniStep)
		{
			pFirstMiniStep->pPreviousWithSameOption(pNewMiniStep);
		}

		this->lfirstMiniStepPerOption[*pNewMiniStep->pOption()] =
			pNewMiniStep;
	}
	else
	{
		MiniStep * pLeftMiniStep = pFirstMiniStep;

		while (pLeftMiniStep->pNextWithSameOption() &&
			pLeftMiniStep->pNextWithSameOption()->orderingKey() <
				pNewMiniStep->orderingKey())
		{
			pLeftMiniStep = pLeftMiniStep->pNextWithSameOption();
		}

		MiniStep * pRightMiniStep = pLeftMiniStep->pNextWithSameOption();

		pLeftMiniStep->pNextWithSameOption(pNewMiniStep);
		pNewMiniStep->pPreviousWithSameOption(pLeftMiniStep);

		if (pRightMiniStep)
		{
			pNewMiniStep->pNextWithSameOption(pRightMiniStep);
			pRightMiniStep->pPreviousWithSameOption(pNewMiniStep);
		}
	}
}


/**
 * Updates the vector of consecutive canceling pairs (CCPs) if the
 * CCP status of the given ministep has changed.
 */
void Chain::updateCCPs(MiniStep * pMiniStep)
{
	if (pMiniStep)
	{
		if (pMiniStep->firstOfConsecutiveCancelingPair() &&
			pMiniStep->consecutiveCancelingPairIndex() == -1)
		{
			pMiniStep->consecutiveCancelingPairIndex(
				this->lccpMiniSteps.size());
			this->lccpMiniSteps.push_back(pMiniStep);
		}
		else if (!pMiniStep->firstOfConsecutiveCancelingPair() &&
			pMiniStep->consecutiveCancelingPairIndex() != -1)
		{
			MiniStep * pLastCCPMiniStep =
				this->lccpMiniSteps[this->lccpMiniSteps.size() - 1];
			int index = pMiniStep->consecutiveCancelingPairIndex();
			this->lccpMiniSteps[index] = pLastCCPMiniStep;
			pLastCCPMiniStep->consecutiveCancelingPairIndex(index);
			this->lccpMiniSteps.pop_back();
			pMiniStep->consecutiveCancelingPairIndex(-1);
		}
	}
}


// ----------------------------------------------------------------------------
// Section: Accessors
// ----------------------------------------------------------------------------

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


/**
 * Returns the number of ministeps of this chain excluding the first,
 * but including the last dummy ministep.
 */
int Chain::ministepCount() const
{
	return this->lminiSteps.size();
}


/**
 * Returns the number of diagonal ministeps in this chain.
 */
int Chain::diagonalMinistepCount() const
{
	return this->ldiagonalMiniSteps.size();
}


/**
 * Returns the sum of reciprocal rates over all non-dummy ministeps.
 */
double Chain::mu() const
{
	return this->lmu;
}


/**
 * Returns the sum of squared reciprocal rates over all non-dummy ministeps.
 */
double Chain::sigma2() const
{
	return this->lsigma2;
}


/**
 * Returns the number of CCPs in this chain.
 */
int Chain::consecutiveCancelingPairCount() const
{
	return this->lccpMiniSteps.size();
}


/**
 * Returns the number of network ministeps with missing observed data
 * at either end of the current period.
 */
int Chain::missingNetworkMiniStepCount() const
{
	return this->lmissingNetworkMiniSteps.size();
}


/**
 * Returns the number of behavior ministeps with missing observed data
 * at either end of the current period.
 */
int Chain::missingBehaviorMiniStepCount() const
{
	return this->lmissingBehaviorMiniSteps.size();
}


// ----------------------------------------------------------------------------
// Section: Random draws
// ----------------------------------------------------------------------------

/**
 * Returns a random ministep excluding the first dummy ministep.
 */
MiniStep * Chain::randomMiniStep() const
{
	return this->lminiSteps[nextInt(this->lminiSteps.size())];
}


/**
 * Returns a random diagonal ministep.
 */
MiniStep * Chain::randomDiagonalMiniStep() const
{
	return this->ldiagonalMiniSteps[nextInt(this->ldiagonalMiniSteps.size())];
}


/**
 * Returns a random ministep from the given interval.
 */
MiniStep * Chain::randomMiniStep(MiniStep * pFirstMiniStep,
	MiniStep * pLastMiniStep) const
{
	int length = this->intervalLength(pFirstMiniStep, pLastMiniStep);
	int index = nextInt(length);
	MiniStep * pMiniStep = pFirstMiniStep;

	while (index > 0)
	{
		pMiniStep = pMiniStep->pNext();
		index--;
	}

	return pMiniStep;
}


/**
 * Returns the first mini step of a random consecutive canceling pair.
 */
MiniStep * Chain::randomConsecutiveCancelingPair() const
{
	return this->lccpMiniSteps[nextInt(this->lccpMiniSteps.size())];
}


/**
 * Returns a random network change with missing observed data.
 */
MiniStep * Chain::randomMissingNetworkMiniStep() const
{
	return this->lmissingNetworkMiniSteps[
		nextInt(this->lmissingNetworkMiniSteps.size())];
}


/**
 * Returns a random behavior change with missing observed data.
 */
MiniStep * Chain::randomMissingBehaviorMiniStep() const
{
	return this->lmissingBehaviorMiniSteps[
		nextInt(this->lmissingBehaviorMiniSteps.size())];
}

/**
 * Prints the list of ccps.
 */
void  Chain::printConsecutiveCancelingPairs() const
{
	// rr debug code for ccps
		Rprintf("\nStart\n ");
	for (unsigned i = 0; i < this->lccpMiniSteps.size(); i++)
	{
		PrintValue(getMiniStepDF(*this->lccpMiniSteps[i]));
	}
		Rprintf("\nend\n ");
}

// ----------------------------------------------------------------------------
// Section: Intervals
// ----------------------------------------------------------------------------

/**
 * Returns the length of the given interval of ministeps.
 */
int Chain::intervalLength(const MiniStep * pFirstMiniStep,
	const MiniStep * pLastMiniStep) const
{
	int length = 1;
	const MiniStep * pMiniStep = pFirstMiniStep;

	while (pMiniStep != pLastMiniStep)
	{
		pMiniStep = pMiniStep->pNext();
		length++;
	}

	return length;
}


// ----------------------------------------------------------------------------
// Section: Same option related
// ----------------------------------------------------------------------------

/**
 * Returns the first ministep of the given option in the subchain starting
 * with the given ministep, or 0, if there is no such a ministep.
 */
MiniStep * Chain::nextMiniStepForOption(const Option & rOption,
	const MiniStep * pFirstMiniStep) const
{
	MiniStep * pMiniStep = 0;
	map<const Option, MiniStep *>::const_iterator iter =
		this->lfirstMiniStepPerOption.find(rOption);

	if (iter != this->lfirstMiniStepPerOption.end())
	{
		pMiniStep = iter->second;

		while (pMiniStep &&
			pMiniStep->orderingKey() < pFirstMiniStep->orderingKey())
		{
			pMiniStep = pMiniStep->pNextWithSameOption();
		}
	}

	return pMiniStep;
}
// ----------------------------------------------------------------------------
// Section: copy a chain (too difficult (for Ruth!) to do a copy constructor!
// ----------------------------------------------------------------------------

/**
 * Returns a copy of the given chain
 */
Chain * Chain::copyChain()
{
	Chain * pChain = new Chain(this->lpData);

	pChain->lperiod = this->lperiod;

	for (unsigned i = 1; i < this->lminiSteps.size(); i++)
	{
		MiniStep *pMiniStep = this->lminiSteps[i]->createCopyMiniStep();

		pChain->insertBefore(pMiniStep, pChain->lpLast);
	}

	pChain->lmu = this->lmu;
	pChain->lsigma2 = this->lsigma2;
//	Rprintf("%x\n", pChain);
	return pChain;
}
}
