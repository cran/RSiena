/******************************************************************************
 * SIENA: Simulation Investigation for Empirical Network Analysis
 *
 * Web: http://www.stats.ox.ac.uk/~snijders/siena/
 *
 * File: TwoPathTable.cpp
 *
 * Description: This file contains the implementation of the TwoPathTable
 * class.
 *****************************************************************************/

#include "TwoPathTable.h"
#include "data/IncidentTieIterator.h"
#include "data/CommonNeighborIterator.h"
#include "data/OneModeNetwork.h"
#include "model/variables/NetworkVariable.h"

namespace siena
{

// ----------------------------------------------------------------------------
// Section: Initialization
// ----------------------------------------------------------------------------

/**
 * Creates a new table for storing two-paths with the specified directions for
 * the first and the second tie.
 */
TwoPathTable::TwoPathTable(NetworkVariable * pVariable,
	Direction firstStepDirection,
	Direction secondStepDirection) : ConfigurationTable(pVariable)
{
	this->lfirstStepDirection = firstStepDirection;
	this->lsecondStepDirection = secondStepDirection;
}


// ----------------------------------------------------------------------------
// Section: ConfigurationTable implementation
// ----------------------------------------------------------------------------

/**
 * Calculates the number of generalized two-paths between the ego and all
 * other actors.
 */
void TwoPathTable::vCalculate()
{
	// Reset the counters to zeroes
	this->reset();

	// Choose between the outgoing and incoming or reciprocated ties of the
	// ego, depending on the direction of the first step.

	if (this->lfirstStepDirection == FORWARD)
	{
		this->performFirstStep(
			this->pVariable()->pNetwork()->outTies(this->pVariable()->ego()));
	}
	else if (this->lfirstStepDirection == BACKWARD)
	{
		this->performFirstStep(
			this->pVariable()->pNetwork()->inTies(this->pVariable()->ego()));
	}
	else
	{
		OneModeNetwork * pOneModeNetwork =
			dynamic_cast<OneModeNetwork *>(this->pVariable()->pNetwork());

		this->performFirstStep(
			pOneModeNetwork->reciprocatedTies(this->pVariable()->ego()));
	}
}


/**
 * Performs the first step by iterating over the actors of the given
 * iterator and invoking the method for the second step.
 */
template<class Iterator>
void TwoPathTable::performFirstStep(Iterator iter)
{
	// TODO: Using templates here is a bad design. It's because
	// IncidentTieIterator and CommonNeighborIterator have no
	// common base class.

	// Try out all possible first steps

	while (iter.valid())
	{
		int middleActor = iter.actor();
		iter.next();

		// Choose the right iterator for the second step

		if (this->lsecondStepDirection == FORWARD)
		{
			this->performSecondStep(
				this->pVariable()->pNetwork()->outTies(middleActor));
		}
		else if (this->lsecondStepDirection == BACKWARD)
		{
			this->performSecondStep(
				this->pVariable()->pNetwork()->inTies(middleActor));
		}
		else
		{
			OneModeNetwork * pOneModeNetwork =
				dynamic_cast<OneModeNetwork *>(this->pVariable()->pNetwork());
			this->performSecondStep(
				pOneModeNetwork->reciprocatedTies(middleActor));
		}
	}
}


/**
 * Performs the second step by iterating over the actors of the given
 * iterator and incrementing their values that are stored in this table.
 */
template<class Iterator>
void TwoPathTable::performSecondStep(Iterator iter)
{
	while (iter.valid())
	{
		this->set(iter.actor(), this->get(iter.actor()) + 1);
		iter.next();
	}
}

}
