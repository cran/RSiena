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
#include "model/variables/NetworkVariable.h"
#include "data/Network.h"

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

	// Choose between the outgoing and incoming ties of the ego, depending
	// on the direction of the first step.
	
	IncidentTieIterator iter;
	
	if (this->lfirstStepDirection == FORWARD)
	{
		iter =
			this->pVariable()->pNetwork()->outTies(this->pVariable()->ego());
	}
	else
	{
		iter = this->pVariable()->pNetwork()->inTies(this->pVariable()->ego());
	}
	
	// Try out all possible first steps

	while (iter.valid())
	{
		int middleActor = iter.actor();
		iter.next();
		
		// Choose the right iterator for the second step
		
		IncidentTieIterator iter1;
		
		if (this->lsecondStepDirection == FORWARD)
		{
			iter1 = this->pVariable()->pNetwork()->outTies(middleActor);
		}
		else
		{
			iter1 = this->pVariable()->pNetwork()->inTies(middleActor);
		}
		
		// Try out all possibilities for the second step and increase
		// the counters of each destination.
		
		while (iter1.valid())
		{
			this->set(iter1.actor(), this->get(iter1.actor()) + 1);
			iter1.next();
		}
	}
}

}
