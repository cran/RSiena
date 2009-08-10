/******************************************************************************
 * SIENA: Simulation Investigation for Empirical Network Analysis
 * 
 * Web: http://www.stats.ox.ac.uk/~snijders/siena/
 * 
 * File: CriticalInStarTable.cpp
 * 
 * Description: This file contains the implementation of the class
 * CriticalInStarTable.
 *****************************************************************************/

#include "CriticalInStarTable.h"
#include "model/variables/NetworkVariable.h"
#include "data/IncidentTieIterator.h"
#include "data/Network.h"

namespace siena
{

// ----------------------------------------------------------------------------
// Section: Initialization
// ----------------------------------------------------------------------------

/**
 * Creates a new table for storing critical in-stars.
 */
CriticalInStarTable::CriticalInStarTable(NetworkVariable * pVariable) :
	ConfigurationTable(pVariable)
{
}


// ----------------------------------------------------------------------------
// Section: ConfigurationTable implementation
// ----------------------------------------------------------------------------

/**
 * Calculates the number of critical in-stars between the ego and all
 * other actors.
 */
void CriticalInStarTable::vCalculate()
{
	// Reset the counters to zeroes
	this->reset();

	// We will need the two-path counts, so make sure they are calculated.
	this->pVariable()->pTwoPathTable()->calculate();
	
	const Network * pNetwork = this->pVariable()->pNetwork();
	
	// Consider each outgoing tie of the ego in turn.
	
	for (IncidentTieIterator iter =
			pNetwork->outTies(this->pVariable()->ego());
		iter.valid();
		iter.next())
	{
		// Get the receiver of the outgoing tie.
		int h = iter.actor();
		
		// If there are more than one two-paths between i and h, then no
		// in-star pointing to h can be critical.
		
		if (this->pVariable()->pTwoPathTable()->get(h) == 0)
		{
			// If there are no two-paths between i and h, then every in-star
			// pointing to h is critical. Just iterate over incoming ties of
			// h to reach them all.
			
			for (IncidentTieIterator iter1 = pNetwork->inTies(h);
				iter1.valid();
				iter1.next())
			{
				this->set(iter1.actor(), this->get(iter1.actor()) + 1);
			}
		}
		else if (this->pVariable()->pTwoPathTable()->get(h) == 1)
		{
			// If there is exactly one two-path between i and h, say
			// i -> j -> h, then <(i,h), (j,h)> is the only critical in-star
			// pointing to h. Just iterate over incoming ties of h to
			// find the actor j.

			// This can be optimized by storing one of the intermediary
			// actors between ego and each of the other actors when
			// counting two-paths. The the only intermediary actor could
			// be accessed in constant time here.
			
			bool found = false;
			
			for (IncidentTieIterator iter1 = pNetwork->inTies(h);
				iter1.valid() && !found;
				iter1.next())
			{
				int j = iter1.actor();
				
				if (this->pVariable()->outTieExists(j))
				{
					this->set(j, this->get(j) + 1);
					found = true;
				} 
			}			
		}
	}
}

}
