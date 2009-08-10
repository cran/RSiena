/******************************************************************************
 * SIENA: Simulation Investigation for Empirical Network Analysis
 * 
 * Web: http://www.stats.ox.ac.uk/~snijders/siena/
 * 
 * File: ConfigurationTable.cpp
 * 
 * Description: This file contains the implementation of the ConfigurationTable
 * class.
 *****************************************************************************/

#include "ConfigurationTable.h"
#include "model/variables/NetworkVariable.h"
#include "data/Network.h"

namespace siena
{

// ----------------------------------------------------------------------------
// Section: Constructors, destructor, and initializers
// ----------------------------------------------------------------------------

/**
 * Creates an empty configuration table for the given network variable, which
 * becomes the owner of the table.
 */
ConfigurationTable::ConfigurationTable(NetworkVariable * pVariable)
{
	this->lpVariable = pVariable;
	this->lpTable = new int[pVariable->n()];
	this->lvalid = false;
}


/**
 * Deallocates this configuration table.
 */
ConfigurationTable::~ConfigurationTable()
{
	delete[] this->lpTable;
	this->lpTable = 0;
}


// ----------------------------------------------------------------------------
// Section: Public interface
// ----------------------------------------------------------------------------

/**
 * Calculates the table for the current network and ego of the owner variable.
 */
void ConfigurationTable::calculate()
{
	// Avoid dupplicate calculation
	
	if (!this->lvalid)
	{
		// Let the derived classes do the actual calculation.
		this->vCalculate();
		
		// Remember that the table has been calculated.
		this->lvalid = true;
	}
}


/**
 * Marks the table to be invalid such that it is recalculated by the next call
 * to calculate().
 */
void ConfigurationTable::invalidate()
{
	this->lvalid = false;
}


// ----------------------------------------------------------------------------
// Section: Accessors
// ----------------------------------------------------------------------------

/**
 * Returns the owner variable.
 */
NetworkVariable * ConfigurationTable::pVariable() const
{
	return this->lpVariable;
}


/**
 * Stores the number of configurations corresponding to the given actor.
 */
void ConfigurationTable::set(int i, int value)
{
	this->lpTable[i] = value;
}


// ----------------------------------------------------------------------------
// Section: Protected methods
// ----------------------------------------------------------------------------

/**
 * Resets the internal array to zeroes.
 */
void ConfigurationTable::reset()
{
	for (int i = 0; i < this->pVariable()->n(); i++)
	{
		this->lpTable[i] = 0;
	}
}

}
