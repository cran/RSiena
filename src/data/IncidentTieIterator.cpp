/******************************************************************************
 * SIENA: Simulation Investigation for Empirical Network Analysis
 * 
 * Web: http://www.stats.ox.ac.uk/~snijders/siena/
 * 
 * File: IncidentTieIterator.cpp
 * 
 * Description: This module defines the class IncidentTieIterator for
 * convenient iteration over incoming or outgoing ties of an actor.
 *****************************************************************************/

#include "IncidentTieIterator.h"
#include "../utils/Utils.h"

namespace siena
{

/**
 * Creates a dummy iterator with no underlying collection of ties.
 */
IncidentTieIterator::IncidentTieIterator()
{
	// The internal iterators lcurrent and lend are not initialized,
	// and consequently this iterator is not valid.
	
	this->linitialized = false;
}


//
// Creates an iterator over a collection of ties represented by the
// given map. The values of the pairs in the map represent the values
// of ties, and the keys represent the corresponding neighbors.
//
IncidentTieIterator::IncidentTieIterator(std::map<int, int> & ties)
{
	this->lcurrent = ties.begin();
	this->lend = ties.end();
	this->linitialized = true;
}


/**
 * Returns the neighbor incident to the current tie.
 */
int IncidentTieIterator::actor() const
{
	if (!this->valid())
	{
		throw InvalidIteratorException();
	}
	
	return this->lcurrent->first;
}


/**
 * Returns the value of the current tie.
 */
int IncidentTieIterator::value() const
{
	if (!this->valid())
	{
		throw InvalidIteratorException();
	}
	
	return this->lcurrent->second;	
}


/**
 * Indicates if the iterator still points to a valid tie.
 */
bool IncidentTieIterator::valid() const
{
	return this->linitialized && this->lcurrent != this->lend;
}


/**
 * Moves the iterator to the next tie.
 */
void IncidentTieIterator::next()
{
	this->lcurrent++;
}

}
