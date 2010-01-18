#include "CommonNeighborIterator.h"
#include "../utils/Utils.h"

namespace siena
{

// ----------------------------------------------------------------------------
// Section: Construction area
// ----------------------------------------------------------------------------

/**
 * Creates an iterator over actors common to both of the given incident
 * tie iterators.
 */
CommonNeighborIterator::CommonNeighborIterator(IncidentTieIterator iter1,
	IncidentTieIterator iter2)
{
	// Store the iterators.
	
	this->liter1 = iter1;
	this->liter2 = iter2;
	
	// Make sure they point to the first common actor.
	this->skipMismatches();
}


// ----------------------------------------------------------------------------
// Section: Public interface
// ----------------------------------------------------------------------------

/**
 * Indicates if there are still some common actors to be reported.
 */
bool CommonNeighborIterator::valid() const
{
	return this->liter1.valid() && this->liter2.valid();
}


/**
 * Returns the current common actor.
 */
int CommonNeighborIterator::actor() const
{
	if (!this->valid())
	{
		throw InvalidIteratorException();
	}

	// Both iterators point to the same actor, so we can use any iterator.
	return this->liter1.actor();
}


/**
 * Moves on to the next common actor if any.
 */
void CommonNeighborIterator::next()
{
	if (!this->valid())
	{
		throw InvalidIteratorException();
	}

	// Advance both iterators until they point to the next common actor
	// or we run out of actors.
	
	this->liter1.next();
	this->liter2.next();
	this->skipMismatches();
}


// ----------------------------------------------------------------------------
// Section: Private methods
// ----------------------------------------------------------------------------

/**
 * Advances both iterators until they point to the same actor. If there is
 * no such an actor, one of the iterators becomes invalid and we stop.
 */
void CommonNeighborIterator::skipMismatches()
{
	while (this->liter1.valid() &&
		this->liter2.valid() &&
		this->liter1.actor() != this->liter2.actor())
	{
		if (this->liter1.actor() < this->liter2.actor())
		{
			this->liter1.next();
		}
		else
		{
			this->liter2.next();
		}
	}
}

}
