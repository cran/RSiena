/******************************************************************************
 * SIENA: Simulation Investigation for Empirical Network Analysis
 *
 * Web: http://www.stats.ox.ac.uk/~snijders/siena/
 *
 * File: IncidentTieIterator.h
 *
 * Description: This module defines the class IncidentTieIterator for
 * convenient iteration over incoming or outgoing ties of an actor.
 *****************************************************************************/

#ifndef INCIDENTTIEITERATOR_H_
#define INCIDENTTIEITERATOR_H_

#include <map>

namespace siena
{

/**
 * This class defines an iterator over incoming or outgoing ties of a specific
 * actor <i>i</i>. The ties are sorted in an increasing order of the neighbors
 * of <i>i</i>.
 */
class IncidentTieIterator
{
	// The class Network needs access to the private constructor.
	friend class Network;

public:
	IncidentTieIterator();

	int actor() const;
	int value() const;
	bool valid() const;
	void next();

private:
	IncidentTieIterator(std::map<int, int> & ties);
	IncidentTieIterator(std::map<int, int> & ties, int lowerBound);

	// Points to the current element in the underlying map
	std::map<int, int>::const_iterator lcurrent;

	// Points to the end of the underlying map
	std::map<int, int>::const_iterator lend;

	// Indicates if the internal iterators have been initialized to
	// point to elements of some underlying map.

	bool linitialized;
};

}

#endif /*INCIDENTTIEITERATOR_H_*/
