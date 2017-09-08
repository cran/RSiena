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

namespace siena {

/**
 * Creates a dummy iterator with no underlying collection of ties.
 */
IncidentTieIterator::IncidentTieIterator() :
		ITieIterator(), //
		lcurrent(), //
		lend(lcurrent) {
}


//
// Creates an iterator over a collection of ties represented by the
// given map. The values of the pairs in the map represent the values
// of ties, and the keys represent the corresponding neighbors.
//
IncidentTieIterator::IncidentTieIterator(const std::map<int, int> & ties) :
		ITieIterator(), //
		lcurrent(ties.begin()), //
		lend(ties.end()) {
}

IncidentTieIterator::IncidentTieIterator(const IncidentTieIterator& rhs) :
		ITieIterator(rhs), //
		lcurrent(rhs.lcurrent), //
		lend(rhs.lend) {
}

IncidentTieIterator* IncidentTieIterator::clone() const {
	return new IncidentTieIterator(*this);
}

//
// Creates an iterator over a collection of ties represented by the
// given map. The values of the pairs in the map represent the values
// of ties, and the keys represent the corresponding neighbors. Only
// neighbors that are greater or equal with the given bound are returned.
//
IncidentTieIterator::IncidentTieIterator(const std::map<int, int> & ties,
		int lowerBound) :
		ITieIterator(), //
		lcurrent(ties.lower_bound(lowerBound)), //
		lend(ties.end()) {
}

}
