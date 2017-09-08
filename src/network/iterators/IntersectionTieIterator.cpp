/******************************************************************************
 * SIENA: Simulation Investigation for Empirical Network Analysis
 *
 * Web: http://www.stats.ox.ac.uk/~snijders/siena/
 *
 * File: IntersectionTieIterator.cpp
 *
 * Description: This module defines an iterator that only accesses the
 * intersection of two iterators.
 *****************************************************************************/

#include "IntersectionTieIterator.h"

namespace siena {

/**
 * @copydoc CombinedTieIterator::CombinedTieIterator()
 * Constructor.
 */
IntersectionTieIterator::IntersectionTieIterator(const ITieIterator& iter1,
		const ITieIterator& iter2) :
		CombinedTieIterator(iter1, iter2) {
	if (valid() && !isCommon()) {
		skip();
	}
}

/**
 * @copydoc CombinedTieIterator::~CombinedTieIterator()
 */
IntersectionTieIterator::~IntersectionTieIterator() {
}

/**
 * @copydoc ITieIterator::actor()
 */
int IntersectionTieIterator::actor() const {
	if (valid()) {
		return lpIter1->actor();
	}
	throw InvalidIteratorException();
}

/**
 * Moves the iterator to the next common actor.
 * @sa CombinedTieIterator::isCommon()
 */
void IntersectionTieIterator::next() {
	lpIter1->next();
	lpIter2->next();
	skip();
}

/**
 * Moves the iterator until one of them becomes invalid or they
 * point to the same actor.
 */
void IntersectionTieIterator::skip() {
	while (valid() && !isCommon()) {
		while (lpIter1->valid() && lpIter1->actor() < lpIter2->actor()) {
			lpIter1->next();
		}
		if (!lpIter1->valid()) {
			return;
		}
		while (lpIter2->valid() && lpIter2->actor() < lpIter1->actor()) {
			lpIter2->next();
		}
	}
}

/**
 * @copydoc ITieIterator::valid()
 */
bool IntersectionTieIterator::valid() const {
	return lpIter1->valid() && lpIter2->valid();
}

/**
 * @copydoc ITieIterator::clone()
 */
IntersectionTieIterator* IntersectionTieIterator::clone() const {
	return new IntersectionTieIterator(*this);
}

/**
 * Copy Constructor.
 */
IntersectionTieIterator::IntersectionTieIterator(
		const IntersectionTieIterator& rhs) :
		CombinedTieIterator(rhs) {
	// note there is no need to skip as rhs is an IntersectionTieIterator itself
}

} /* namespace siena */
