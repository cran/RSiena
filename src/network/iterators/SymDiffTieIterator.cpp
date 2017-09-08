/******************************************************************************
 * SIENA: Simulation Investigation for Empirical Network Analysis
 *
 * Web: http://www.stats.ox.ac.uk/~snijders/siena/
 *
 * File: SymDiffTieIterator.cpp
 *
 * Description: This module defines an iterator that only accesses the
 * symmetric difference of two iterators.
 *****************************************************************************/

#include "SymDiffTieIterator.h"

namespace siena {

/**
 * @copydoc UnionTieIterator::UnionTieIterator(ITieIterator,ITieIterator)
 */
SymDiffTieIterator::SymDiffTieIterator(const ITieIterator& iter1,
		const ITieIterator& iter2) :
		UnionTieIterator(iter1, iter2) {
	// move the iterator to the first valid position
	init();
}

/**
 * @copydoc UnionTieIterator::UnionTieIterator(ITieIterator,ITieIterator,int,int)
 */
SymDiffTieIterator::SymDiffTieIterator(const ITieIterator& iter1,
		const ITieIterator& iter2, int idIter1, int idIter2) :
		UnionTieIterator(iter1, iter2, idIter1, idIter2) {
	// move the iterator to the first valid position
	init();
}

/**
 * @copydoc UnionTieIterator::~UnionTieIterator()
 */
SymDiffTieIterator::~SymDiffTieIterator() {
}

/**
 * @copydoc UnionTieIterator::clone()
 */
SymDiffTieIterator* SymDiffTieIterator::clone() const {
	return new SymDiffTieIterator(*this);
}

/**
 * Moves the iterator to the next position where both iterators
 * point to different actors.
 */
void SymDiffTieIterator::next() {
	do {
		UnionTieIterator::next();
	} while (lpIter1->valid() && lpIter2->valid() && isCommon());
}

/**
 * Copy Constructor.
 */
SymDiffTieIterator::SymDiffTieIterator(const SymDiffTieIterator& rhs) :
		UnionTieIterator(rhs) {
}

/**
 * Moves the iterator to the first position where both iterators
 * point to different actors.
 */
void SymDiffTieIterator::init() {
	if (lpIter1->valid() && lpIter2->valid() && isCommon()) {
		next();
	}
}

} /* namespace siena */
