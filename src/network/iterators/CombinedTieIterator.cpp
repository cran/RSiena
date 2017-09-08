/******************************************************************************
 * SIENA: Simulation Investigation for Empirical Network Analysis
 *
 * Web: http://www.stats.ox.ac.uk/~snijders/siena/
 *
 * File: CombinedTieIterator.cpp
 *
 * Description: This abstract iterator allows to combine two ITieIterators.
 *****************************************************************************/

#include "CombinedTieIterator.h"

namespace siena {

/**
 * Constructor cloning the both input iterators.
 * @param[in] iter1 The first iterator.
 * @param[in] iter2 The second iterator.
 */
CombinedTieIterator::CombinedTieIterator(const ITieIterator& iter1,
		const ITieIterator& iter2) :
		ITieIterator(), //
		lpIter1(iter1.clone()), //
		lpIter2(iter2.clone()) {
}

/**
 * Destructor.
 */
CombinedTieIterator::~CombinedTieIterator() {
	delete lpIter1;
	delete lpIter2;
}

/**
 * Copy Constructor.
 */
CombinedTieIterator::CombinedTieIterator(const CombinedTieIterator& rhs) :
		ITieIterator(rhs), //
		lpIter1(rhs.lpIter1->clone()), //
		lpIter2(rhs.lpIter2->clone()) {
}

} /* namespace siena */
