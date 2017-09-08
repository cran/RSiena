/******************************************************************************
 * SIENA: Simulation Investigation for Empirical Network Analysis
 *
 * Web: http://www.stats.ox.ac.uk/~snijders/siena/
 *
 * File: UnionTieIterator.cpp
 *
 * Description: This module defines an iterator that only the union
 * of two iterators.
 *****************************************************************************/

#include "UnionTieIterator.h"

namespace siena {

/**
 * Constructor cloning the both input iterators and additionally
 * adding identifiers to both iterators
 * @param[in] iter1 The first iterator.
 * @param[in] iter2 The second iterator.
 */
UnionTieIterator::UnionTieIterator(const ITieIterator& iter1,
		const ITieIterator& iter2) :
		CombinedTieIterator(iter1, iter2) {
	// adds default id's to the iterators
	init(1, 2);
}

/**
 * Constructor cloning the both input iterators and additionally
 * adding the given identifiers to both iterators
 * @param[in] iter1 The first iterator.
 * @param[in] iter2 The second iterator.
 * @param[in] idIter1 The ID of iter1.
 * @param[in] idIter2 The ID of iter2.
 */
UnionTieIterator::UnionTieIterator(const ITieIterator& iter1,
		const ITieIterator& iter2, int idIter1, int idIter2) :
		CombinedTieIterator(iter1, iter2) {
	// adds the given id's to the iterators
	init(idIter1, idIter2);
}

/**
 * @copydoc CombinedTieIterator::~CombinedTieIterator()
 */
UnionTieIterator::~UnionTieIterator() {
}

/**
 * @copydoc ITieIterator::next()
 */
void UnionTieIterator::next() {
	// at least one of the iterators has to be valid
	if (!valid()) {
		return;
	}
	// if both iterators are valid move to the next position
	if (lpIter1->valid() && lpIter2->valid()) {
		if (lpIter1->actor() < lpIter2->actor()) {
			lpIter1->next();
		} else if (lpIter1->actor() > lpIter2->actor()) {
			lpIter2->next();
		} else {
			lpIter1->next();
			lpIter2->next();
		}
	} else if (lpIter1->valid()) {
		lpIter1->next();
	} else {
		lpIter2->next();
	}
}

/**
 * @copydoc ITieIterator::actor()
 */
int UnionTieIterator::actor() const {
	// if both iterators are valid return the actor with the lower value
	if (lpIter1->valid() && lpIter2->valid()) {
		if (lpIter1->actor() <= lpIter2->actor()) {
			return lpIter1->actor();
		} else {
			return lpIter2->actor();
		}
		// else return the actor of the valid iterator
	} else if (lpIter1->valid()) {
		return lpIter1->actor();
	} else if (lpIter2->valid()) {
		return lpIter2->actor();
	}
	throw InvalidIteratorException();
}

/**
 * @copydoc ITieIterator::valid()
 */
bool UnionTieIterator::valid() const {
	return lpIter1->valid() || lpIter2->valid();
}

/**
 * @copydoc ITieIterator::clone()
 */
UnionTieIterator * UnionTieIterator::clone() const {
	return new UnionTieIterator(*this);
}

/**
 * Returns the ID of the iterator not adjacent to the current actor.
 * Note that if the current actor is common this method returns the
 * ID of the first iterator.
 * @return The ID of the iterator not adjacent to the current node
 */
int UnionTieIterator::getInactiveIterID() {
	if (valid()) {
		if (!lpIter1->valid()) {
			return lIdIter1;
		}
		if (!lpIter2->valid() || lpIter1->actor() < lpIter2->actor()) {
			return lIdIter2;
		}
		return lIdIter1;
	}
	throw InvalidIteratorException();
}

/**
 * Returns the ID of the iterator adjacent to the current actor. Note
 * that if the current actor is common this will always return the
 * ID of the second iterator.
 * @return The ID of the iterator adjacent to the current node
 */
int UnionTieIterator::getActiveIterID() {
	if (getInactiveIterID() == lIdIter1) {
		return lIdIter2;
	}
	return lIdIter1;
}

/**
 * Tells whether the current actor is adjacent to both iterators
 * or not.
 * @returns <code>True</code> if both iterators are adjacent to
 * the current actor, <code>False</code> otherwise.
 */
bool UnionTieIterator::isCommonNeighbor() const {
	return lpIter1->valid() && lpIter2->valid() && isCommon();
}

/**
 * Copy Constructor.
 */
UnionTieIterator::UnionTieIterator(const UnionTieIterator& rhs) :
		CombinedTieIterator(rhs), //
		lIdIter1(rhs.lIdIter1), //
		lIdIter2(rhs.lIdIter2) {
}

/**
 * Intializes the ID's of both iterators.
 * @param[in] idIter1 The ID of the first iterator.
 * @param[in] idIter2 The ID of the second iterator.
 */
void UnionTieIterator::init(int idIter1, int idIter2) {
	lIdIter1 = idIter1;
	lIdIter2 = idIter2;
}

} /* namespace siena */
