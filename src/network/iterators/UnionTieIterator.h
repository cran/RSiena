/******************************************************************************
 * SIENA: Simulation Investigation for Empirical Network Analysis
 *
 * Web: http://www.stats.ox.ac.uk/~snijders/siena/
 *
 * File: UnionTieIterator.h
 *
 * Description: This module defines an iterator that only the union
 * of two iterators.
 *****************************************************************************/

#ifndef UNIONTIEITERATOR_H_
#define UNIONTIEITERATOR_H_

#include "CombinedTieIterator.h"

namespace siena {

// ----------------------------------------------------------------------------
// Section: UnionTieIterator class
// ----------------------------------------------------------------------------

class UnionTieIterator: public CombinedTieIterator {
public:

	UnionTieIterator(const ITieIterator& iter1, const ITieIterator& iter2);

	UnionTieIterator(const ITieIterator& iter1, const ITieIterator& iter2,
			int idIter1, int idIter2);

	virtual ~UnionTieIterator();

	virtual void next();

	int actor() const;

	bool valid() const;

	virtual UnionTieIterator * clone() const;

	int getInactiveIterID();

	int getActiveIterID();

	bool isCommonNeighbor() const;

protected:

	UnionTieIterator(const UnionTieIterator& rhs);

private:

	// The ID of the first iterator.
	int lIdIter1;

	// The ID of the second iterator.
	int lIdIter2;

	void init(int idIter1, int idIter2);

	// Disable assignment operator.
	UnionTieIterator& operator=(const UnionTieIterator&);
};

} /* namespace siena */

#endif /* UNIONTIEITERATOR_H_ */
