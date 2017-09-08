/******************************************************************************
 * SIENA: Simulation Investigation for Empirical Network Analysis
 *
 * Web: http://www.stats.ox.ac.uk/~snijders/siena/
 *
 * File: IntersectionTieIterator.h
 *
 * Description: This module defines an iterator that only accesses the
 * intersection of two iterators.
 *****************************************************************************/

#ifndef INTERSECTIONTIEITERATOR_H_
#define INTERSECTIONTIEITERATOR_H_

#include "CombinedTieIterator.h"

namespace siena {

// ----------------------------------------------------------------------------
// Section: IntersectionTieIterator class
// ----------------------------------------------------------------------------

class IntersectionTieIterator: public CombinedTieIterator {
public:

	IntersectionTieIterator(const ITieIterator& iter1,
			const ITieIterator& iter2);

	virtual ~IntersectionTieIterator();

	virtual void next();

	virtual int actor() const;

	virtual bool valid() const;

	virtual IntersectionTieIterator * clone() const;

	void skip();

protected:

	IntersectionTieIterator(const IntersectionTieIterator& rhs);

private:

	/// Disable assignment operator.
	IntersectionTieIterator& operator=(const IntersectionTieIterator&);

};

} /* namespace siena */

#endif /* INTERSECTIONTIEITERATOR_H_ */
