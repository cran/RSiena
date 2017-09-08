/******************************************************************************
 * SIENA: Simulation Investigation for Empirical Network Analysis
 *
 * Web: http://www.stats.ox.ac.uk/~snijders/siena/
 *
 * File: SymDiffTieIterator.h
 *
 * Description: This module defines an iterator that only accesses the
 * symmetric difference of two iterators.
 *****************************************************************************/

#ifndef SYMDIFFTIEITERATOR_H_
#define SYMDIFFTIEITERATOR_H_

#include "UnionTieIterator.h"

namespace siena {

// ----------------------------------------------------------------------------
// Section: SymDiffTieIterator class
// ----------------------------------------------------------------------------

class SymDiffTieIterator: public UnionTieIterator {

public:

	SymDiffTieIterator(const ITieIterator& iter1, const ITieIterator& iter2);

	SymDiffTieIterator(const ITieIterator& iter1, const ITieIterator& iter2,
			int idIter1, int idIter2);

	virtual ~SymDiffTieIterator();

	virtual SymDiffTieIterator* clone() const;

	virtual void next();

protected:

	SymDiffTieIterator(const SymDiffTieIterator& rhs);

private:

	void init();

};

} /* namespace siena */

#endif /* SYMDIFFTIEITERATOR_H_ */
