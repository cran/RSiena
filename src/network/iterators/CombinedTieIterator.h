/******************************************************************************
 * SIENA: Simulation Investigation for Empirical Network Analysis
 *
 * Web: http://www.stats.ox.ac.uk/~snijders/siena/
 *
 * File: CombinedTieIterator.h
 *
 * Description: This abstract iterator allows to combine two ITieIterators.
 *****************************************************************************/

#ifndef COMBINEDTIEITERATOR_H_
#define COMBINEDTIEITERATOR_H_

#include "ITieIterator.h"

namespace siena {

// ----------------------------------------------------------------------------
// Section: CombinedTieIterator abstract class
// ----------------------------------------------------------------------------

class CombinedTieIterator: public ITieIterator {
public:

	virtual ~CombinedTieIterator();

protected:

	CombinedTieIterator(const ITieIterator& iter1, const ITieIterator& iter2);

	CombinedTieIterator(const CombinedTieIterator& rhs);

	/**
	 * Tells whether both iterators point to the same actor.
	 * @return <code>True</code> if both iterators point to the same
	 * actor and <code>False</code> otherwise.
	 */
	inline bool isCommon() const {
		return lpIter1->actor() == lpIter2->actor();
	}

protected:

	/**
	 * The first iterator.
	 */
	ITieIterator* const lpIter1;

	/**
	 * The second iterator.
	 */
	ITieIterator* const lpIter2;

private:

	/// Disable assignment operator.
	CombinedTieIterator& operator=(const CombinedTieIterator&);
};

} /* namespace siena */
#endif /* COMBINEDTIEITERATOR_H_ */
