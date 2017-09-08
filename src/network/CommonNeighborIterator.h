#ifndef COMMONNEIGHBORITERATOR_H_
#define COMMONNEIGHBORITERATOR_H_

#include "IncidentTieIterator.h"
#include "iterators/ITieIterator.h"

namespace siena {

/**
 * This class defines an iterator over actors that are common to a pair
 * of incident tie iterators.
 */
class CommonNeighborIterator: ITieIterator {
public:
	CommonNeighborIterator(const IncidentTieIterator& iter1,
			const IncidentTieIterator& iter2);
	CommonNeighborIterator(const CommonNeighborIterator& rhs);
	bool valid() const;
	int actor() const;
	void next();
	CommonNeighborIterator* clone() const;

private:
	void skipMismatches();

	IncidentTieIterator liter1;
	IncidentTieIterator liter2;

	CommonNeighborIterator& operator=(const CommonNeighborIterator&);
};

}

#endif /*COMMONNEIGHBORITERATOR_H_*/
