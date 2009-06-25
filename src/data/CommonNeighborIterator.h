#ifndef COMMONNEIGHBORITERATOR_H_
#define COMMONNEIGHBORITERATOR_H_

#include "IncidentTieIterator.h"

namespace siena
{

/**
 * This class defines an iterator over actors that are common to a pair
 * of incident tie iterators.
 */
class CommonNeighborIterator
{
public:
	CommonNeighborIterator(IncidentTieIterator iter1,
		IncidentTieIterator iter2);
	
	bool valid() const;
	int neighbor() const;
	void next();
	
private:
	void skipMismatches();
	
	IncidentTieIterator liter1;
	IncidentTieIterator liter2;
};

}

#endif /*COMMONNEIGHBORITERATOR_H_*/
