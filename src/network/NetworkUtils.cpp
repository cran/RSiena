/******************************************************************************
 * SIENA: Simulation Investigation for Empirical Network Analysis
 *
 * Web: http://www.stats.ox.ac.uk/~snijders/siena/
 *
 * File: DataUtils.cpp
 *
 * Description: This module contains some utilities specific to the
 * 'data' library.
 *****************************************************************************/

#include "NetworkUtils.h"
#include "IncidentTieIterator.h"
#include "network/CommonNeighborIterator.h"
#include "network/Network.h"

namespace siena
{

/**
 * Returns the number of actors iterated over by both of the given iterators.
 */
int commonActorCount(IncidentTieIterator iter1, IncidentTieIterator iter2)
{
	// Note that ties incident to an actor are sorted in an increasing order
	// of its neighbors.

	int count = 0;

	CommonNeighborIterator iter(iter1, iter2);

	while (iter.valid())
	{
		count++;
		iter.next();
	}

	return count;
}


/**
 * Creates a new network representing the symmetric difference of the
 * two given networks.
 */
Network * symmetricDifference(const Network * pNetwork1,
	const Network * pNetwork2)
{
	Network * pDifference = new Network(pNetwork1->n(), pNetwork1->m());

	for (int i = 0; i < pNetwork1->n(); i++)
	{
		IncidentTieIterator iter1 = pNetwork1->outTies(i);
		IncidentTieIterator iter2 = pNetwork2->outTies(i);

		while (iter1.valid() && iter2.valid())
		{
			if (iter1.actor() < iter2.actor())
			{
				pDifference->setTieValue(i, iter1.actor(), 1);
				iter1.next();
			}
			else if (iter1.actor() > iter2.actor())
			{
				pDifference->setTieValue(i, iter2.actor(), 1);
				iter2.next();
			}
			else
			{
				iter1.next();
				iter2.next();
			}
		}

		while (iter1.valid())
		{
			pDifference->setTieValue(i, iter1.actor(), 1);
			iter1.next();
		}

		while (iter2.valid())
		{
			pDifference->setTieValue(i, iter2.actor(), 1);
			iter2.next();
		}
	}

	return pDifference;
}

}
