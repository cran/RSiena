/******************************************************************************
 * SIENA: Simulation Investigation for Empirical Network Analysis
 *
 * Web: http://www.stats.ox.ac.uk/~snijders/siena/
 *
 * File: DataUtils.h
 *
 * Description: This module contains some utilities specific to the
 * 'data' library.
 *****************************************************************************/

#ifndef DATAUTILS_H_
#define DATAUTILS_H_

namespace siena
{

// ----------------------------------------------------------------------------
// Section: Enums
// ----------------------------------------------------------------------------

/**
 * This enumeration defines two possible directions of traversing a tie.
 */
enum Direction {FORWARD, BACKWARD, RECIPROCAL};


// ----------------------------------------------------------------------------
// Section: Forward declarations
// ----------------------------------------------------------------------------

class IncidentTieIterator;
class Network;


// ----------------------------------------------------------------------------
// Section: Utility functions
// ----------------------------------------------------------------------------

int commonActorCount(IncidentTieIterator iter1, IncidentTieIterator iter2);
Network * symmetricDifference(const Network * pNetwork1,
	const Network * pNetwork2);
}

#endif /*DATAUTILS_H_*/
