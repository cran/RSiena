/******************************************************************************
 * SIENA: Simulation Investigation for Empirical Network Analysis
 *
 * Web: http://www.stats.ox.ac.uk/~snijders/siena/
 *
 * File: PermittedChangeFilter.h
 *
 * Description: This file contains the definition of the
 * PermittedChangeFilter class.
 *****************************************************************************/

#ifndef PERMITTEDCHANGEFILTER_H_
#define PERMITTEDCHANGEFILTER_H_

namespace siena
{

// ----------------------------------------------------------------------------
// Section: Forward declarations
// ----------------------------------------------------------------------------

class NetworkVariable;


// ----------------------------------------------------------------------------
// Section: Class definition
// ----------------------------------------------------------------------------

/**
 * Defines a filter that may forbid some tie changes for a given ego, when
 * the ego is given the opportunity to change a tie in a certain network.
 */
class PermittedChangeFilter
{
public:
	PermittedChangeFilter(const NetworkVariable * pVariable);

	/**
	 * Forbids tie changes between the given ego and some alters
	 * by setting the permitted flag to false for these alters.
	 */
	virtual void filterPermittedChanges(int ego, bool * permitted) = 0;

protected:
	inline const NetworkVariable * pVariable() const;

private:
	// The network variable, whose changes are controlled by this
	// filter.

	const NetworkVariable * lpVariable;
};


// ----------------------------------------------------------------------------
// Section: Inline methods
// ----------------------------------------------------------------------------

/**
 * Returns the network variable, whose changes are controlled by this filter.
 */
const NetworkVariable * PermittedChangeFilter::pVariable() const
{
	return this->lpVariable;
}

}

#endif /* PERMITTEDCHANGEFILTER_H_ */
