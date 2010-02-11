/******************************************************************************
 * SIENA: Simulation Investigation for Empirical Network Analysis
 *
 * Web: http://www.stats.ox.ac.uk/~snijders/siena/
 *
 * File: NetworkChange.h
 *
 * Description: This file contains the definition of the NetworkChange class.
 *****************************************************************************/

#ifndef NETWORKCHANGE_H_
#define NETWORKCHANGE_H_

#include "MiniStep.h"

namespace siena
{

/**
 * Defines a ministep changing a network variable.
 */
class NetworkChange: public MiniStep
{
public:
	NetworkChange(int ego,
		int alter,
		string variableName,
		int difference);
	virtual ~NetworkChange();

	inline int alter() const;

	virtual void makeChange(DependentVariable * pVariable);

private:
	// The alter whose incoming tie is changed
	int lalter;
};


// ----------------------------------------------------------------------------
// Section: Inline methods
// ----------------------------------------------------------------------------

/**
 * Returns the alter of this ministep.
 */
int NetworkChange::alter() const
{
	return this->lalter;
}

}

#endif /* NETWORKCHANGE_H_ */
