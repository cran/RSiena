/******************************************************************************
 * SIENA: Simulation Investigation for Empirical Network Analysis
 *
 * Web: http://www.stats.ox.ac.uk/~snijders/siena/
 *
 * File: MixedNetworkAlterFunction.h
 *
 * Description: This file contains the definition of the
 * MixedNetworkAlterFunction class.
 *****************************************************************************/


#ifndef MIXEDNETWORKALTERFUNCTION_H_
#define MIXEDNETWORKALTERFUNCTION_H_

#include <string>
#include "AlterFunction.h"
#include "utils/NamedObject.h"

using namespace std;

namespace siena
{

class Network;
class TwoNetworkCache;


class MixedNetworkAlterFunction: public AlterFunction
{
public:
	MixedNetworkAlterFunction(string firstNetworkName,
		string secondNetworkName);
	virtual ~MixedNetworkAlterFunction();

	virtual void initialize(const Data * pData,
		State * pState,
		int period,
		Cache * pCache);

protected:
	inline const Network * pFirstNetwork() const;
	inline const Network * pSecondNetwork() const;
	inline TwoNetworkCache * pTwoNetworkCache() const;

private:
	string lname1;
	string lname2;
	const Network * lpFirstNetwork;
	const Network * lpSecondNetwork;
	TwoNetworkCache * lpTwoNetworkCache;
};


// ----------------------------------------------------------------------------
// Section: Inline methods
// ----------------------------------------------------------------------------

const Network * MixedNetworkAlterFunction::pFirstNetwork() const
{
	return this->lpFirstNetwork;
}

const Network * MixedNetworkAlterFunction::pSecondNetwork() const
{
	return this->lpSecondNetwork;
}

TwoNetworkCache * MixedNetworkAlterFunction::pTwoNetworkCache() const
{
	return this->lpTwoNetworkCache;
}

}

#endif /* MIXEDNETWORKALTERFUNCTION_H_ */
