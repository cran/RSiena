/******************************************************************************
 * SIENA: Simulation Investigation for Empirical Network Analysis
 *
 * Web: http://www.stats.ox.ac.uk/~snijders/siena/
 *
 * File: NetworkAlterFunction.h
 *
 * Description: This file contains the definition of the
 * NetworkAlterFunction class.
 *****************************************************************************/


#ifndef NETWORKALTERFUNCTION_H_
#define NETWORKALTERFUNCTION_H_

#include <string>
#include "AlterFunction.h"
#include "utils/NamedObject.h"

namespace siena
{

class Network;
class NetworkCache;

class NetworkAlterFunction: public AlterFunction, public NamedObject
{
public:
	NetworkAlterFunction(std::string networkName);
	virtual ~NetworkAlterFunction();

	virtual void initialize(const Data * pData,
		State * pState, int period, Cache * pCache);

protected:
	bool inTieExists(int alter) const;
	bool outTieExists(int alter) const;
	inline const Network * pNetwork() const;
	inline NetworkCache * pNetworkCache() const;

private:
	const Network * lpNetwork;
	NetworkCache * lpNetworkCache;
};


// ----------------------------------------------------------------------------
// Section: Inline methods
// ----------------------------------------------------------------------------

const Network * NetworkAlterFunction::pNetwork() const
{
	return this->lpNetwork;
}

NetworkCache * NetworkAlterFunction::pNetworkCache() const
{
	return this->lpNetworkCache;
}

}

#endif /* NETWORKALTERFUNCTION_H_ */
