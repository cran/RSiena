/******************************************************************************
 * SIENA: Simulation Investigation for Empirical Network Analysis
 *
 * Web: http://www.stats.ox.ac.uk/~snijders/siena/
 *
 * File: Cache.h
 *
 * Description: This file contains the definition of the Cache class.
 *****************************************************************************/

#ifndef CACHE_H_
#define CACHE_H_

#include <map>

using namespace std;

namespace siena
{

class TwoNetworkCache;
class NetworkCache;
class Network;


class Cache
{
public:
	Cache();
	virtual ~Cache();

	NetworkCache * pNetworkCache(const Network * pNetwork);
	TwoNetworkCache * pTwoNetworkCache(const Network * pFirstNetwork,
		const Network * pSecondNetwork);
	void initialize(int ego);

private:
	map<const Network *, NetworkCache *> lnetworkCaches;
	map<const Network *, map <const Network *, TwoNetworkCache *> >
		ltwoNetworkCaches;
	int lego;
};

}

#endif /* CACHE_H_ */
