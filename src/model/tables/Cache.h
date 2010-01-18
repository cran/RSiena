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

class NetworkCache;
class Network;


class Cache
{
public:
	Cache();
	virtual ~Cache();

	NetworkCache * pNetworkCache(const Network * pNetwork);
	void initialize(int ego);

private:
	map<const Network *, NetworkCache *> lnetworkCaches;
	int lego;
};

}

#endif /* CACHE_H_ */
