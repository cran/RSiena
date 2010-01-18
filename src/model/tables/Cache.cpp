/******************************************************************************
 * SIENA: Simulation Investigation for Empirical Network Analysis
 *
 * Web: http://www.stats.ox.ac.uk/~snijders/siena/
 *
 * File: Cache.cpp
 *
 * Description: This file contains the implementation of the class Cache.
 *****************************************************************************/

#include "Cache.h"
#include "utils/Utils.h"
#include "network/Network.h"
#include "model/tables/NetworkCache.h"


namespace siena
{

Cache::Cache()
{
	this->lego = -1;
}


Cache::~Cache()
{
	clearMap(this->lnetworkCaches, false, true);
}


NetworkCache * Cache::pNetworkCache(const Network * pNetwork)
{
	NetworkCache * pNetworkCache = 0;
	map<const Network *, NetworkCache *>::iterator iter =
		this->lnetworkCaches.find(pNetwork);

	if (iter != this->lnetworkCaches.end())
	{
		pNetworkCache = iter->second;
	}
	else
	{
		pNetworkCache = new NetworkCache(pNetwork);
		pNetworkCache->initialize(this->lego);
		this->lnetworkCaches[pNetwork] = pNetworkCache;
	}

	return pNetworkCache;
}


void Cache::initialize(int ego)
{
	this->lego = ego;

	for (map<const Network *, NetworkCache *>::iterator iter =
			this->lnetworkCaches.begin();
		iter != this->lnetworkCaches.end();
		iter++)
	{
		NetworkCache * pNetworkCache = iter->second;
		pNetworkCache->initialize(ego);
	}
}

}
