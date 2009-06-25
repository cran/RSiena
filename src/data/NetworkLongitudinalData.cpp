/******************************************************************************
 * SIENA: Simulation Investigation for Empirical Network Analysis
 * 
 * Web: http://www.stats.ox.ac.uk/~snijders/siena/
 * 
 * File: NetworkLongitudinalData.cpp
 * 
 * Description: This file contains the implementation of the
 * NetworkLongitudinalData class.
 *****************************************************************************/

#include <limits>
#include "NetworkLongitudinalData.h"
#include "data/Network.h"
#include "data/OneModeNetwork.h"
#include "data/ActorSet.h"
#include "model/variables/NetworkVariable.h"

namespace siena
{

// ----------------------------------------------------------------------------
// Section: Constructors, destructors.
// ----------------------------------------------------------------------------

/**
 * Creates a data object for storing the given number of observations of a
 * network. Initially the networks are empty at each observation.
 * @param[in] name the name of the corresponding network variable
 * @param[in] pSenders the set of actors acting as senders of ties
 * @param[in] pReceivers the set of actors acting as receivers of ties
 * @param[in] observationCount the number of observations to be stored
 */
NetworkLongitudinalData::NetworkLongitudinalData(
	std::string name,
	const ActorSet * pSenders,
	const ActorSet * pReceivers,
	int observationCount) : LongitudinalData(name, pSenders, observationCount)
{
	this->lpReceivers = pReceivers;
	this->lnetworks = new Network * [observationCount];
	this->lstructuralTieNetworks = new Network * [observationCount];
	this->lmissingTieNetworks = new Network * [observationCount];
	this->lmaxDegree = std::numeric_limits<int>::max();
	
	for (int i = 0; i < observationCount; i++)
	{
		if (pSenders == pReceivers)
		{
			this->lnetworks[i] = new OneModeNetwork(pSenders->n(), false);
			this->lstructuralTieNetworks[i] =
				new OneModeNetwork(pSenders->n(), false);
			this->lmissingTieNetworks[i] =
				new OneModeNetwork(pSenders->n(), false);
		}
		else
		{
			this->lnetworks[i] = new Network(pSenders->n(), pReceivers->n());
			this->lstructuralTieNetworks[i] =
				new Network(pSenders->n(), pReceivers->n());
			this->lmissingTieNetworks[i] =
				new Network(pSenders->n(), pReceivers->n());			
		}
	}
}


/**
 * Deallocates this data object including all its networks.
 */
NetworkLongitudinalData::~NetworkLongitudinalData()
{
	for (int i = 0; i < this->observationCount(); i++)
	{
		delete this->lnetworks[i];
		delete this->lstructuralTieNetworks[i];
		delete this->lmissingTieNetworks[i];
	}
	
	delete[] this->lnetworks;
	delete[] this->lstructuralTieNetworks;
	delete[] this->lmissingTieNetworks;
	
	this->lnetworks = 0;
	this->lstructuralTieNetworks = 0;
	this->lmissingTieNetworks = 0;
}


// ----------------------------------------------------------------------------
// Section: Accessors
// ----------------------------------------------------------------------------

/**
 * Returns the set of actors acting as tie senders.
 */
const ActorSet * NetworkLongitudinalData::pSenders() const
{
	return this->pActorSet();
}


/**
 * Returns the set of actors acting as tie receivers.
 */
const ActorSet * NetworkLongitudinalData::pReceivers() const
{
	return this->lpReceivers;
}


/**
 * Returns the network as of the given observation moment. This method
 * should be used for populating the data object, as it provides write
 * access to the network.
 */
Network * NetworkLongitudinalData::pNetwork(int observation)
{
	return this->lnetworks[observation];
}


/**
 * Returns the observed network as of the given observation moment.
 */
const Network * NetworkLongitudinalData::pNetwork(int observation) const
{
	return this->lnetworks[observation];
}


/**
 * Returns the network storing the structural tie indicators for the given
 * observation. This method should be used for populating the data object,
 * as it provides write access to the network.
 */
Network * NetworkLongitudinalData::pStructuralTieNetwork(int observation)
{
	return this->lstructuralTieNetworks[observation];
}


/**
 * Returns the network storing the structural tie indicators for the given
 * observation.
 */
const Network * NetworkLongitudinalData::pStructuralTieNetwork(int observation)
	const
{
	return this->lstructuralTieNetworks[observation];
}


/**
 * Returns the network storing the missing tie indicators for the given
 * observation. This method should be used for populating the data object,
 * as it provides write access to the network.
 */
Network * NetworkLongitudinalData::pMissingTieNetwork(int observation)
{
	return this->lmissingTieNetworks[observation];
}


/**
 * Returns the network storing the missing tie indicators for the given
 * observation.
 */
const Network * NetworkLongitudinalData::pMissingTieNetwork(int observation)
	const
{
	return this->lmissingTieNetworks[observation];
}


/**
 * Returns the observed value of a tie from <i>i</i> to <i>j</i> at the given
 * observation.
 */
int NetworkLongitudinalData::tieValue(int i, int j, int observation) const
{
	return this->lnetworks[observation]->tieValue(i, j);
}


/**
 * Returns if the tie value between the given actors is missing at the
 * given observation.
 */
bool NetworkLongitudinalData::missing(int i, int j, int observation) const
{
	return this->lmissingTieNetworks[observation]->tieValue(i, j);
}


/**
 * Returns the number of structurally determined tie variables from the given
 * actor at the given observation.
 */
int NetworkLongitudinalData::structuralTieCount(int actor, int observation)
	const
{
	return this->lstructuralTieNetworks[observation]->outDegree(actor);
}


/**
 * Stores the maximum permitted out-degree of an actor.
 */
void NetworkLongitudinalData::maxDegree(int degree)
{
	this->lmaxDegree = degree;
}


/**
 * Returns the maximum permitted out-degree of an actor.
 */
int NetworkLongitudinalData::maxDegree() const
{
	return this->lmaxDegree;
}


// ----------------------------------------------------------------------------
// Section: Overrides
// ----------------------------------------------------------------------------

/**
 * Returns a new network variable based on this longitudinal data. The
 * caller must take care of the deallocation of the returned variable.
 * @param[in] pSimulation the simulation using the new variable
 */
DependentVariable * NetworkLongitudinalData::createVariable(
	EpochSimulation * pSimulation)
{
	return new NetworkVariable(this, pSimulation);
}

}
