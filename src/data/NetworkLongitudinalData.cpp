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
#include "network/Network.h"
#include "network/OneModeNetwork.h"
#include "data/ActorSet.h"

namespace siena
{

// ----------------------------------------------------------------------------
// Section: Constructors, destructors.
// ----------------------------------------------------------------------------

/**
 * Creates a data object for storing the given number of observations of a
 * network. Initially the networks are empty at each observation.
 * @param[in] id the ID that is unique among all longitudinal data object
 * of the parent Data instance
 * @param[in] name the name of the corresponding network variable
 * @param[in] pSenders the set of actors acting as senders of ties
 * @param[in] pReceivers the set of actors acting as receivers of ties
 * @param[in] observationCount the number of observations to be stored
 */
NetworkLongitudinalData::NetworkLongitudinalData(int id,
	std::string name,
	const ActorSet * pSenders,
	const ActorSet * pReceivers,
	int observationCount) :
		LongitudinalData(id, name, pSenders, observationCount)
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
// Section: Preprocessing
// ----------------------------------------------------------------------------

/**
 * Calculates various statistical properties from the stored network data.
 */
void NetworkLongitudinalData::calculateProperties()
{
	// Calculate the average indegree and outdegree.

	this->laverageInDegree = 0;
	this->laverageOutDegree = 0;

	for (int observation = 0;
		observation < this->observationCount();
		observation++)
	{
		Network * pNetwork = this->lnetworks[observation];

		for (int i = 0; i < this->lpReceivers->n(); i++)
		{
			this->laverageInDegree += pNetwork->inDegree(i);
		}

		for (int i = 0; i < this->pActorSet()->n(); i++)
		{
			this->laverageOutDegree += pNetwork->outDegree(i);
		}
	}

	this->laverageInDegree /=
		this->lpReceivers->n() * this->observationCount();
	this->laverageOutDegree /=
		this->pActorSet()->n() * this->observationCount();
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
 * Returns the observed network as of the given observation moment.
 */
const Network * NetworkLongitudinalData::pNetwork(int observation) const
{
	return this->lnetworks[observation];
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
 * observation.
 */
const Network * NetworkLongitudinalData::pMissingTieNetwork(int observation)
	const
{
	return this->lmissingTieNetworks[observation];
}


/**
 * Returns the observed value of the tie from <i>i</i> to <i>j</i> at the given
 * observation.
 */
int NetworkLongitudinalData::tieValue(int i, int j, int observation) const
{
	return this->lnetworks[observation]->tieValue(i, j);
}


/**
 * Stores the observed value of the tie from <i>i</i> to <i>j</i> at the given
 * observation.
 */
void NetworkLongitudinalData::tieValue(int i,
	int j,
	int observation,
	int value)
{
	this->lnetworks[observation]->setTieValue(i, j, value);
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
 * Stores if the tie value between the given actors is missing at the
 * given observation.
 */
void NetworkLongitudinalData::missing(int i, int j, int observation, bool flag)
{
	if (flag)
	{
		this->lmissingTieNetworks[observation]->setTieValue(i, j, 1);
	}
	else
	{
		this->lmissingTieNetworks[observation]->setTieValue(i, j, 0);
	}
}


/**
 * Returns if the tie value between the given actors is structurally determined
 * at the given observation.
 */
bool NetworkLongitudinalData::structural(int i, int j, int observation) const
{
	return this->lstructuralTieNetworks[observation]->tieValue(i, j);
}


/**
 * Stores if the tie value between the given actors is structurally determined
 * at the given observation.
 */
void NetworkLongitudinalData::structural(int i,
	int j,
	int observation,
	bool flag)
{
	if (flag)
	{
		this->lstructuralTieNetworks[observation]->setTieValue(i, j, 1);
	}
	else
	{
		this->lstructuralTieNetworks[observation]->setTieValue(i, j, 0);
	}
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
/**
 * Stores the average in-degree over all receivers and observations.
 */
void NetworkLongitudinalData::averageInDegree(double val)
{
	this->laverageInDegree = val;
}


/**
 * Store the average out-degree over all senders and observations.
 */
void NetworkLongitudinalData::averageOutDegree(double val)
{
	this->laverageOutDegree = val;
}


/**
 * Returns the average in-degree over all receivers and observations.
 */
double NetworkLongitudinalData::averageInDegree() const
{
	return this->laverageInDegree;
}


/**
 * Returns the average out-degree over all senders and observations.
 */
double NetworkLongitudinalData::averageOutDegree() const
{
	return this->laverageOutDegree;
}

}
