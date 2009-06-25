/******************************************************************************
 * SIENA: Simulation Investigation for Empirical Network Analysis
 *
 * Web: http://www.stats.ox.ac.uk/~snijders/siena/
 *
 * File: OneModeNetworkLongitudinalData.cpp
 *
 * Description: This file contains the implementation of the
 * OneModeNetworkLongitudinalData class.
 *****************************************************************************/

#include "OneModeNetworkLongitudinalData.h"

namespace siena
{

// ----------------------------------------------------------------------------
// Section: Construction, destruction
// ----------------------------------------------------------------------------

/**
 * Constructs a data object for storing observations of a one-mode network.
 * @param[in] name the name of the network
 * @param[in] pActors the set of actors underlying the network
 * @param[in] observationCount the number of observations to be stored
 */
OneModeNetworkLongitudinalData::OneModeNetworkLongitudinalData(
	std::string name,
	const ActorSet * pActors,
	int observationCount) :
		NetworkLongitudinalData(name, pActors, pActors, observationCount)
{
	this->lsymmetric = false;
	this->lbalanceMean = 0;
}


/**
 * Destroys this data object.
 */
OneModeNetworkLongitudinalData::~OneModeNetworkLongitudinalData()
{
}


// ----------------------------------------------------------------------------
// Section: Accessors
// ----------------------------------------------------------------------------

/**
 * Returns if the network is symmetric at all observations. It is just a
 * storage of the flag computed in R, and it is not computed dynamically.
 */
bool OneModeNetworkLongitudinalData::symmetric() const
{
	return this->lsymmetric;
}


/**
 * Stores if the network is supposed to be symmetric at all observations.
 */
void OneModeNetworkLongitudinalData::symmetric(bool flag)
{
	this->lsymmetric = flag;
}


/**
 * Returns the centering constant for the balance effect. This class just
 * stores the value computed in R.
 */
double OneModeNetworkLongitudinalData::balanceMean() const
{
	return this->lbalanceMean;
}


/**
 * Stores the centering constant for the balance effect.
 */
void OneModeNetworkLongitudinalData::balanceMean(double value)
{
	this->lbalanceMean = value;
}

}
