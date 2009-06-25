/******************************************************************************
 * SIENA: Simulation Investigation for Empirical Network Analysis
 *
 * Web: http://www.stats.ox.ac.uk/~snijders/siena/
 *
 * File: OneModeNetworkLongitudinalData.h
 *
 * Description: This file contains the definition of the
 * OneModeNetworkLongitudinalData class.
 *****************************************************************************/

#ifndef ONEMODENETWORKLONGITUDINALDATA_H_
#define ONEMODENETWORKLONGITUDINALDATA_H_

#include "NetworkLongitudinalData.h"

namespace siena
{

/**
 * This class stores one or more observations of a one-mode network variable.
 */
class OneModeNetworkLongitudinalData : public NetworkLongitudinalData
{
public:
	OneModeNetworkLongitudinalData(std::string name,
		const ActorSet * pActors,
		int observationCount);
	virtual ~OneModeNetworkLongitudinalData();

	bool symmetric() const;
	void symmetric(bool flag);
	double balanceMean() const;
	void balanceMean(double value);

private:
	// Indicates is the networks are supposed to be symmetric at all
	// observations.

	bool lsymmetric;

	// The centering constant for the balance effect.
	double lbalanceMean;
};

}

#endif /*ONEMODENETWORKLONGITUDINALDATA_H_*/
