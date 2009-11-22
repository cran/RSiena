/******************************************************************************
 * SIENA: Simulation Investigation for Empirical Network Analysis
 *
 * Web: http://www.stats.ox.ac.uk/~snijders/siena/
 *
 * File: DyadicCovariateDependentNetworkEffect.cpp
 *
 * Description: This file contains the implementation of the
 * DyadicCovariateDependentNetworkEffect class.
 *****************************************************************************/

#include <stdexcept>

#include "DyadicCovariateDependentNetworkEffect.h"
#include "data/ConstantDyadicCovariate.h"
#include "data/ChangingDyadicCovariate.h"
#include "data/DyadicCovariateValueIterator.h"
#include "model/State.h"
#include "model/EffectInfo.h"
#include "model/EpochSimulation.h"

namespace siena
{

/**
 * Constructor.
 */
DyadicCovariateDependentNetworkEffect::DyadicCovariateDependentNetworkEffect(
	const EffectInfo * pEffectInfo) : NetworkEffect(pEffectInfo)
{
	this->lpConstantCovariate = 0;
	this->lpChangingCovariate = 0;
}


/**
 * Initializes this effect.
 * @param[in] pData the observed data
 * @param[in] pState the current state of the dependent variables
 * @param[in] period the period of interest
 * @param[in] pCache the cache object to be used to speed up calculations
 */
void DyadicCovariateDependentNetworkEffect::initialize(const Data * pData,
	State * pState,
	int period,
	Cache * pCache)
{
	NetworkEffect::initialize(pData, pState, period, pCache);
	string name = this->pEffectInfo()->interactionName1();

	this->lpConstantCovariate =	pData->pConstantDyadicCovariate(name);
	this->lpChangingCovariate =	pData->pChangingDyadicCovariate(name);

	if (!this->lpConstantCovariate && !this->lpChangingCovariate)
	{
		throw logic_error(
			"Dyadic covariate variable '" + name + "' expected.");
	}
}


/**
 * Returns the covariate value for the given pair of actors.
 */
double DyadicCovariateDependentNetworkEffect::value(int i, int j) const
{
	double value = 0;

	if (this->lpConstantCovariate)
	{
		value = this->lpConstantCovariate->value(i, j) -
			this->lpConstantCovariate->mean();
	}
	else
	{
		value = this->lpChangingCovariate->value(i, j, this->period()) -
			this->lpChangingCovariate->mean();
	}

	return value;
}


/**
 * Returns if the covariate value for the given pair of actors is missing.
 */
bool DyadicCovariateDependentNetworkEffect::missing(int i, int j) const
{
	bool missing = false;

	if (this->lpConstantCovariate)
	{
		missing = this->lpConstantCovariate->missing(i, j);
	}
	else
	{
		missing = this->lpChangingCovariate->missing(i, j, this->period());
	}

	return missing;
}


/**
 * Returns an iterator over non-zero non-missing values of the given row
 * of the covariate.
 */
DyadicCovariateValueIterator
	DyadicCovariateDependentNetworkEffect::rowValues(int i) const
{
	if (this->lpConstantCovariate)
	{
		return this->lpConstantCovariate->rowValues(i);
	}
	else
	{
		return this->lpChangingCovariate->rowValues(i, this->period());
	}
}


/**
 * Returns an iterator over non-zero non-missing values of the given column
 * of the covariate.
 */
DyadicCovariateValueIterator
	DyadicCovariateDependentNetworkEffect::columnValues(int j) const
{
	if (this->lpConstantCovariate)
	{
		return this->lpConstantCovariate->columnValues(j);
	}
	else
	{
		return this->lpChangingCovariate->columnValues(j, this->period());
	}
}

}