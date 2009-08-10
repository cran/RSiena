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
 * Initializes this effect for calculating contributions of changes to
 * evaluation or endowment functions.
 */
void DyadicCovariateDependentNetworkEffect::initialize(
	EpochSimulation * pSimulation)
{
	NetworkEffect::initialize(pSimulation);

	this->lpConstantCovariate =
		pSimulation->pData()->pConstantDyadicCovariate(
			this->pEffectInfo()->interactionName1());
	this->lpChangingCovariate =
		pSimulation->pData()->pChangingDyadicCovariate(
			this->pEffectInfo()->interactionName1());
//	this->lperiod = pSimulation->period();
//	Rprintf("init dyadic %x %d\n", this, this->lperiod);

	if (!this->lpConstantCovariate && !this->lpChangingCovariate)
	{
		throw logic_error("Dyadic covariate variable '" +
			this->pEffectInfo()->interactionName1() +
			"' expected.");
	}
}


/**
 * Initializes this effect for calculating the corresponding statistics.
 * @param[in] pData the observed data
 * @param[in] pState the current state of the dependent variables
 * @param[in] period the period of interest
 */
void DyadicCovariateDependentNetworkEffect::initialize(const Data * pData,
	State * pState,
	int period)
{
	NetworkEffect::initialize(pData, pState, period);

	this->lpConstantCovariate =
		pData->pConstantDyadicCovariate(
			this->pEffectInfo()->interactionName1());
	this->lpChangingCovariate =
		pData->pChangingDyadicCovariate(
			this->pEffectInfo()->interactionName1());

	if (!this->lpConstantCovariate && !this->lpChangingCovariate)
	{
		throw logic_error("Dyadic covariate variable '" +
			this->pEffectInfo()->interactionName1() +
			"' expected.");
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
