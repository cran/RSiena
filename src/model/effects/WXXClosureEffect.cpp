/******************************************************************************
 * SIENA: Simulation Investigation for Empirical Network Analysis
 *
 * Web: http://www.stats.ox.ac.uk/~snijders/siena/
 *
 * File: WXXClosureEffect.cpp
 *
 * Description: This file contains the implementation of the
 * WXXClosureEffect class.
 *****************************************************************************/

#include "WXXClosureEffect.h"
#include "data/Network.h"
#include "data/IncidentTieIterator.h"
#include "data/NetworkLongitudinalData.h"
#include "data/DyadicCovariateValueIterator.h"
#include "model/variables/NetworkVariable.h"

namespace siena
{

/**
 * Constructor.
 */
WXXClosureEffect::WXXClosureEffect(const EffectInfo * pEffectInfo) :
	DyadicCovariateDependentNetworkEffect(pEffectInfo)
{
	this->lsums = 0;
}


/**
 * Destructor.
 */
WXXClosureEffect::~WXXClosureEffect()
{
	delete[] this->lsums;
	this->lsums = 0;
}


/**
 * Initializes this effect for the use with the given epoch simulation.
 */
void WXXClosureEffect::initialize(EpochSimulation * pSimulation)
{
	DyadicCovariateDependentNetworkEffect::initialize(pSimulation);

	delete[] this->lsums;
	this->lsums = new double[this->pVariable()->n()];
}


/**
 * Initializes this effect for calculating the corresponding statistics.
 * @param[in] pData the observed data
 * @param[in] pState the current state of the dependent variables
 * @param[in] period the period of interest
 */
void WXXClosureEffect::initialize(const Data * pData,
	State * pState,
	int period)
{
	DyadicCovariateDependentNetworkEffect::initialize(pData, pState, period);
}


/**
 * Does the necessary preprocessing work for calculating the tie flip
 * contributions for a specific ego. This method must be invoked before
 * calling NetworkEffect::calculateTieFlipContribution(...).
 */
void WXXClosureEffect::preprocessEgo()
{
	this->calculateSums(this->pVariable()->ego(),
		this->pVariable()->pNetwork(),
		this->lsums);
}


/**
 * For each j and the given i, this method calculates the sum
 * sum_h w_{ih} x_{hj}.
 */
void WXXClosureEffect::calculateSums(int i, Network * pNetwork, double * sums)
	const
{
	int n = pNetwork->n();

	// Initialize

	for (int j = 0; j < n; j++)
	{
		sums[j] = 0;
	}

	// Iterate over all h with non-zero non-missing w_{ih}

	for (DyadicCovariateValueIterator iterH = this->rowValues(i);
		iterH.valid();
		iterH.next())
	{
		int h = iterH.actor();

		// Iterate over all j with a tie from h

		for (IncidentTieIterator iterJ = pNetwork->outTies(h);
			iterJ.valid();
			iterJ.next())
		{
			int j = iterJ.actor();

			// Add the term w_{ih} x_{hj} (= w_{ih})
			sums[j] += iterH.value();
		}
	}
}


/**
 * Calculates the contribution of a tie flip to the given actor.
 */
double WXXClosureEffect::calculateTieFlipContribution(int alter) const
{
	double change = this->lsums[alter];

	if (this->pVariable()->outTieExists(alter))
	{
		change = -change;
	}

	return change;
}


/**
 * Detailed comment in the base class.
 */
double WXXClosureEffect::statistic(Network * pNetwork,
	Network * pSummationTieNetwork) const
{
	double statistic = 0;
	int n = pNetwork->n();
	double * sums = new double[n];

	for (int i = 0; i < n; i++)
	{
		this->calculateSums(i, pNetwork, sums);

		for (IncidentTieIterator iter = pSummationTieNetwork->outTies(i);
			iter.valid();
			iter.next())
		{
			statistic += sums[iter.actor()];
		}
	}

	delete[] sums;
	return statistic;
}

}
