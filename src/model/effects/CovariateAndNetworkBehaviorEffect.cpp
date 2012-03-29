/******************************************************************************
 * SIENA: Simulation Investigation for Empirical Network Analysis
 *
 * Web: http://www.stats.ox.ac.uk/~snijders/siena/
 *
 * File: CovariateAndNetworkBehaviorEffect.cpp
 *
 * Description: This file contains the implementation of the
 * CovariateAndNetworkBehaviorEffect class.
 *****************************************************************************/

#include <stdexcept>
//#include "R_ext/Print.h"

#include "CovariateAndNetworkBehaviorEffect.h"
#include "network/Network.h"
#include "network/IncidentTieIterator.h"
#include "model/State.h"
#include "model/EffectInfo.h"

namespace siena
{

/**
 * Constructor.
 */
CovariateAndNetworkBehaviorEffect::CovariateAndNetworkBehaviorEffect(
	const EffectInfo * pEffectInfo) :
	CovariateDependentBehaviorEffect(pEffectInfo)
{
	// set up the extras if any
	this->laverageAlterValues = 0;
	this->laverageAlterMissing = 0;
}

/**
 * Deallocates this effect object;
 */
CovariateAndNetworkBehaviorEffect::~CovariateAndNetworkBehaviorEffect()
{
	delete [] this->laverageAlterValues;
	delete [] this->laverageAlterMissing;
}

/**
 * Initializes this effect.
 * @param[in] pData the observed data
 * @param[in] pState the current state of the dependent variables
 * @param[in] period the period of interest
 * @param[in] pCache the cache object to be used to speed up calculations
 */
void CovariateAndNetworkBehaviorEffect::initialize(const Data * pData,
	State * pState,
	int period,
	Cache * pCache)
{
	CovariateDependentBehaviorEffect::initialize(pData, pState, period, pCache);

	string networkName = this->pEffectInfo()->interactionName2();

	this->lpNetwork = pState->pNetwork(networkName);

	if (!this->lpNetwork)
	{
		throw logic_error("Network '" + networkName + "' expected.");
	}
	// initialize extras if any
	if (this->laverageAlterValues)
	{
		delete [] this->laverageAlterValues;
	}
	if (this->laverageAlterMissing)
	{
		delete [] this->laverageAlterMissing;
	}
	this->laverageAlterValues = new double[this->lpNetwork->n()];
	this->laverageAlterMissing = new bool[this->lpNetwork->n()];
}

/**
 * Returns if the dummy covariate value for the given actor is based on
 * all missing values.
 */
bool CovariateAndNetworkBehaviorEffect::missingDummy(int i) const
{
	return this->laverageAlterMissing[i];
}

/**
 * Returns the average alter covariate value for the given actor.
 */
double CovariateAndNetworkBehaviorEffect::averageAlterValue(int i) const
{
	return this->laverageAlterValues[i];
}


/**
 * Does the necessary preprocessing work for calculating the tie flip
 * contributions for a specific ego. This method must be invoked before
 * calling NetworkEffect::calculateTieFlipContribution(...).
 */
void CovariateAndNetworkBehaviorEffect::preprocessEgo(int ego)
{
	//CovariateDependentBehaviorEffect::preprocessEgo(ego);

	// set up the covariate based on current values of the network
	const Network * pNetwork = this->pNetwork();


	for (int i = 0; i < pNetwork->n(); i++)
	{
		this->laverageAlterMissing[i] = false;
		int numberNonMissing = 0;
		this->laverageAlterValues[i] = 0;
		if (pNetwork->outDegree(i) > 0)
		{
			for (IncidentTieIterator iter = pNetwork->outTies(i);
				 iter.valid();
				 iter.next())
			{
				int j = iter.actor();
				this->laverageAlterValues[i] += this->covariateValue(j);
				if (!this->missingCovariate(j, this->period()))
				{
					numberNonMissing++;
				}
// 				Rprintf("%d %f %d %d %d %d\n",
// 					j,
// 					this->covariateValue(j),
// 					this->period(),
// 					this->missingCovariate(j, this->period()),
// 					numberNonMissing, i);
			}
			this->laverageAlterValues[i] /= pNetwork->outDegree(i);
			if (numberNonMissing == 0)
			{
				this->laverageAlterMissing[i] = true;
			}
		}
		else
		{
			this->laverageAlterValues[i] = 0;
		}
//		Rprintf("%d %f\n", i,this->laverageAlterValues[i]);
	}
}


}
