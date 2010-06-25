/******************************************************************************
 * SIENA: Simulation Investigation for Empirical Network Analysis
 *
 * Web: http://www.stats.ox.ac.uk/~snijders/siena/
 *
 * File: CovariateDistance2DependentNetworkEffect.cpp
 *
 * Description: This file contains the implementation of the
 * CovariateDistance2DependentNetworkEffect class.
 *****************************************************************************/

#include <stdexcept>
#include "R_ext/Print.h"
#include "CovariateDistance2NetworkEffect.h"
#include "data/ConstantCovariate.h"
#include "data/ChangingCovariate.h"
#include "data/BehaviorLongitudinalData.h"
#include "network/IncidentTieIterator.h"
#include "network/Network.h"
#include "model/State.h"
#include "model/EffectInfo.h"
#include "model/EpochSimulation.h"
#include "model/variables/BehaviorVariable.h"

namespace siena
{

/**
 * Constructor.
 */
CovariateDistance2NetworkEffect::CovariateDistance2NetworkEffect(
	const EffectInfo * pEffectInfo) :
	CovariateDependentNetworkEffect(pEffectInfo)
{
	// set up the extras if any
	this->laverageAlterValues = 0;
	this->laverageAlterMissing = 0;
}

/**
 * Deallocates this effect object;
 */
CovariateDistance2NetworkEffect::~CovariateDistance2NetworkEffect()
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
void CovariateDistance2NetworkEffect::initialize(const Data * pData,
	State * pState,
	int period,
	Cache * pCache)
{
	CovariateDependentNetworkEffect::initialize(pData, pState, period, pCache);
	// initialize extras if any
	this->laverageAlterValues = new double[this->pNetwork()->n()];
	this->laverageAlterMissing = new bool[this->pNetwork()->n()];
}

/**
 * Returns if the dummy covariate value for the given actor is based on
 * all missing values.
 */
bool CovariateDistance2NetworkEffect::missingDummy(int i) const
{
	return this->laverageAlterMissing[i];
}

/**
 * Returns the average alter covariate value for the given actor.
 */
double CovariateDistance2NetworkEffect::averageAlterValue(int i) const
{
	return this->laverageAlterValues[i];
}


/**
 * Does the necessary preprocessing work for calculating the tie flip
 * contributions for a specific ego. This method must be invoked before
 * calling NetworkEffect::calculateTieFlipContribution(...).
 */
void CovariateDistance2NetworkEffect::preprocessEgo(int ego)
{
	NetworkEffect::preprocessEgo(ego);

	// set up the covariate based on current values of the network

	const Network * pNetwork = this->pNetwork();


	for (int i = 0; i < pNetwork->n(); i++)
	{
		int numberNonMissing = 0;
		this->laverageAlterMissing[i] = false;
		this->laverageAlterValues[i] = 0;
		if (pNetwork->outDegree(i) > 0)
		{
			for (IncidentTieIterator iter = pNetwork->outTies(i);
				 iter.valid();
				 iter.next())
			{
				int j = iter.actor();
				this->laverageAlterValues[i] += this->value(j);
				if (!this->missing(j))
				{
					numberNonMissing++;
				}
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
	}
}
/**
 * Returns the centered similarity of the given actors wrt to the given network.
 */
double CovariateDistance2NetworkEffect::similarityNetwork(int i, int j,
	std::string networkName) const
{
	double similarity = 0;

	if (this->pConstantCovariate())
	{
		similarity =
			this->pConstantCovariate()->similarityNetwork(
				this->averageAlterValue(i),
				this->averageAlterValue(j), networkName);
	}
	else if (this->pChangingCovariate())
	{
		similarity =
			this->pChangingCovariate()->similarityNetwork(
				this->averageAlterValue(i),
				this->averageAlterValue(j), networkName);
	}
	else
	{
		similarity =
			this->pBehaviorData()->similarityNetwork(
				this->averageAlterValue(i),
				this->averageAlterValue(j), networkName);
	}

	return similarity;
}
}
