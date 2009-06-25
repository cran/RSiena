/******************************************************************************
 * SIENA: Simulation Investigation for Empirical Network Analysis
 *
 * Web: http://www.stats.ox.ac.uk/~snijders/siena/
 *
 * File: BehaviorLongitudinalData.cpp
 *
 * Description: This module implements the class BehaviorLongitudinalData.
 *****************************************************************************/

#include <limits>
#include <stdexcept>
#include <cmath>

#include "BehaviorLongitudinalData.h"
#include "utils/Utils.h"
#include "data/ActorSet.h"
#include "model/variables/BehaviorVariable.h"

namespace siena
{

/**
 * Constructs a data object for storing the observed values of a behavioral
 * variable for the given set of actors at the given number of observations.
 * Initially, all values are set to 0.
 * @param[in] name the name of the corresponding behavior variable
 */
BehaviorLongitudinalData::BehaviorLongitudinalData(std::string name,
	const ActorSet * pActorSet,
	int observationCount) :
		LongitudinalData(name, pActorSet, observationCount)
{
	this->lvalues = new int * [observationCount];
	this->lmissing = new bool * [observationCount];

	for (int i = 0; i < observationCount; i++)
	{
		this->lvalues[i] = new int[pActorSet->n()];
		this->lmissing[i] = new bool[pActorSet->n()];

		for (int actor = 0; actor < pActorSet->n(); actor++)
		{
			this->lvalues[i][actor] = 0;
			this->lmissing[i][actor] = false;
		}
	}
}


/**
 * Dealloacates this data object.
 */
BehaviorLongitudinalData::~BehaviorLongitudinalData()
{
	for (int i = 0; i < this->observationCount(); i++)
	{
		delete[] this->lvalues[i];
		delete[] this->lmissing[i];
	}

	delete[] this->lvalues;
	delete[] this->lmissing;
	this->lvalues = 0;
	this->lmissing = 0;
}


/**
 * Returns the observed value of the behavioral variable for the given actor
 * at the specified observation.
 */
int BehaviorLongitudinalData::value(int observation, int actor) const
{
	return this->lvalues[observation][actor];
}


/**
 * Stores the observed value of the behavioral variable for the given actor
 * at the specified observation.
 */
void BehaviorLongitudinalData::value(int observation,
	int actor,
	int value)
{
	this->lvalues[observation][actor] = value;
}


/**
 * Returns the whole array of observed values for the given observation.
 */
const int * BehaviorLongitudinalData::values(int observation) const
{
	return this->lvalues[observation];
}


/**
 * Returns if the value of the behavioral variable is missing for the given
 * actor at the specified observation.
 */
bool BehaviorLongitudinalData::missing(int observation, int actor) const
{
	return this->lmissing[observation][actor];
}


/**
 * Stores if the value of the behavioral variable is missing for the given
 * actor at the specified observation.
 */
void BehaviorLongitudinalData::missing(int observation,
	int actor,
	bool missing)
{
	this->lmissing[observation][actor] = missing;
}


/**
 * Returns a new behavior variable based on this longitudinal data. The
 * caller must take care of the deallocation of the returned variable.
 * @param[in] pSimulation the simulation using the new variable
 */
DependentVariable * BehaviorLongitudinalData::createVariable(
	EpochSimulation * pSimulation)
{
	return new BehaviorVariable(this, pSimulation);
}


// ----------------------------------------------------------------------------
// Section: Various statistics on the observed values
// ----------------------------------------------------------------------------

/**
 * Returns the smallest observed value.
 */
int BehaviorLongitudinalData::min() const
{
	return this->lmin;
}


/**
 * Returns the largest observed value.
 */
int BehaviorLongitudinalData::max() const
{
	return this->lmax;
}


/**
 * Returns the overall mean value over all observations.
 */
double BehaviorLongitudinalData::overallMean() const
{
	return this->loverallMean;
}


/**
 * Returns the centered similarity for the given values defined as
 * 1 - |a - b| / range - similarityMean.
 */
double BehaviorLongitudinalData::similarity(double a, double b) const
{
	return 1.0 - fabs(a - b) / this->lrange - this->lsimilarityMean;
}


/**
 * Returns the similarity mean value over all observations.
 */
double BehaviorLongitudinalData::similarityMean() const
{
	return this->lsimilarityMean;
}


/**
 * Stores the similarity mean value over all observations.
 */
void BehaviorLongitudinalData::similarityMean(double similarityMean)
{
	this->lsimilarityMean = similarityMean;
}


/**
 * Returns the range of observed values.
 */
int BehaviorLongitudinalData::range() const
{
	return this->lrange;
}


/**
 * Calculates various statistics from the observed data.
 */
void BehaviorLongitudinalData::calculateStatistics()
{
	this->lmin = numeric_limits<int>::max();
	this->lmax = numeric_limits<int>::min();
	this->loverallMean = 0;

	for (int i = 0; i < this->observationCount(); i++)
	{
		int nonMissingValueCount = 0;
		int sum = 0;

		for (int actor = 0; actor < this->n(); actor++)
		{
			if (!this->lmissing[i][actor])
			{
				this->lmin = std::min(this->lmin, this->lvalues[i][actor]);
				this->lmax = std::max(this->lmax, this->lvalues[i][actor]);
				sum += this->lvalues[i][actor];
				nonMissingValueCount++;
			}
		}

		if (nonMissingValueCount == 0)
		{
			throw logic_error(
				"No valid data for dependent actor variable '" +
				this->name() +
				"', period " + toString(i));
		}

		this->loverallMean += ((double) sum) / nonMissingValueCount;
	}

	this->loverallMean /= this->observationCount();
	this->lrange = this->lmax - this->lmin;

	if (this->lrange == 0)
	{
		throw logic_error(
			"All observed values are equal for the behavior variable " +
			this->name());
	}
}

}
