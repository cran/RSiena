/******************************************************************************
 * SIENA: Simulation Investigation for Empirical Network Analysis
 * 
 * Web: http://www.stats.ox.ac.uk/~snijders/siena/
 * 
 * File: Covariate.cpp
 * 
 * Description: This file contains the implementation of the
 * Covariate class.
 *****************************************************************************/

#include <cmath>

#include "Covariate.h"

namespace siena
{

// ----------------------------------------------------------------------------
// Section: Construction
// ----------------------------------------------------------------------------

/**
 * Creates a new covariate with the given name for the given set of actors.
 */
Covariate::Covariate(std::string name, const ActorSet * pActorSet) :
	NamedObject(name)
{
	this->lpActorSet = pActorSet;
	this->lrange = 0;
	this->lsimilarityMean = 0;
}


// ----------------------------------------------------------------------------
// Section: Accessors
// ----------------------------------------------------------------------------

/**
 * Returns the domain set of actors for this covariate.
 */
const ActorSet * Covariate::pActorSet() const
{
	return this->lpActorSet;
}


/**
 * Stores the range of values of this covariate, which is calculated in R.
 */
void Covariate::range(double range)
{
	this->lrange = range;
}


/**
 * Stores the similarity mean of this covariate, which is calculated in R.
 */
void Covariate::similarityMean(double similarityMean)
{
	this->lsimilarityMean = similarityMean;
}


// ----------------------------------------------------------------------------
// Section: Similarity
// ----------------------------------------------------------------------------

/**
 * Returns the centered similarity for the given values defined as
 * 1 - |a - b| / range - similarityMean.
 */
double Covariate::similarity(double a, double b) const
{
	return 1.0 - fabs(a - b) / this->lrange - this->lsimilarityMean;
}

}
