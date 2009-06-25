/******************************************************************************
 * SIENA: Simulation Investigation for Empirical Network Analysis
 * 
 * Web: http://www.stats.ox.ac.uk/~snijders/siena/
 * 
 * File: Covariate.h
 * 
 * Description: This file contains the definition of the
 * Covariate class.
 *****************************************************************************/

#ifndef COVARIATE_H_
#define COVARIATE_H_

#include "utils/NamedObject.h"

namespace siena
{

// ----------------------------------------------------------------------------
// Section: Forward declarations
// ----------------------------------------------------------------------------

class ActorSet;


// ----------------------------------------------------------------------------
// Section: Covariate class
// ----------------------------------------------------------------------------

/**
 * This is the base class for constant and changing covariates.
 */
class Covariate : public NamedObject
{
public:
	Covariate(std::string name, const ActorSet * pActorSet);
	
	const ActorSet * pActorSet() const;
	
	inline double range() const;
	void range(double range);
	inline double similarityMean() const;
	void similarityMean(double similarityMean);
	
	double similarity(double a, double b) const;
	
private:
	// The underlying set of actors
	const ActorSet * lpActorSet;
	
	// The overall range of values (to be passed from outside)
	double lrange;
	
	// The similarity mean (to be passed from outside)
	double lsimilarityMean;
};


// ----------------------------------------------------------------------------
// Section: Inline methods
// ----------------------------------------------------------------------------

/**
 * Returns the range of values of this covariate.
 */
double Covariate::range() const
{
	return this->lrange;
}


/**
 * Returns the similarity mean of this covariate.
 */
double Covariate::similarityMean() const
{
	return this->lsimilarityMean;
}

}

#endif /*COVARIATE_H_*/
