/******************************************************************************
 * SIENA: Simulation Investigation for Empirical Network Analysis
 *
 * Web: http://www.stats.ox.ac.uk/~snijders/siena/
 *
 * File: BehaviorLongitudinalData.h
 *
 * Description: This module defines the class BehaviorLongitudinalData.
 *****************************************************************************/

#ifndef BEHAVIORLONGITUDINALDATA_H_
#define BEHAVIORLONGITUDINALDATA_H_

#include "data/LongitudinalData.h"

namespace siena
{

/**
 * This class stores the observed values of a behavior variable for one or more
 * observation moments.
 */
class BehaviorLongitudinalData : public LongitudinalData
{
public:
	BehaviorLongitudinalData(std::string name,
		const ActorSet * pActorSet,
		int observationCount);
	virtual ~BehaviorLongitudinalData();

	int value(int observation, int actor) const;
	void value(int observation, int actor, int value);
	bool missing(int observation, int actor) const;
	void missing(int observation, int actor, bool missing);
	const int * values(int observation) const;

	int min() const;
	int max() const;
	double overallMean() const;
	int range() const;
	double similarity(double a, double b) const;
	double similarityMean() const;
	void similarityMean(double similarityMean);
	void calculateProperties();

private:
	// An array of values per each observation
	int ** lvalues;

	// Missingness indicators
	bool ** lmissing;

	// The smallest non-missing value
	int lmin;

	// The largest non-missing value
	int lmax;

	// The overall mean value over all observations (zamean in Pascal)
	double loverallMean;

	// The range of observed values
	int lrange;

	// The similarity mean
	double lsimilarityMean;
};

}

#endif /*BEHAVIORLONGITUDINALDATA_H_*/
