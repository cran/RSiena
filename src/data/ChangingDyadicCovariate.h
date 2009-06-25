/******************************************************************************
 * SIENA: Simulation Investigation for Empirical Network Analysis
 * 
 * Web: http://www.stats.ox.ac.uk/~snijders/siena/
 * 
 * File: ChangingDyadicCovariate.h
 * 
 * Description: This file contains the definition of the
 * ChangingDyadicCovariate class.
 *****************************************************************************/

#ifndef CHANGINGDYADICCOVARIATE_H_
#define CHANGINGDYADICCOVARIATE_H_

#include <map>
#include <set>
#include "DyadicCovariate.h"

using namespace std;

namespace siena
{

/**
 * This class defines a dyadic covariate that changes over time.
 */
class ChangingDyadicCovariate : public DyadicCovariate
{
public:
	ChangingDyadicCovariate(std::string name,
		const ActorSet * pFirstActorSet,
		const ActorSet * pSecondActorSet,
		int observationCount);
	virtual ~ChangingDyadicCovariate();

	double value(int i, int j, int observation) const;
	void value(int i, int j, int observation, double value);
	bool missing(int i, int j, int observation) const;
	void missing(int i, int j, int observation, bool flag);

private:
	// An array of arrays of maps storing the non-zero values of the covariate.
	// The value at observation k for a pair (i,j) is stored in
	// lpValues[k][i][j].
	
	map<int, double> ** lpValues;
	
	// An array of arrays of actor sets. Actor j belongs to lpMissings[k][i]
	// if and only if the covariate value for the pair (i,j) is missing at
	// observation k.
	
	set<int> ** lpMissings;

	// The number of observations
	int lobservationCount;
};

}

#endif /*CHANGINGDYADICCOVARIATE_H_*/
