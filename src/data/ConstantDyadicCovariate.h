/******************************************************************************
 * SIENA: Simulation Investigation for Empirical Network Analysis
 * 
 * Web: http://www.stats.ox.ac.uk/~snijders/siena/
 * 
 * File: ConstantDyadicCovariate.h
 * 
 * Description: This file contains the definition of the
 * ConstantDyadicCovariate class.
 *****************************************************************************/

#ifndef CONSTANTDYADICCOVARIATE_H_
#define CONSTANTDYADICCOVARIATE_H_

#include <map>
#include <set>
#include "DyadicCovariate.h"

using namespace std;

namespace siena
{

/**
 * This class defines a constant dyadic covariate.
 */
class ConstantDyadicCovariate : public DyadicCovariate
{
public:
	ConstantDyadicCovariate(std::string name,
		const ActorSet * pFirstActorSet,
		const ActorSet * pSecondActorSet);
	virtual ~ConstantDyadicCovariate();
	
	double value(int i, int j) const;
	void value(int i, int j, double value);
	bool missing(int i, int j) const;
	void missing(int i, int j, bool flag);
	
private:
	// An array of maps storing the non-zero values of the covariate.
	// A value for a pair (i,j) is stored in lpValues[i][j].
	
	map<int, double> * lpValues;

	// An array of sets of actors. Actor j belongs to lpMissings[i]
	// if and only if the covariate value for the pair (i,j) is missing.
	
	set<int> * lpMissings;
};

}

#endif /*CONSTANTDYADICCOVARIATE_H_*/
