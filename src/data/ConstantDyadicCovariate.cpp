/******************************************************************************
 * SIENA: Simulation Investigation for Empirical Network Analysis
 * 
 * Web: http://www.stats.ox.ac.uk/~snijders/siena/
 * 
 * File: ConstantDyadicCovariate.cpp
 * 
 * Description: This file contains the implementation of the
 * ConstantDyadicCovariate class.
 *****************************************************************************/

#include "ConstantDyadicCovariate.h"
#include "data/ActorSet.h"

namespace siena
{

/**
 * Creates a constant dyadic covariate between a pair of actor sets.
 * @param[in] name the name of the covariate
 * @param[in] pFirstActorSet one of the involved actor sets
 * @param[in] pSecondActorSet the other of the involved actor sets
 */
ConstantDyadicCovariate::ConstantDyadicCovariate(std::string name,
	const ActorSet * pFirstActorSet,
	const ActorSet * pSecondActorSet) :
		DyadicCovariate(name, pFirstActorSet, pSecondActorSet)
{
	this->lpValues = new map<int, double>[pFirstActorSet->n()];
	this->lpMissings = new set<int>[pFirstActorSet->n()];
}


/**
 * Deallocates this covariate object.
 */
ConstantDyadicCovariate::~ConstantDyadicCovariate()
{
	delete[] this->lpValues;
	delete[] this->lpMissings;
	this->lpValues = 0;
	this->lpMissings = 0;
}


/**
 * Stores the value for the given pair of actors.
 * @param[in] i the first actor of the pair
 * @param[in] j the second actor of the pair
 * @param[in] value the value to be stored
 */
void ConstantDyadicCovariate::value(int i, int j, double value)
{
	if (value)
	{
		this->lpValues[i][j] = value;
	}
	else
	{
		this->lpValues[i].erase(j);
	}
}


/**
 * Returns the value of the covariate for the given pair of actors.
 */
double ConstantDyadicCovariate::value(int i, int j) const
{
	map<int, double>::const_iterator iter = this->lpValues[i].find(j);
	double value = 0;
	
	if (iter != this->lpValues[i].end())
	{
		value = iter->second;
	}
	
	return value;
}


/**
 * Stores if the value for the given pair of actors is missing.
 * @param[in] i the first actor of the pair
 * @param[in] j the second actor of the pair
 * @param[in] flag indicates if the value is missing
 */
void ConstantDyadicCovariate::missing(int i, int j, bool flag)
{
	if (flag)
	{
		this->lpMissings[i].insert(j);
	}
	else
	{
		this->lpMissings[i].erase(j);
	}
}


/**
 * Returns if the value for the given pair of actors is missing.
 */
bool ConstantDyadicCovariate::missing(int i, int j) const
{
	return this->lpMissings[i].find(j) != this->lpMissings[i].end();
}

}
