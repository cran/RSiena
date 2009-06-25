/******************************************************************************
 * SIENA: Simulation Investigation for Empirical Network Analysis
 * 
 * Web: http://www.stats.ox.ac.uk/~snijders/siena/
 * 
 * File: ChangingDyadicCovariate.cpp
 * 
 * Description: This file contains the implementation of the
 * ChangingDyadicCovariate class.
 *****************************************************************************/

#include "ChangingDyadicCovariate.h"
#include "data/ActorSet.h"

namespace siena
{

/**
 * Creates a changing dyadic covariate between a pair of actor sets.
 * @param[in] name the name of the covariate
 * @param[in] pFirstActorSet one of the involved actor sets
 * @param[in] pSecondActorSet the other of the involved actor sets
 * @param[in] observationCount the number of observations of the covariate
 */
ChangingDyadicCovariate::ChangingDyadicCovariate(std::string name,
	const ActorSet * pFirstActorSet,
	const ActorSet * pSecondActorSet,
	int observationCount) :
		DyadicCovariate(name, pFirstActorSet, pSecondActorSet)
{
	this->lobservationCount = observationCount;
	this->lpValues = new map<int, double> * [observationCount];
	this->lpMissings = new set<int> * [observationCount];
	
	for (int k = 0; k < observationCount; k++)
	{
		this->lpValues[k] = new map<int, double>[pFirstActorSet->n()];
		this->lpMissings[k] = new set<int>[pFirstActorSet->n()];
	}
}


/**
 * Deallocates this covariate object.
 */
ChangingDyadicCovariate::~ChangingDyadicCovariate()
{
	for (int k = 0; k < this->lobservationCount; k++)
	{
		delete[] this->lpValues[k];
		delete[] this->lpMissings[k];
	}
	
	delete[] this->lpValues;
	delete[] this->lpMissings;
	this->lpValues = 0;	
	this->lpMissings = 0;
}


/**
 * Stores the value for the given pair of actors at the given observation.
 * @param[in] i the first actor of the pair
 * @param[in] j the second actor of the pair
 * @param[in] observation the number of the observation
 * @param[in] value the value to be stored
 */
void ChangingDyadicCovariate::value(int i,
	int j,
	int observation,
	double value)
{
	if (value)
	{
		this->lpValues[observation][i][j] = value;
	}
	else
	{
		this->lpValues[observation][i].erase(j);
	}
}


/**
 * Returns the value of the covariate for the given pair of actors at the
 * given observation.
 */
double ChangingDyadicCovariate::value(int i, int j, int observation) const
{
	map<int, double>::const_iterator iter =
		this->lpValues[observation][i].find(j);
	double value = 0;
	
	if (iter != this->lpValues[observation][i].end())
	{
		value = iter->second;
	}
	
	return value;
}


/**
 * Stores if the value for the given pair of actors is missing
 * at the given observation.
 * @param[in] i the first actor of the pair
 * @param[in] j the second actor of the pair
 * @param[in] observation the number of the observation
 * @param[in] flag indicates if the value is missing
 */
void ChangingDyadicCovariate::missing(int i,
	int j,
	int observation,
	bool flag)
{
	if (flag)
	{
		this->lpMissings[observation][i].insert(j);
	}
	else
	{
		this->lpMissings[observation][i].erase(j);
	}
}


/**
 * Returns if the value is missing for the given pair of actors at the
 * given observation.
 */
bool ChangingDyadicCovariate::missing(int i, int j, int observation) const
{
	return this->lpMissings[observation][i].find(j) !=
		this->lpMissings[observation][i].end();
}

}
