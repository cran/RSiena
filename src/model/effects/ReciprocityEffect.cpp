/******************************************************************************
 * SIENA: Simulation Investigation for Empirical Network Analysis
 * 
 * Web: http://www.stats.ox.ac.uk/~snijders/siena/
 * 
 * File: ReciprocityEffect.cpp
 * 
 * Description: This file contains the implementation of the class
 * ReciprocityEffect.
 *****************************************************************************/

#include "ReciprocityEffect.h"
#include "data/OneModeNetwork.h"
#include "data/IncidentTieIterator.h"
#include "model/variables/NetworkVariable.h"

namespace siena
{

/**
 * Constructor.
 */
ReciprocityEffect::ReciprocityEffect(const EffectInfo * pEffectInfo) :
	NetworkEffect(pEffectInfo)
{
}


/**
 * Calculates the contribution of a tie flip to the given actor.
 */
double ReciprocityEffect::calculateTieFlipContribution(int alter) const
{
	double change = 0;
	
	// This tie flip has an effect only if there is a tie from the alter
	// to the ego.
	
	if (this->pVariable()->inTieExists(alter))
	{
		// Check if we gain or loose a reciprocated tie by doing this flip.
		
		change = 1;

		if (this->pVariable()->outTieExists(alter))
		{
			change = -1;
		}
	}
	
	return change;
}


/**
 * Returns the statistic corresponding to this effect as part of
 * the evaluation function with respect to the given network.
 */
double ReciprocityEffect::evaluationStatistic(Network * pNetwork) const
{
	OneModeNetwork * pOneModeNetwork = (OneModeNetwork *) pNetwork;
	int reciprocatedTieCount = 0;
	
	for (int i = 0; i < pOneModeNetwork->n(); i++)
	{
		reciprocatedTieCount += pOneModeNetwork->reciprocalDegree(i);
	}
	
	return reciprocatedTieCount;
}


/**
 * Returns the statistic corresponding to this effect as part of
 * the endowment function with respect to an initial network
 * and a network of lost ties. The current network is implicit as
 * the introduced ties are not relevant for calculating
 * endowment statistics.
 */
double ReciprocityEffect::endowmentStatistic(Network * pInitialNetwork,
	Network * pLostTieNetwork) const
{
	int n = pInitialNetwork->n();
	int * marks = new int[n];
	
	for (int i = 0; i < n; i++)
	{
		marks[i] = -1;
	}
	
	int counter = 0;
	
	for (int i = 0; i < n; i++)
	{
		for (IncidentTieIterator iter = pInitialNetwork->inTies(i);
			iter.valid();
			iter.next())
		{
			marks[iter.actor()] = i;
		}

		for (IncidentTieIterator iter = pLostTieNetwork->outTies(i);
			iter.valid();
			iter.next())
		{
			if (marks[iter.actor()] == i)
			{
				counter++;
			}
		}
	}
	
	return counter;
}

}
