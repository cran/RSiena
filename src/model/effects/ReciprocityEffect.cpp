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
#include "network/OneModeNetwork.h"
#include "network/IncidentTieIterator.h"
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
double ReciprocityEffect::calculateContribution(int alter) const
{
	double change = 0;

	// This tie flip has an effect only if there is a tie from the alter
	// to the ego.

	if (this->inTieExists(alter))
	{
		change = 1;
	}

	return change;
}


/**
 * See base class.
 */
double ReciprocityEffect::statistic(const Network * pSummationTieNetwork) const
{
	const Network * pNetwork = this->pNetwork();
	int n = pNetwork->n();
	int * marks = new int[n];

	for (int i = 0; i < n; i++)
	{
		marks[i] = -1;
	}

	int counter = 0;

	for (int i = 0; i < n; i++)
	{
		for (IncidentTieIterator iter = pNetwork->inTies(i);
			iter.valid();
			iter.next())
		{
			marks[iter.actor()] = i;
		}

		for (IncidentTieIterator iter = pSummationTieNetwork->outTies(i);
			iter.valid();
			iter.next())
		{
			if (marks[iter.actor()] == i)
			{
				counter++;
			}
		}
	}

	delete[] marks;
	return counter;
}

}
