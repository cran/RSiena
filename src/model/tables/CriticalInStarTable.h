/******************************************************************************
 * SIENA: Simulation Investigation for Empirical Network Analysis
 * 
 * Web: http://www.stats.ox.ac.uk/~snijders/siena/
 * 
 * File: CriticalInStarTable.h
 * 
 * Description: This file defines the class CriticalInStarTable.
 *****************************************************************************/

#ifndef CRITICALINSTARTABLE_H_
#define CRITICALINSTARTABLE_H_

#include "ConfigurationTable.h"

namespace siena
{

/**
 * This class defines a configuration table storing the number of critical
 * in-stars between a fixed actor i, the ego, and each of the other actors.
 * An in-star <(i,h), (j,h)> is called critical, if there are no two-paths
 * from i to h through an intermediary actor other than j.
 */
class CriticalInStarTable : public ConfigurationTable
{
public:
	CriticalInStarTable(NetworkVariable * pVariable);

protected:
	virtual void vCalculate();
};

}

#endif /*CRITICALINSTARTABLE_H_*/
