/******************************************************************************
 * SIENA: Simulation Investigation for Empirical Network Analysis
 *
 * Web: http://www.stats.ox.ac.uk/~snijders/siena/
 *
 * File: BalanceEffect.h
 *
 * Description:  This file contains the declaration of the class
 * BalanceEffect.
 *****************************************************************************/

#ifndef BALANCEEFFECT_H_
#define BALANCEEFFECT_H_

#include "NetworkEffect.h"
#include "data/IncidentTieIterator.h"

namespace siena
{

/**
 * This class defines the balance effect.
 */
class BalanceEffect : public NetworkEffect
{
public:
	BalanceEffect(const EffectInfo * pEffectInfo);

	virtual void initialize(EpochSimulation * pSimulation);
	virtual void initialize(const Data * pData, State * pState, int period);

	virtual double calculateTieFlipContribution(int alter) const;
	virtual bool usesTable(const ConfigurationTable * pTable) const;

protected:
	virtual double statistic(Network * pNetwork,
		Network * pSummationTieNetwork) const;

private:
	void markInvalidActors(IncidentTieIterator iter,
		int * flag,
		int markValue,
		int & validActorCount) const;

	// The centering constant b_0.
	double lbalanceMean;
};

}

#endif /*BALANCEEFFECT_H_*/
