/******************************************************************************
 * SIENA: Simulation Investigation for Empirical Network Analysis
 *
 * Web: http://www.stats.ox.ac.uk/~snijders/siena/
 *
 * File: NetworkEffect.h
 *
 * Description: This file contains the definition of the
 * NetworkEffect class.
 *****************************************************************************/

#ifndef NETWORKEFFECT_H_
#define NETWORKEFFECT_H_

#include "Effect.h"

namespace siena
{

// ----------------------------------------------------------------------------
// Section: Forward declarations
// ----------------------------------------------------------------------------

class Network;
class NetworkVariable;
class NetworkLongitudinalData;
class ConfigurationTable;


// ----------------------------------------------------------------------------
// Section: NetworkEffect class
// ----------------------------------------------------------------------------

/**
 * The base class for all network effects.
 */
class NetworkEffect : public Effect
{
public:
	NetworkEffect(const EffectInfo * pEffectInfo);

	virtual void initialize(EpochSimulation * pSimulation);
	virtual void initialize(const Data * pData, State * pState, int period);

	inline const NetworkVariable * pVariable() const;
	inline const NetworkLongitudinalData * pData() const;
	virtual void preprocessEgo();

	/**
	 * Assuming that the ego would flip the tie to the given actor,
	 * this method calculates the change in the statistic corresponding
	 * to this effect. The method has to be overriden by all concrete
	 * effect classes.
	 */
	virtual double calculateTieFlipContribution(int alter) const = 0;

	virtual bool usesTable(const ConfigurationTable * pTable) const;

	virtual double evaluationStatistic(Network * pNetwork) const;
	virtual double endowmentStatistic(Network * pInitialNetwork,
		Network * pLostTieNetwork) const;

protected:
	virtual double statistic(Network * pNetwork,
		Network * pSummationTieNetwork) const;

private:
	// The network variable this effect is associated with
	const NetworkVariable * lpVariable;

	// The observed network data underlying the network variable.
	const NetworkLongitudinalData * lpNetworkData;
};


// ----------------------------------------------------------------------------
// Section: Inline methods
// ----------------------------------------------------------------------------

/**
 * Returns the network variable this effect is associated with.
 */
const NetworkVariable * NetworkEffect::pVariable() const
{
	return this->lpVariable;
}


/**
 * Returns the observed network data corresponding to the network variable of
 * this effect.
 */
const NetworkLongitudinalData * NetworkEffect::pData() const
{
	return this->lpNetworkData;
}

}

#endif /*NETWORKEFFECT_H_*/
