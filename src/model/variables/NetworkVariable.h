/******************************************************************************
 * SIENA: Simulation Investigation for Empirical Network Analysis
 *
 * Web: http://www.stats.ox.ac.uk/~snijders/siena/
 *
 * File: NetworkVariable.h
 *
 * Description: This file contains the definition of the
 * NetworkVariable class.
 *****************************************************************************/

#ifndef NETWORKVARIABLE_H_
#define NETWORKVARIABLE_H_

#include <vector>

#include "DependentVariable.h"

namespace siena
{

// ----------------------------------------------------------------------------
// Section: Forward declarations
// ----------------------------------------------------------------------------

class Network;
class NetworkLongitudinalData;
class NetworkCache;


// ----------------------------------------------------------------------------
// Section: NetworkVariable class
// ----------------------------------------------------------------------------

/**
 * This class represents the state of a one-mode network variable.
 * @see DependentVariable
 */
class NetworkVariable : public DependentVariable
{
public:
	NetworkVariable(NetworkLongitudinalData * pData,
		EpochSimulation * pSimulation);
	virtual ~NetworkVariable();

	const SimulationActorSet * pSenders() const;
	const SimulationActorSet * pReceivers() const;
	virtual int m() const;
	virtual LongitudinalData * pData() const;
	bool oneModeNetwork() const;

	virtual void initialize(int period);
	virtual bool canMakeChange(int actor) const;
	virtual void makeChange(int actor);
	virtual void actOnJoiner(const SimulationActorSet * pActorSet, int actor);
	virtual void actOnLeaver(const SimulationActorSet * pActorSet, int actor);
	virtual void setLeaverBack(const SimulationActorSet * pActorSet,
		int actor);
	Network * pNetwork() const;

	int ego() const;

private:
	void preprocessEgo();
	void calculatePermissibleChanges();
	void calculateTieFlipContributions();
	void calculateTieFlipProbabilities();
	void accumulateScores(int alter) const;

	// The current state of the network
	Network * lpNetwork;

	// The observed data for this network variable
	NetworkLongitudinalData * lpData;

	// The set of actors acting as tie senders
	const SimulationActorSet * lpSenders;

	// The set of actors acting as tie receivers
	const SimulationActorSet * lpReceivers;

	// The number of structural tie variables to active alters per each actor.
	int * lactiveStructuralTieCount;

	// The current ego in the method makeChange
	int lego;

	// Indicates if a tie flip to a certain actor is permitted.
	bool * lpermitted;

	// A two-dimensional array of tie flip contributions to effects, where
	// rows correspond to alters and columns correspond to effects in the
	// evaluation function.

	double ** levaluationEffectContribution;

	// A two-dimensional array of tie flip contributions to effects, where
	// rows correspond to alters and columns correspond to effects in the
	// endowment function.

	double ** lendowmentEffectContribution;

	// Selection probability per each alter
	double * lprobabilities;

	// The cache object for repeated access of various structural properties
	// of this network during the simulation.

	NetworkCache * lpNetworkCache;
};


// ----------------------------------------------------------------------------
// Section: Inline accessors
// ----------------------------------------------------------------------------

/**
 * Returns the current ego making the change of the network.
 */
inline int NetworkVariable::ego() const
{
	return this->lego;
}

}

#endif /*NETWORKVARIABLE_H_*/
