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
class ConfigurationTable;


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

	const ActorSet * pSenders() const;
	const ActorSet * pReceivers() const;
	virtual int m() const;
	virtual LongitudinalData * pData() const;
	bool oneModeNetwork() const;

	virtual void initialize(int period);
	virtual bool canMakeChange(int actor) const;
	virtual void makeChange(int actor);
	virtual void actOnJoiner(const ActorSet * pActorSet, int actor);
	virtual void actOnLeaver(const ActorSet * pActorSet, int actor);
	virtual void setLeaverBack(const ActorSet * pActorSet, int actor);
	Network * pNetwork() const;
	Network * pPredictorNetwork() const;
	void pPredictorNetwork(Network *);

	int ego() const;
	bool outTieExists(int alter) const;
	bool inTieExists(int alter) const;
	inline ConfigurationTable * pTwoPathTable() const;
	inline ConfigurationTable * pReverseTwoPathTable() const;
	inline ConfigurationTable * pInStarTable() const;
	inline ConfigurationTable * pOutStarTable() const;
	inline ConfigurationTable * pCriticalInStarTable() const;

private:
	void preprocessEgo();
	bool required(const ConfigurationTable * pTable) const;
	void calculatePermissibleChanges();
	void calculateTieFlipContributions();
	void calculateTieFlipProbabilities();
	void accumulateScores(int alter) const;
	void initConfigurationTables();
	void deleteConfigurationTables();

	// The current state of the network
	Network * lpNetwork;

	// The observed data for this network variable
	NetworkLongitudinalData * lpData;

	// The number of structural tie variables to active alters per each actor.
	int * lactiveStructuralTieCount;

	// The current ego in the method makeChange
	int lego;

	// A list of various tables described below
	std::vector<ConfigurationTable *> lconfigurationTables;

	// The number of two-paths from the ego to each of the alters
	ConfigurationTable * lpTwoPathTable;

	// The number of two-paths from each of the alters to the ego
	ConfigurationTable * lpReverseTwoPathTable;

	// The number of in-stars between the ego and each of the alters.
	ConfigurationTable * lpInStarTable;

	// The number of out-stars between the ego and each of the alters.
	ConfigurationTable * lpOutStarTable;

	// The number of in-stars <(i,h), (j,h)> between the ego i and each
	// of the alters j, such that there are no two paths i -> h' -> h for
	// h' != j.

	ConfigurationTable * lpCriticalInStarTable;

	// Indicates if there is a tie from ego to each of the alters.
	bool * lpHasOutTie;

	// Indicates if there is a tie to ego from each of the alters.
	bool * lpHasInTie;

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

	// The network as at the start of the period, after removing missing
	// values. Used when this variable is a predictor in an effect for another
	// dependent variable.

	Network * lpPredictorNetwork;

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


/**
 * Returns the table storing the number of two-paths from the ego to
 * each of the alters.
 */
inline ConfigurationTable * NetworkVariable::pTwoPathTable() const
{
	return this->lpTwoPathTable;
}


/**
 * Returns the table storing the number of two-paths from each of the
 * alters to the ego.
 */
inline ConfigurationTable * NetworkVariable::pReverseTwoPathTable() const
{
	return this->lpReverseTwoPathTable;
}


/**
 * Returns the table storing the number of in-stars between the ego and
 * each of the alters.
 */
inline ConfigurationTable * NetworkVariable::pInStarTable() const
{
	return this->lpInStarTable;
}


/**
 * Returns the table storing the number of out-stars between the ego and
 * each of the alters.
 */
inline ConfigurationTable * NetworkVariable::pOutStarTable() const
{
	return this->lpOutStarTable;
}


/**
 * Returns the table storing the number of critical in-stars between the
 * ego and each of the alters. An in-star <(i,h), (j,h)> is critical if
 * there are no two paths i -> h' -> h for h' != j.
 */
inline ConfigurationTable * NetworkVariable::pCriticalInStarTable() const
{
	return this->lpCriticalInStarTable;
}


/**
 * Indicates if there is a tie from the ego to the given alter. This
 * method runs in constant time.
 */
inline bool NetworkVariable::outTieExists(int alter) const
{
	return this->lpHasOutTie[alter];
}


/**
 * Indicates if there is a tie from the given alter to the ego. This
 * method runs in constant time.
 */
inline bool NetworkVariable::inTieExists(int alter) const
{
	return this->lpHasInTie[alter];
}

}

#endif /*NETWORKVARIABLE_H_*/
