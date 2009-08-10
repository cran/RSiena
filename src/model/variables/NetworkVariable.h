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
	inline ConfigurationTable * pRRTable() const;
	inline ConfigurationTable * pRFTable() const;
	inline ConfigurationTable * pRBTable() const;
	inline ConfigurationTable * pFRTable() const;
	inline ConfigurationTable * pBRTable() const;

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

	// The set of actors acting as tie senders
	const SimulationActorSet * lpSenders;

	// The set of actors acting as tie receivers
	const SimulationActorSet * lpReceivers;

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

	// The number of actors h with reciprocated ties to both i and j.
	ConfigurationTable * lpRRTable;

	// The number of actors h with a reciprocated tie to i and a tie to j.
	ConfigurationTable * lpRFTable;

	// The number of actors h with a reciprocated tie to i and a tie from j.
	ConfigurationTable * lpRBTable;

	// The number of actors h with a tie to i and a reciprocated tie to j.
	ConfigurationTable * lpFRTable;

	// The number of actors h with a tie from i and a reciprocated tie to j.
	ConfigurationTable * lpBRTable;

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
 * Returns the table storing the number of actors with reciprocated ties
 * to both i and j.
 */
inline ConfigurationTable * NetworkVariable::pRRTable() const
{
	return this->lpRRTable;
}


/**
 * Returns the table storing the number of actors with a reciprocated tie
 * to i and a tie to j.
 */
inline ConfigurationTable * NetworkVariable::pRFTable() const
{
	return this->lpRFTable;
}


/**
 * Returns the table storing the number of actors with a reciprocated tie
 * to i and a tie from j.
 */
inline ConfigurationTable * NetworkVariable::pRBTable() const
{
	return this->lpRBTable;
}


/**
 * Returns the table storing the number of actors with a tie to i and a
 * reciprocated tie to j.
 */
inline ConfigurationTable * NetworkVariable::pFRTable() const
{
	return this->lpFRTable;
}


/**
 * Returns the table storing the number of actors with a tie from i and a
 * reciprocated tie to j.
 */
inline ConfigurationTable * NetworkVariable::pBRTable() const
{
	return this->lpBRTable;
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
