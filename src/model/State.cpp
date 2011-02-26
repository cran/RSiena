#include <stdexcept>

#include "State.h"
#include "data/Data.h"
#include "network/Network.h"
#include "data/NetworkLongitudinalData.h"
#include "data/BehaviorLongitudinalData.h"
#include "model/EpochSimulation.h"
#include "model/variables/NetworkVariable.h"
#include "model/variables/BehaviorVariable.h"

namespace siena
{

/**
 * Creates a state of variables as of the given observation of the given
 * Data object. The values may be copied of referenced directly as indicated
 * by the parameter <code>copyValues</code>.
 */
State::State(const Data * pData, int observation, bool copyValues)
{
	const vector<LongitudinalData *> & rVariables =
		pData->rDependentVariableData();

	for (unsigned i = 0; i < rVariables.size(); i++)
	{
		NetworkLongitudinalData * pNetworkData =
			dynamic_cast<NetworkLongitudinalData *>(rVariables[i]);
		BehaviorLongitudinalData * pBehaviorData =
			dynamic_cast<BehaviorLongitudinalData *>(rVariables[i]);

		if (pNetworkData)
		{
			const Network * pNetwork = pNetworkData->pNetwork(observation);

			if (copyValues)
			{
				pNetwork = pNetwork->clone();
			}

			this->lnetworks[pNetworkData->name()] = pNetwork;
		}
		else if (pBehaviorData)
		{
			const int * values = pBehaviorData->values(observation);

			if (copyValues)
			{
				int * copies = new int[pBehaviorData->n()];

				for (int actor = 0; actor < pBehaviorData->n(); actor++)
				{
					copies[actor] = values[actor];
				}

				values = copies;
			}

			this->lbehaviors[pBehaviorData->name()] = values;
		}
		else
		{
			throw domain_error("Unexpected class of longitudinal data");
		}
	}

	this->lownedValues = copyValues;
}


State::State(EpochSimulation * pSimulation)
{
	const vector<DependentVariable *> & rVariables =
		pSimulation->rVariables();

	for (unsigned i = 0; i < rVariables.size(); i++)
	{
		NetworkVariable * pNetworkVariable =
			dynamic_cast<NetworkVariable *>(rVariables[i]);
		BehaviorVariable * pBehaviorVariable =
			dynamic_cast<BehaviorVariable *>(rVariables[i]);

		if (pNetworkVariable)
		{
			this->lnetworks[pNetworkVariable->name()] =
				pNetworkVariable->pNetwork();
		}
		else if (pBehaviorVariable)
		{
			this->lbehaviors[pBehaviorVariable->name()] =
				pBehaviorVariable->values();
		}
		else
		{
			throw domain_error("Unexpected class of dependent variable");
		}
	}

	this->lownedValues = false;
}


/**
 * Default constructor creating an empty state. The current values of dependent
 * variables can be stored later with the appropriate setters.
 */
State::State()
{
	this->lownedValues = false;
}


/**
 * Deallocates this state.
 */
State::~State()
{
	if (this->lownedValues)
	{
		this->deleteValues();
	}
}


/**
 * Returns the network for the given name or 0 if there is no such a network
 * stored in this state.
 */
const Network * State::pNetwork(string name) const
{
	const Network * pNetwork = 0;
	map<string, const Network *>::const_iterator iter =
		this->lnetworks.find(name);

	if (iter != this->lnetworks.end())
	{
		pNetwork = iter->second;
	}

	return pNetwork;
}


/**
 * Stores the network for the given name.
 */
void State::pNetwork(string name, const Network * pNetwork)
{
	this->lnetworks[name] = pNetwork;
}


/**
 * Returns the values of the behavior variable with the given name, or 0
 * if no such values are stored in this state.
 */
const int * State::behaviorValues(string name) const
{
	const int * values = 0;
	map<string, const int *>::const_iterator iter =
		this->lbehaviors.find(name);

	if (iter != this->lbehaviors.end())
	{
		values = iter->second;
	}

	return values;
}


/**
 * Stores the values of a behavior variable with the given name.
 */
void State::behaviorValues(string name, const int * values)
{
	this->lbehaviors[name] = values;
}


/**
 * Deletes the values stored in this state.
 */
void State::deleteValues()
{
	// Cannot use clearMap as the keys are not pointers.

	while (!this->lnetworks.empty())
	{
		const Network * pNetwork = this->lnetworks.begin()->second;
		this->lnetworks.erase(this->lnetworks.begin());
		delete pNetwork;
	}

	// Cannot use clearMap as the values are arrays and should be deleted with
	// delete[] operator.

	while (!this->lbehaviors.empty())
	{
		const int * values = this->lbehaviors.begin()->second;
		this->lbehaviors.erase(this->lbehaviors.begin());
		delete[] values;
	}
}

}
