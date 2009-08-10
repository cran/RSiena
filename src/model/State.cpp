#include <stdexcept>

#include "State.h"
#include "data/Data.h"
#include "data/NetworkLongitudinalData.h"
#include "data/BehaviorLongitudinalData.h"
#include "model/EpochSimulation.h"
#include "model/variables/NetworkVariable.h"
#include "model/variables/BehaviorVariable.h"

namespace siena
{

State::State(Data * pData, int observation)
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
			this->lnetworks[pNetworkData->name()] =
				pNetworkData->pNetwork(observation);
		}
		else if (pBehaviorData)
		{
			this->lbehaviors[pBehaviorData->name()] =
				pBehaviorData->values(observation);
		}
		else
		{
			throw domain_error("Unexpected class of longitudinal data");
		}
	}
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
	lpSimulation = pSimulation;
}


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

const EpochSimulation * State::pSimulation() const
{

	return lpSimulation;
}

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

}
