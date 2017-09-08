#ifndef STATE_H_
#define STATE_H_

#include <string>
#include <map>

namespace siena
{

class Network;
class Data;
class EpochSimulation;


/**
 * This class represents a state of dependent variables.
 */
class State
{
public:
	State(const Data * pData, int observation, bool copyValues = false);
	State(EpochSimulation * pSimulation);
	State();
	virtual ~State();

	const Network * pNetwork(std::string name) const;
	void pNetwork(std::string name, const Network * pNetwork);

	const int * behaviorValues(std::string name) const;
	void behaviorValues(std::string name, const int * values);

	void deleteValues();

private:
	std::map<std::string, const Network *> lnetworks;
	std::map<std::string, const int *> lbehaviors;
	bool lownedValues;
};

}

#endif /*STATE_H_*/
