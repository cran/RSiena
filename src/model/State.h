#ifndef STATE_H_
#define STATE_H_

#include <string>
#include <map>

using namespace std;

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
	State(Data * pData, int observation);
	State(EpochSimulation * pSimulation);

	const Network * pNetwork(string name) const;
	const int * behaviorValues(string name) const;
	const EpochSimulation * pSimulation() const;

private:
	map<string, const Network *> lnetworks;
	map<string, const int *> lbehaviors;
	EpochSimulation * lpSimulation;
};

}

#endif /*STATE_H_*/
