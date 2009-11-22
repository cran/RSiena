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
	State();

	const Network * pNetwork(string name) const;
	void pNetwork(string name, const Network * pNetwork);

	const int * behaviorValues(string name) const;
	void behaviorValues(string name, const int * values);

	void deleteValues();

private:
	map<string, const Network *> lnetworks;
	map<string, const int *> lbehaviors;
};

}

#endif /*STATE_H_*/
