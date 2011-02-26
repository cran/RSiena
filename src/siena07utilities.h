/******************************************************************************
 * SIENA: Simulation Investigation for Empirical Network Analysis
 *
 * Web: http://www.stats.ox.ac.uk/~snijders/siena/
 *
 * File: siena07utilities.h
 *
 * Description: This file contains prototypes for various utilities used
 * for profiling, debugging, and creating R objects from C ones.
 *****************************************************************************/

#ifndef SIENA07UTILITIES_H_
#define SIENA07UTILITIES_H_
namespace siena
{
	class Data;
	class BehaviorVariable;
    class Network;
	class MiniStep;
	class Chain;
}
using namespace std;
using namespace siena;

//--------------------------------------
// utility functions to process groups
//-------------------------------------

/* Calculate the period number of this group and period, to access
 * correct stored chain.
 */
int periodFromStart(vector<Data *> & pGroupData, int group, int period);

/* Calculate the total number of periods in all groups, which is the dimension
 * of some returned arrays.
 */

int totalPeriods(vector<Data *> & pGroupData);
/**
 * Traps errors so R can stop the function rather than being stoppped itself.
 *
 */
void Rterminate();

/**
 * print out the data for profiling with gprof
 *
 */
void printOutData(Data *pData);

SEXP getBehaviorValues(const BehaviorVariable & behavior);
SEXP getAdjacency(const Network& net);
SEXP getEdgeList(const Network& net);
/**
 * utilities to access chains and ministeps
 *
 */
namespace siena
{
	SEXP getMiniStepDF(const MiniStep& miniStep);
	SEXP getChainDF(const Chain& chain);
	SEXP getMiniStepList(const MiniStep& miniStep, int period,
		const EpochSimulation& epochSimulation);
	SEXP getChainList(const Chain& chain, const EpochSimulation& epochSimulation);
}
#endif /*SIENA07UTILITIES_H_*/
