/******************************************************************************
 * SIENA: Simulation Investigation for Empirical Network Analysis
 *
 * Web: http://www.stats.ox.ac.uk/~snijders/siena/
 *
 * File: EpochSimulation.h
 *
 * Description: This file contains the definition of the
 * EpochSimulation class.
 *****************************************************************************/

#ifndef EPOCHSIMULATION_H_
#define EPOCHSIMULATION_H_

#include <vector>
#include <map>
#include <string>
#include "data/Data.h"

using namespace std;

namespace siena
{

// ----------------------------------------------------------------------------
// Section: Forward declarations
// ----------------------------------------------------------------------------

class DependentVariable;
class Model;
class ActorSet;
class EffectInfo;
class SimulationActorSet;
class State;
class Cache;
class Chain;
class MiniStep;


// ----------------------------------------------------------------------------
// Section: EpochSimulation class
// ----------------------------------------------------------------------------

/**
 * This class provides the functionality necessary for simulating a model
 * between two observations.
 */
class EpochSimulation
{
public:
    EpochSimulation(Data * pData, Model * pModel);
    virtual ~EpochSimulation();

    void initialize(int period);

    // Method of moments related
    void runEpoch(int period);

    // Accessors

    const Data * pData() const;
    const Model * pModel() const;
    const DependentVariable * pVariable(string name) const;
    const vector<DependentVariable *> & rVariables() const;
    const SimulationActorSet * pSimulationActorSet(
    	const ActorSet * pOriginalActorSet) const;
    int period() const;
    double time() const;
    Cache * pCache() const;

    // Scores

    double score(const EffectInfo * pEffect) const;
    void score(const EffectInfo * pEffect, double value);
	map<const EffectInfo *, double>
		derivative(const EffectInfo * pEffect1) const;
	double derivative(const EffectInfo * pEffect1,
		const EffectInfo * pEffect2) const;
	void derivative(const EffectInfo * pEffect1, const EffectInfo * pEffect2,
		double value);
	Chain * pChain();
	double calculateChainProbabilities(Chain * chain);
	void updateParameters(int period);

protected:
    void calculateRates();
    double totalRate() const;
    DependentVariable * chooseVariable() const;
    int chooseActor(const DependentVariable * pVariable) const;

    // A vector of dependent variables with their current values
    vector<DependentVariable *> lvariables;

private:
    void runStep();
    void drawTimeIncrement();
    bool reachedCompositionChange() const;
    void makeNextCompositionChange();
	void setLeaversBack();
    void accumulateRateScores(double tau,
    	const DependentVariable * pSelectedVariable = 0,
    	int selectedActor = 0);

    // The observed data the model is based on
    Data * lpData;

    // The actor-based model to be simulated
    Model * lpModel;

    // A wrapper object per actor set for simulation purposes
    vector<SimulationActorSet *> lsimulationActorSets;

    // Stores the wrappers of each original actor set
    map<const ActorSet *, SimulationActorSet *> lactorSetMap;

    // The dependent variable for look-ups by variable names
    map<string, DependentVariable *> lvariableMap;

    // The current period to be simulated
    int lperiod;

    // An array of cummulative rates used for the random selection of
    // the dependent variable to change and the actor to make the change.

    double * lcummulativeRates;

    // The total rate over all dependent variables
    double ltotalRate;

    // The current time of the simmulation
    double ltime;

    // The current increment of time of the simmulation
    double ltau;

    // A sorted set of exogenous events of composition change
    const EventSet * lpEvents;

    // An iterator to the next event still to be processed.
    EventSet::const_iterator lnextEvent;

    // Target amount of change for this period if we are using conditional simulation
    int ltargetChange;

    // The dependent variable the simulation is conditioned upon
    DependentVariable * lpConditioningVariable;

    // Values of scores in this simulation: one for each selected effect,
    // including the rate effects, but excluding the basic rate effect.

    map<const EffectInfo *, double> lscores;
    map<const EffectInfo *, map <const EffectInfo *, double> > lderivatives;

    State * lpState;
    Cache * lpCache;

	Chain * lpChain;

};

}

#endif /*EPOCHSIMULATION_H_*/
