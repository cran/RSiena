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

    void runEpoch(int period);

    const Data * pData() const;
    const Model * pModel() const;
    const DependentVariable * pVariable(string name) const;
    const vector<DependentVariable *> & rVariables() const;
    const SimulationActorSet * pSimulationActorSet(
    	const ActorSet * pOriginalActorSet) const;
    int period() const;

    double time() const;

    double score(const EffectInfo * pEffect) const;
    void score(const EffectInfo * pEffect, double value);

    Cache * pCache() const;

private:
    void runStep();
    void calculateRates();
    void drawTimeIncrement();
    bool reachedCompositionChange() const;
    void makeNextCompositionChange();
	void setLeaversBack();
    DependentVariable * chooseVariable() const;
    int chooseActor(const DependentVariable * pVariable) const;
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

    // A list of dependent variables with their current values
    vector<DependentVariable *> lvariables;

    // The current period to be simulated
    int lperiod;

    // An array of cummulative rates used for the random selection of
    // the dependent variable to change and the actor to make the change.

    double * lcummulativeRates;

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

    State * lpState;
    Cache * lpCache;
};

}

#endif /*EPOCHSIMULATION_H_*/
