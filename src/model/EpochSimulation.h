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
class Function;


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

    bool active(const ActorSet * pActorSet, int actor);
    int activeActorCount(const ActorSet * pActorSet);

    const Data * pData() const;
    const Model * pModel() const;
    const DependentVariable * pVariable(string name) const;
    const vector<DependentVariable *> & rVariables() const;
    int period() const;

    double time() const;

    double score(const EffectInfo * pEffect) const;
    void score(const EffectInfo * pEffect, double value);

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

    // A list of dependent variables with their current values
    vector<DependentVariable *> lvariables;
    
    // The current period to be simulated
    int lperiod;

    // An array of cummulative rates used for the random selection of
    // the dependent variable to change and the actor to make the change.

    double * lcummulativeRates;

    // A flag per each actor of each actor set indicating if the actor
    // has joined or not
    map<const ActorSet *, bool *> lactive;

    // The number of currently active actors per actor set
    map<const ActorSet *, int> lactiveActorCount;

    // The current time of the simmulation
    double ltime;

    // The current increment of time of the simmulation
    double ltau;

    // A sorted set of exogenous events of composition change
    const EventSet * lpEvents;

    // An iterator to the next event still to be processed.
    EventSet::iterator lnextEvent;

    // Target amount of change for this period if we are using conditional simulation
    int ltargetChange;

    // The dependent variable the simulation is conditioned upon
    DependentVariable * lpConditioningVariable;

    // Observed values of statistics in this simulation: one value for each
    // selected effect, including the rate effects, except basic rate for
    // conditioning variable, if any.
    vector<double> lsimulatedStatistics;

    // Values of scores in this simulation: one for each selected effect,
    // including the rate effects, but excluding the basic rate effect.

    map<const EffectInfo *, double> lscores;
};

}

#endif /*EPOCHSIMULATION_H_*/
