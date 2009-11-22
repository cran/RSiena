/******************************************************************************
 * SIENA: Simulation Investigation for Empirical Network Analysis
 *
 * Web: http://www.stats.ox.ac.uk/~snijders/siena/
 *
 * File: DependentVariable.h
 *
 * Description: This file contains the definition of the
 * DependentVariable class.
 *****************************************************************************/

#ifndef DEPENDENTVARIABLE_H_
#define DEPENDENTVARIABLE_H_

#include <map>
#include <string>
#include "utils/NamedObject.h"
#include "model/Function.h"

using namespace std;

namespace siena
{

// ----------------------------------------------------------------------------
// Section: Forward declarations
// ----------------------------------------------------------------------------

class ConstantCovariate;
class ChangingCovariate;
class NetworkVariable;
class BehaviorVariable;
class EffectValueTable;
class EpochSimulation;
class ActorSet;
class SimulationActorSet;
class LongitudinalData;
class BehaviorLongitudinalData;
class NetworkLongitudinalData;
class Network;
class EffectInfo;
class StructuralRateEffect;


// ----------------------------------------------------------------------------
// Section: DependentVariable class
// ----------------------------------------------------------------------------

/**
 * This class represents a certain dependent variable. It is the base class of
 * NetworkVariable and BehaviorVariable.
 * The class stores the current state of the variable and provides methods
 * supporting simulations of actor-based models.
 */
class DependentVariable : public NamedObject
{
public:
	DependentVariable(string name,
		const ActorSet * pActorSet,
		int observationCount,
		EpochSimulation * pSimulation);
	virtual ~DependentVariable();

	void initializeRateFunction();
	void initializeEvaluationFunction();
	void initializeEndowmentFunction();

	inline const SimulationActorSet * pActorSet() const;
	int n() const;
	virtual int m() const = 0;
	virtual LongitudinalData * pData() const = 0;

	inline const Function * pEvaluationFunction() const;
	inline const Function * pEndowmentFunction() const;

	virtual void initialize(int period);
	inline int period() const;
	virtual bool canMakeChange(int actor) const;
	virtual void makeChange(int actor) = 0;

	virtual void actOnJoiner(const SimulationActorSet * pActorSet,
		int actor) = 0;
	virtual void actOnLeaver(const SimulationActorSet * pActorSet,
		int actor) = 0;
	virtual void setLeaverBack(const SimulationActorSet * pActorSet,
		int actor) = 0;

	void calculateRates();
	double totalRate() const;
	double rate(int actor) const;

	int simulatedDistance() const;

	void accumulateRateScores(double tau,
		const DependentVariable * pSelectedVariable = 0,
		int selectedActor = 0);
	double basicRateScore() const;
	double constantCovariateScore(const ConstantCovariate * pCovariate) const;
	double changingCovariateScore(const ChangingCovariate * pCovariate) const;
	double behaviorVariableScore(const BehaviorVariable * pBehavior) const;
	double outDegreeScore(const NetworkVariable * pNetwork) const;
	double inDegreeScore(const NetworkVariable * pNetwork) const;
	double reciprocalDegreeScore(const NetworkVariable * pNetwork) const;
	double inverseOutDegreeScore(const NetworkVariable * pNetwork) const;

protected:
	inline EpochSimulation * pSimulation() const;
	void simulatedDistance(int distance);

private:
	void initializeFunction(Function * pFunction,
		const vector<EffectInfo *> & rEffects) const;

	virtual double calculateRate(int i);
	double structuralRate(int i) const;
	void updateCovariateRates();
	inline double basicRate() const;

	// A simulation of the actor-based model, which owns this variable
	EpochSimulation * lpSimulation;

	// The underlying set of actors
	const SimulationActorSet * lpActorSet;

	// The current period (in [0, observations - 2])
	int lperiod;

	// The total rate of change summed over all actors
	double ltotalRate;

	// The rate of change for each actor
	double * lrate;

	// The basic rate parameter for the current period
	double lbasicRate;

	// The covariate-based component of the rate function per each actor
	double * lcovariateRates;

	// Parameters for rate effects depending on constant covariates
	map<const ConstantCovariate *, double> lconstantCovariateParameters;

	// Parameters for rate effects depending on changing covariates
	map<const ChangingCovariate *, double> lchangingCovariateParameters;

	// Parameters for rate effects depending on behavior variables
	map<const BehaviorVariable *, double> lbehaviorVariableParameters;

	// The structural rate effects. Currently, there are four types of
	// structural rate effects, namely, the out-degree, in-degree,
	// reciprocal degree, and inverse out-degree effects.

	vector<StructuralRateEffect *> lstructuralRateEffects;

	// The evaluation function for this variable
	Function * lpEvaluationFunction;

	// The endowment function for this variable
	Function * lpEndowmentFunction;

	// The distance of this variable from the observed data at the beginning
	// of the current period

	int lsimulatedDistance;

	// The score for the basic rate parameter for this variable for this period
	double lbasicRateScore;

	// Scores for rate effects depending on constant covariates
	map<const ConstantCovariate *, double> lconstantCovariateScores;

	// Scores for rate effects depending on changing covariates
	map<const ChangingCovariate *, double> lchangingCovariateScores;

	// Scores for rate effects depending on behavior variables
	map<const BehaviorVariable *, double> lbehaviorVariableScores;

	// Scores for rate effects depending on out degree
	map<const NetworkVariable *, double> loutDegreeScores;

	// Scores for rate effects depending on in degree
	map<const NetworkVariable *, double> linDegreeScores;

	// Scores for rate effects depending on reciprocal degree
	map<const NetworkVariable *, double> lreciprocalDegreeScores;

	// Scores for rate effects depending on inverse degree
	map<const NetworkVariable *, double> linverseOutDegreeScores;
};


// ----------------------------------------------------------------------------
// Section: Inline methods
// ----------------------------------------------------------------------------

/**
 * Returns the actor-based model owning this variable.
 */
EpochSimulation * DependentVariable::pSimulation() const
{
	return this->lpSimulation;
}


/**
 * Returns the set of actors underlying this dependent variable.
 */
const SimulationActorSet * DependentVariable::pActorSet() const
{
	return this->lpActorSet;
}


/**
 * Returns the evaluation function of this variable.
 */
const Function * DependentVariable::pEvaluationFunction() const
{
	return this->lpEvaluationFunction;
}


/**
 * Returns the endowment function of this variable.
 */
const Function * DependentVariable::pEndowmentFunction() const
{
	return this->lpEndowmentFunction;
}


/**
 * Returns the index of the current period.
 */
int DependentVariable::period() const
{
	return this->lperiod;
}


/**
 * Returns the basic rate parameter for the current period.
 */
double DependentVariable::basicRate() const
{
	return this->lbasicRate;
}

}

#endif /*DEPENDENTVARIABLE_H_*/
