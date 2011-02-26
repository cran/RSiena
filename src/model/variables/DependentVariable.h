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
class MiniStep;


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
		EpochSimulation * pSimulation);
	virtual ~DependentVariable();

	void initializeRateFunction();
	void initializeEvaluationFunction();
	void initializeEndowmentFunction();

	inline const SimulationActorSet * pActorSet() const;
	int n() const;
	virtual int m() const = 0;
	virtual LongitudinalData * pData() const = 0;
	int id() const;
	virtual bool networkVariable() const;
	virtual bool behaviorVariable() const;
	virtual bool symmetric() const;
	virtual bool constrained() const;
	virtual int alter() const;

	inline const Function * pEvaluationFunction() const;
	inline const Function * pEndowmentFunction() const;

	virtual void initialize(int period);
	inline int period() const;
	virtual bool canMakeChange(int actor) const;
	virtual void makeChange(int actor) = 0;
	bool successfulChange() const;

	virtual void actOnJoiner(const SimulationActorSet * pActorSet,
		int actor);
	virtual void actOnLeaver(const SimulationActorSet * pActorSet,
		int actor);
	virtual void setLeaverBack(const SimulationActorSet * pActorSet,
		int actor) = 0;

	void calculateRates();
	double totalRate() const;
	double rate(int actor) const;
	inline double basicRate() const;

	int simulatedDistance() const;

	void accumulateRateScores(double tau,
		const DependentVariable * pSelectedVariable = 0,
		int selectedActor = 0);
	void accumulateRateScores(double tau,
		const DependentVariable * pSelectedVariable,
		int selectedActor, int alter);
	double basicRateScore() const;
	double constantCovariateScore(const ConstantCovariate * pCovariate) const;
	double changingCovariateScore(const ChangingCovariate * pCovariate) const;
	double behaviorVariableScore(const BehaviorVariable * pBehavior) const;
	double outDegreeScore(const NetworkVariable * pNetwork) const;
	double inDegreeScore(const NetworkVariable * pNetwork) const;
	double reciprocalDegreeScore(const NetworkVariable * pNetwork) const;
	double inverseOutDegreeScore(const NetworkVariable * pNetwork) const;

	// Maximum likelihood related

	/**
	 * Calculates the probability of the given ministep assuming that the
	 * ego of the ministep will change this variable.
	 */
	virtual double probability(MiniStep * pMiniStep) = 0;

	virtual bool validMiniStep(const MiniStep * pMiniStep,
		bool checkUpOnlyDownOnlyConditions = true) const;

	void updateEffectParameters();
	void updateEffectInfoParameters();

	/**
	 * Returns if the observed value for the option of the given ministep
	 * is missing at either end of the period.
	 */
	virtual bool missing(const MiniStep * pMiniStep) const = 0;

	/**
	 * Returns if the given ministep is structurally determined for the period
	 */
	virtual bool structural(const MiniStep * pMiniStep) const = 0;

	/**
	 * Generates a random ministep for the given ego.
	 */
	virtual MiniStep * randomMiniStep(int ego) = 0;

	void calculateMaximumLikelihoodRateScores(int activeMiniStepCount);
	void calculateMaximumLikelihoodRateDerivatives(int activeMiniStepCount);
	double basicRateDerivative() const;
	virtual double calculateChoiceProbability(const MiniStep * pMiniStep) const = 0;

	// Bayesian related: get new values
	void sampleBasicRate(int miniStepCount);
	double sampleParameters(double scaleFactor);
	// get and set sampled basic rates and shapes
	double sampledBasicRates(unsigned iteration) const;
	void sampledBasicRates(double value);
	int sampledBasicRatesDistributions(unsigned iteration) const;
	void sampledBasicRatesDistributions(int value);
	void clearSampledBasicRates();

protected:
	inline EpochSimulation * pSimulation() const;
	void simulatedDistance(int distance);
	void invalidateRates();
	void successfulChange(bool success);

private:
	void initializeFunction(Function * pFunction,
		const vector<EffectInfo *> & rEffects) const;

	bool constantRates() const;
	double calculateRate(int i);
	double structuralRate(int i) const;
	double behaviorVariableRate(int i) const;
	void updateCovariateRates();
	void calculateScoreSumTerms();

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

	// The derivative for the basic rate parameter for this variable for
	// this period
	double lbasicRateDerivative;

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

	// Sum term for scores for rate effects depending on constant covariates
	map<const ConstantCovariate *, double> lconstantCovariateSumTerm;

	// Sum term for scores for rate effects depending on changing covariates
	map<const ChangingCovariate *, double> lchangingCovariateSumTerm;

	// Sum term for scores for rate effects depending on behavior variables
	map<const BehaviorVariable *, double> lbehaviorVariableSumTerm;

	// Sum term for scores for rate effects depending on out degree
	map<const NetworkVariable *, double> loutDegreeSumTerm;

	// Sum term for scores for rate effects depending on in degree
	map<const NetworkVariable *, double> linDegreeSumTerm;

	// Sum term for scores for rate effects depending on reciprocal degree
	map<const NetworkVariable *, double> lreciprocalDegreeSumTerm;

	// Sum term for scores for rate effects depending on inverse degree
	map<const NetworkVariable *, double> linverseOutDegreeSumTerm;

	// Sum term for model B scores for rate effects depending on constant
	// covariates
	map<const ConstantCovariate *, double> lconstantCovariateModelBSumTerm;

	// Sum term for model B scores for rate effects depending on changing
	// covariates
	map<const ChangingCovariate *, double> lchangingCovariateModelBSumTerm;

	// Sum term for model B scores for rate effects depending on behavior
	// variables
	map<const BehaviorVariable *, double> lbehaviorVariableModelBSumTerm;

	// Sum term for model B scores for rate effects depending on out degree
	map<const NetworkVariable *, double> loutDegreeModelBSumTerm;

	// Sum term for model B scores for rate effects depending on inverse degree
	map<const NetworkVariable *, double> linverseOutDegreeModelBSumTerm;

	// Indicates if the rates are valid and shouldn't be recalculated
	// provided that the rates are constant during the period.

	int lvalidRates;

	// store for number of acceptances and rejections for non basic rate
	// parameters in Bayesian modelling

	int lacceptances;
	int lrejections;

	// flag to indicate we gave up on a step due to uponly and other filters
	bool lsuccessfulChange;
	// store for sampled parameters and the shapes used for basic rate
	// parameters
	vector<double> lsampledBasicRates;
	vector<int> lsampledBasicRatesDistributions;

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
