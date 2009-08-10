/******************************************************************************
 * SIENA: Simulation Investigation for Empirical Network Analysis
 *
 * Web: http://www.stats.ox.ac.uk/~snijders/siena/
 *
 * File: EffectInfo.cpp
 *
 * Description: This file contains the implementation of the
 * EffectInfo class.
 *****************************************************************************/

#include "EffectInfo.h"

namespace siena
{

/**
 * The constructor.
 * @param[in] variableName the name of the variable this effect is associated
 * with
 * @param[in] effectName the name of the effect
 * @param[in] effectType the type of the effect ("rate", "eval", or "endow")
 * @param[in] parameter the multiplicative weight of the effect
 * @param[in] internalEffectParameter the internal effect parameter
 * (if applicable)
 * @param[in] interactionName1 the name of another variable or covariate
 * (if any) the effect interacts with. If the effect interacts with two
 * variables or covariates, this parameter specifies one of them.
 * @param[in] interactionName2 the name of the other interacting variable or
 * covariate, if the effect has two such interactions. We should make sure
 * that the order of these interaction variables is the same both in R and C++.
 * @param[in] rateType distinguishes between structural rate effects and
 * covariate rate effects
 */
EffectInfo::EffectInfo(string variableName,
	string effectName,
	string effectType,
	double parameter,
	double internalEffectParameter,
	string interactionName1,
	string interactionName2,
	string rateType)
{
	this->lvariableName = variableName;
	this->leffectName = effectName;
	this->leffectType = effectType;
	this->lparameter = parameter;
	this->linternalEffectParameter = internalEffectParameter;
	this->linteractionName1 = interactionName1;
	this->linteractionName2 = interactionName2;
	this->lrateType = rateType;
}


// ----------------------------------------------------------------------------
// Section: Accessors
// ----------------------------------------------------------------------------

/**
 * Returns the name of the dependent variable this effect is associated with.
 */
string EffectInfo::variableName() const
{
	return this->lvariableName;
}


/**
 * Returns the parameter of this effect used in the respective function
 * (evaluation, endowment, or rate function).
 */
void EffectInfo::parameter(double value)
{
	this->lparameter = value;
}


/**
 * Sets the value of the parameter to be used as a multiplicative weight
 * for this effect, when calculating the respective function
 * (evaluation, endowment, or rate function).
 */
double EffectInfo::parameter() const
{
	return this->lparameter;
}
/**
 * Returns the value of the internaleffectparameter.
 */
double EffectInfo::internalEffectParameter() const
{
	return this->linternalEffectParameter;
}


/**
 * Returns the name of this effect specifying its semantics.
 */
string EffectInfo::effectName() const
{
	return this->leffectName;
}


/**
 * Returns the name of the variable or covariate this effect interacts with,
 * if any. If the effect has two such interactions, this method returns one
 * of them.
 */
string EffectInfo::interactionName1() const
{
	return this->linteractionName1;
}


/**
 * Returns the type of rate effect of this effect. Blank if not a rate effect,
 * values "structural" or "covariate".
 */
string EffectInfo::rateType() const
{
	return this->lrateType;
}

}
