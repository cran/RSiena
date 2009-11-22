/******************************************************************************
 * SIENA: Simulation Investigation for Empirical Network Analysis
 *
 * Web: http://www.stats.ox.ac.uk/~snijders/siena/
 *
 * File: EffectInfo.h
 *
 * Description: This file contains the definition of the
 * EffectInfo class.
 *****************************************************************************/

#ifndef EFFECTINFO_H_
#define EFFECTINFO_H_

#include <string>

using namespace std;

namespace siena
{

/**
 * Encapsulates the information relevant to an effect, except for the
 * basic rate effects.
 */
class EffectInfo
{
public:
	EffectInfo(string variableName,
			   string effectName,
			   string effectType,
			   double parameter,
			   double internalEffectParameter,
			   string interactionName1,
			   string interactionName2,
			   string rateType);

	void parameter(double value);

	string variableName() const;
	string effectName() const;
	double parameter() const;
	double internalEffectParameter() const;
	string interactionName1() const;
	string rateType() const;

private:
	// The name of the variable this effect is associated with
	string lvariableName;

	// A short name of the effect used to identify the semantics
	// of this effect

	string leffectName;

	// The type of the effect ("rate", "eval", or "endow")
	string leffectType;

	// The multiplicative weight in the respective function
	double lparameter;

	// The internal parameter, if applicable
	double linternalEffectParameter;

	// The name of a variable or covariate this effect is interacting with,
	// if applicable

	string linteractionName1;

	// The name of the other interaction variable or covariate, if the
	// effect has two such interactions

	string linteractionName2;

	// Distinguishes between structural rate effects and covariate rate effects
	string lrateType;
};

}

#endif /*EFFECTINFO_H_*/
