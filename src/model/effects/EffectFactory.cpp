/******************************************************************************
 * SIENA: Simulation Investigation for Empirical Network Analysis
 *
 * Web: http://www.stats.ox.ac.uk/~snijders/siena/
 *
 * File: EffectFactory.cpp
 *
 * Description: This file contains the implementation of the
 * EffectFactory class.
 *****************************************************************************/

#include <stdexcept>

#include "EffectFactory.h"
#include "model/EffectInfo.h"
#include "model/effects/AllEffects.h"

namespace siena
{

/**
 * Creates and returns a concrete effect of the Effect class hierarchy
 * corresponding to the given generic effect descriptor.
 */
Effect * EffectFactory::createEffect(const EffectInfo * pEffectInfo) const
{
	Effect * pEffect = 0;
	string effectName = pEffectInfo->effectName();

	if (effectName == "density")
	{
		pEffect = new DensityEffect(pEffectInfo);
	}
	else if (effectName == "recip")
	{
		pEffect = new ReciprocityEffect(pEffectInfo);
	}
	else if (effectName == "transTrip")
	{
		pEffect = new TransitiveTripletsEffect(pEffectInfo);
	}
	else if (effectName == "transTriads")
	{
		pEffect = new TransitiveTriadsEffect(pEffectInfo);
	}
	else if (effectName == "transMedTrip")
	{
		pEffect = new TransitiveMediatedTripletsEffect(pEffectInfo);
	}
	else if (effectName == "cycle3")
	{
		pEffect = new ThreeCyclesEffect(pEffectInfo);
	}
	else if (effectName == "transTies")
	{
		pEffect = new TransitiveTiesEffect(pEffectInfo);
	}
	else if (effectName == "between")
	{
		pEffect = new BetweennessEffect(pEffectInfo);
	}
	else if (effectName == "balance")
	{
		pEffect = new BalanceEffect(pEffectInfo);
	}
	else if (effectName == "nbrDist2")
	{
		pEffect = new DistanceTwoEffect(pEffectInfo, 1);
	}
	else if (effectName == "nbrDist2twice")
	{
		pEffect = new DistanceTwoEffect(pEffectInfo, 2);
	}
//	else if (effectName == "denseTriads")
//	{
//		pEffect = new DenseTriadsEffect(pEffectInfo);
//	}
	else if (effectName == "inPop")
	{
		pEffect = new IndegreePopularityEffect(pEffectInfo, false);
	}
	else if (effectName == "inPopSqrt")
	{
		pEffect = new IndegreePopularityEffect(pEffectInfo, true);
	}
	else if (effectName == "outPop")
	{
		pEffect = new OutdegreePopularityEffect(pEffectInfo, false);
	}
	else if (effectName == "outPopSqrt")
	{
		pEffect = new OutdegreePopularityEffect(pEffectInfo, true);
	}
	else if (effectName == "inAct")
	{
		pEffect = new IndegreeActivityEffect(pEffectInfo, false);
	}
	else if (effectName == "inActSqrt")
	{
		pEffect = new IndegreeActivityEffect(pEffectInfo, true);
	}
	else if (effectName == "outAct")
	{
		pEffect = new OutdegreeActivityEffect(pEffectInfo);
	}
	else if (effectName == "outActSqrt")
	{
		pEffect = new OutdegreeActivitySqrtEffect(pEffectInfo);
	}
	else if (effectName == "outInv")
	{
		pEffect = new InverseOutdegreeEffect(pEffectInfo);
	}
	else if (effectName == "outSqInv")
	{
		pEffect = new InverseSquaredOutdegreeEffect(pEffectInfo);
	}
	else if (effectName == "outOutAss")
	{
		pEffect = new OutOutDegreeAssortativityEffect(pEffectInfo);
	}
	else if (effectName == "outInAss")
	{
		pEffect = new OutInDegreeAssortativityEffect(pEffectInfo);
	}
	else if (effectName == "inOutAss")
	{
		pEffect = new InOutDegreeAssortativityEffect(pEffectInfo);
	}
	else if (effectName == "inInAss")
	{
		pEffect = new InInDegreeAssortativityEffect(pEffectInfo);
	}
	else if (effectName == "X")
	{
		pEffect = new DyadicCovariateMainEffect(pEffectInfo);
	}
	else if (effectName == "Xrecip")
	{
		pEffect = new DyadicCovariateReciprocityEffect(pEffectInfo);
	}
	else if (effectName == "WWX")
	{
		pEffect = new WWXClosureEffect(pEffectInfo);
	}
	else if (effectName == "WXX")
	{
		pEffect = new WXXClosureEffect(pEffectInfo);
	}
	else if (effectName == "XWX")
	{
		pEffect = new XWXClosureEffect(pEffectInfo);
	}
	else if (effectName == "altX")
	{
		pEffect = new CovariateAlterEffect(pEffectInfo, false);
	}
	else if (effectName == "altSqX")
	{
		pEffect = new CovariateAlterEffect(pEffectInfo, true);
	}
	else if (effectName == "egoX")
	{
		pEffect = new CovariateEgoEffect(pEffectInfo);
	}
	else if (effectName == "simX")
	{
		pEffect = new CovariateSimilarityEffect(pEffectInfo, false);
	}
	else if (effectName == "simRecipX")
	{
		pEffect = new CovariateSimilarityEffect(pEffectInfo, true);
	}
	else if (effectName == "sameX")
	{
		pEffect = new SameCovariateEffect(pEffectInfo, false);
	}
	else if (effectName == "sameXRecip")
	{
		pEffect = new SameCovariateEffect(pEffectInfo, true);
	}
	else if (effectName == "egoXaltX")
	{
		pEffect = new CovariateEgoAlterEffect(pEffectInfo, false);
	}
	else if (effectName == "egoXaltXRecip")
	{
		pEffect = new CovariateEgoAlterEffect(pEffectInfo, true);
	}
//	else if (effectName == "IndTies")
//	{
//		pEffect = new CovariateIndirectTiesEffect(pEffectInfo);
//	}
	else if (effectName == "linear")
	{
		pEffect = new LinearShapeEffect(pEffectInfo);
	}
	else if (effectName == "quad")
	{
		pEffect = new QuadraticShapeEffect(pEffectInfo);
	}
	else if (effectName == "avSim")
	{
		pEffect = new AverageSimilarityEffect(pEffectInfo);
	}
	else if (effectName == "totSim")
	{
		pEffect = new TotalSimilarityEffect(pEffectInfo);
	}
	else if (effectName == "avAlt")
	{
		pEffect = new AverageAlterEffect(pEffectInfo);
	}
	else if (effectName == "indeg")
	{
		pEffect = new IndegreeEffect(pEffectInfo);
	}
	else if (effectName == "outdeg")
	{
		pEffect = new OutdegreeEffect(pEffectInfo);
	}
	else if (effectName == "effFrom")
	{
		pEffect = new ConstantCovariateMainBehaviorEffect(pEffectInfo);
	}
	else if (effectName == "effFromVar")
	{
		pEffect = new ChangingCovariateMainBehaviorEffect(pEffectInfo);
	}
	else if (effectName == "effFromBeh")
	{
		pEffect = new BehaviorMainBehaviorEffect(pEffectInfo);
	}
	else
	{
		throw domain_error("Unexpected effect name: " + effectName);
	}

	return pEffect;
}

}
