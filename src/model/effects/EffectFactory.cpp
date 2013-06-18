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
#include <cmath>

#include "EffectFactory.h"
#include "data/Data.h"
#include "data/NetworkLongitudinalData.h"
#include "model/EffectInfo.h"
#include "model/effects/AllEffects.h"
#include "model/effects/generic/GenericNetworkEffect.h"
#include "model/effects/generic/OutTieFunction.h"
#include "model/effects/generic/InTieFunction.h"
#include "model/effects/generic/ProductFunction.h"
#include "model/effects/generic/ConstantFunction.h"
#include "model/effects/generic/InDegreeFunction.h"
#include "model/effects/generic/IntSqrtFunction.h"
#include "model/effects/generic/DifferenceFunction.h"
#include "model/effects/generic/EgoInDegreeFunction.h"
#include "model/effects/generic/OutDegreeFunction.h"
#include "model/effects/generic/EgoOutDegreeFunction.h"
#include "model/effects/generic/BetweennessFunction.h"
#include "model/effects/generic/GwespFunction.h"
#include "model/effects/generic/InStarFunction.h"
#include "model/effects/generic/OutStarFunction.h"
#include "model/effects/generic/ReciprocatedTwoPathFunction.h"
#include "model/effects/generic/TwoPathFunction.h"
#include "model/effects/generic/ReverseTwoPathFunction.h"
#include "model/effects/generic/MixedTwoPathFunction.h"
#include "model/effects/generic/ConditionalFunction.h"
#include "model/effects/generic/EqualCovariatePredicate.h"
#include "model/effects/generic/MissingCovariatePredicate.h"
#include "model/effects/generic/CovariateDistance2AlterNetworkFunction.h"
#include "model/effects/generic/CovariateDistance2SimilarityNetworkFunction.h"
#include "model/effects/generic/CovariateMixedNetworkAlterFunction.h"
#include "model/effects/generic/SameCovariateTwoPathFunction.h"
#include "model/effects/generic/SameCovariateMixedTwoPathFunction.h"


#include "model/tables/EgocentricConfigurationTable.h"
#include "model/tables/NetworkCache.h"

namespace siena
{

/**
 * Constructor.
 * @param[in] pData the data this factory will create effects for
 */
EffectFactory::EffectFactory(const Data * pData)
{
	this->lpData = pData;
}


/**
 * Creates and returns a concrete effect of the Effect class hierarchy
 * corresponding to the given generic effect descriptor.
 */
Effect * EffectFactory::createEffect(const EffectInfo * pEffectInfo) const
{
	Effect * pEffect = 0;
	string effectName = pEffectInfo->effectName();

	// Handle the user-defined interaction effects first.

	if (pEffectInfo->pEffectInfo1())
	{
		// The info object of the first interacting effect is defined,
		// which means that we have a user-defined interaction effect.

		Effect *pEffect1 = this->createEffect(pEffectInfo->pEffectInfo1());
		Effect *pEffect2 = this->createEffect(pEffectInfo->pEffectInfo2());
		Effect *pEffect3 = 0;
		if (pEffectInfo->pEffectInfo3())
		{
			pEffect3 = this->createEffect(pEffectInfo->pEffectInfo3());
		}


		NetworkEffect * pNetworkEffect1 =
			dynamic_cast<NetworkEffect *>(pEffect1);
		BehaviorEffect * pBehaviorEffect1 =
			dynamic_cast<BehaviorEffect *>(pEffect1);

		if (pNetworkEffect1)
		{
			NetworkEffect * pNetworkEffect2 =
				dynamic_cast<NetworkEffect *>(pEffect2);
			NetworkEffect * pNetworkEffect3 = 0;
			if (pEffect3)
			{
				 pNetworkEffect3 = dynamic_cast<NetworkEffect *>(pEffect3);
			}
			pEffect = new NetworkInteractionEffect(pEffectInfo,
				pNetworkEffect1,
				pNetworkEffect2,
				pNetworkEffect3);
		}
		else
		{
			BehaviorEffect * pBehaviorEffect2 =
				dynamic_cast<BehaviorEffect *>(pEffect2);
			BehaviorEffect * pBehaviorEffect3 = 0;
			if (pEffect3)
			{
				pBehaviorEffect3 = dynamic_cast<BehaviorEffect *>(pEffect3);
			}

			pEffect = new BehaviorInteractionEffect(pEffectInfo,
				pBehaviorEffect1,
				pBehaviorEffect2,
				pBehaviorEffect3);
		}
	}
	else if (effectName == "density")
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
	else if (effectName == "transRecTrip")
	{
		pEffect = new TransitiveReciprocatedTripletsEffect(pEffectInfo);
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
	else if (effectName == "denseTriads")
	{
		pEffect = new DenseTriadsEffect(pEffectInfo);
	}
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
	else if (effectName == "outTrunc")
	{
		pEffect = new TruncatedOutdegreeEffect(pEffectInfo);
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
	else if (effectName == "XRecip")
	{
		pEffect = new DyadicCovariateReciprocityEffect(pEffectInfo);
	}
	else if (effectName == "WWX")
	{
		pEffect = new WWXClosureEffect(pEffectInfo, true, true);
	}
	else if (effectName == "cyWWX")
	{
		pEffect = new WWXClosureEffect(pEffectInfo, false, false);
	}
	else if (effectName == "InWWX")
	{
		pEffect = new WWXClosureEffect(pEffectInfo, false, true);
	}
	else if (effectName == "OutWWX")
	{
		pEffect = new WWXClosureEffect(pEffectInfo, true, false);
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
	else if (effectName == "simXTransTrip")
	{
		pEffect = new SimilarityTransitiveTripletsEffect(pEffectInfo, false);
	}
	else if (effectName == "sameX")
	{
		pEffect = new SameCovariateEffect(pEffectInfo, false);
	}
	else if (effectName == "higher")
	{
		pEffect = new HigherCovariateEffect(pEffectInfo);
	}
	else if (effectName == "sameXRecip")
	{
		pEffect = new SameCovariateEffect(pEffectInfo, true);
	}
	else if (effectName == "sameXTransTrip")
	{
		pEffect = new SameCovariateTransitiveTripletsEffect(pEffectInfo, false);
	}
	else if (effectName == "jumpXTransTrip")
	{
		pEffect = new JumpCovariateTransitiveTripletsEffect(pEffectInfo, false);
	}
	else if (effectName == "egoXaltX")
	{
		pEffect = new CovariateEgoAlterEffect(pEffectInfo, false);
	}
	else if (effectName == "egoXaltXRecip")
	{
		pEffect = new CovariateEgoAlterEffect(pEffectInfo, true);
	}
	else if (effectName == "IndTies")
	{
		pEffect = new CovariateIndirectTiesEffect(pEffectInfo);
	}
	else if (effectName == "cycle4")
	{
		pEffect = new FourCyclesEffect(pEffectInfo);
	}
	else if (effectName == "gwespFF")
	{
		EgocentricConfigurationTable * (NetworkCache::*mytable)() const =
			&NetworkCache::pTwoPathTable;
		GwespFunction * pFunction =
			new GwespFunction(pEffectInfo->variableName(),
				mytable, pEffectInfo->internalEffectParameter());

 		pEffect = new GenericNetworkEffect(pEffectInfo,
			pFunction);
	}
	else if (effectName == "gwespFB")
 	{
		EgocentricConfigurationTable * (NetworkCache::*mytable)() const =
			&NetworkCache::pInStarTable;
		GwespFunction * pFunction =
			new GwespFunction(pEffectInfo->variableName(),
				mytable, pEffectInfo->internalEffectParameter());

		pEffect = new GenericNetworkEffect(pEffectInfo,
			pFunction);
	}
	else if (effectName == "gwespBF")
	{
		EgocentricConfigurationTable * (NetworkCache::*mytable)() const =
			&NetworkCache::pOutStarTable;
		GwespFunction * pFunction =
			new GwespFunction(pEffectInfo->variableName(),
				mytable, pEffectInfo->internalEffectParameter());
 		pEffect = new GenericNetworkEffect(pEffectInfo,
			pFunction);
	}
	else if (effectName == "gwespBB")
	{
		EgocentricConfigurationTable * (NetworkCache::*mytable)() const =
			&NetworkCache::pReverseTwoPathTable;
		GwespFunction * pFunction =
			new GwespFunction(pEffectInfo->variableName(),
				mytable, pEffectInfo->internalEffectParameter());
 		pEffect = new GenericNetworkEffect(pEffectInfo,
			pFunction);
	}
	else if (effectName == "gwespRR")
	{
		EgocentricConfigurationTable * (NetworkCache::*mytable)() const =
			&NetworkCache::pRRTable;
		GwespFunction * pFunction =
			new GwespFunction(pEffectInfo->variableName(),
				mytable, pEffectInfo->internalEffectParameter());
 		pEffect = new GenericNetworkEffect(pEffectInfo,
			pFunction);
	}
	else if (effectName == "inStructEq")
	{
		pEffect = new InStructuralEquivalenceEffect(pEffectInfo);
	}
	else if (effectName == "crprod")
	{
		pEffect = new GenericNetworkEffect(pEffectInfo,
			new OutTieFunction(pEffectInfo->interactionName1()));
	}
	else if (effectName == "crprodRecip")
	{
		pEffect = new GenericNetworkEffect(pEffectInfo,
			new InTieFunction(pEffectInfo->interactionName1()));
	}
	else if (effectName == "crprodMutual")
	{
		pEffect = new GenericNetworkEffect(pEffectInfo,
			new ProductFunction(
				new OutTieFunction(pEffectInfo->interactionName1()),
				new InTieFunction(pEffectInfo->interactionName1())));
	}
	else if (effectName == "inPopIntn")
	{
		AlterFunction * pFirstFunction =
			new InDegreeFunction(pEffectInfo->interactionName1());
		ConstantFunction * pSecondFunction =
			new ConstantFunction(pEffectInfo->interactionName1(),
				AVERAGE_IN_DEGREE);

		if (pEffectInfo->internalEffectParameter() == 2)
		{
			pFirstFunction = new IntSqrtFunction(pFirstFunction);
			pSecondFunction->pFunction(sqrt);
		}

		pEffect = new GenericNetworkEffect(pEffectInfo,
			new DifferenceFunction(pFirstFunction, pSecondFunction));
	}
	else if (effectName == "inActIntn")
	{
		AlterFunction * pFirstFunction =
			new EgoInDegreeFunction(pEffectInfo->interactionName1());
		ConstantFunction * pSecondFunction =
			new ConstantFunction(pEffectInfo->interactionName1(),
				AVERAGE_IN_DEGREE);

		if (pEffectInfo->internalEffectParameter() == 2)
		{
			pFirstFunction = new IntSqrtFunction(pFirstFunction);
			pSecondFunction->pFunction(sqrt);
		}

		pEffect = new GenericNetworkEffect(pEffectInfo,
			new DifferenceFunction(pFirstFunction, pSecondFunction));
	}
	else if (effectName == "outPopIntn")
	{
		AlterFunction * pFirstFunction =
			new OutDegreeFunction(pEffectInfo->interactionName1());
		ConstantFunction * pSecondFunction =
			new ConstantFunction(pEffectInfo->interactionName1(),
				AVERAGE_OUT_DEGREE);

		if (pEffectInfo->internalEffectParameter() == 2)
		{
			pFirstFunction = new IntSqrtFunction(pFirstFunction);
			pSecondFunction->pFunction(sqrt);
		}

		pEffect = new GenericNetworkEffect(pEffectInfo,
			new DifferenceFunction(pFirstFunction, pSecondFunction));
	}
	else if (effectName == "outActIntn")
	{
		AlterFunction * pFirstFunction =
			new EgoOutDegreeFunction(pEffectInfo->interactionName1());
		ConstantFunction * pSecondFunction =
			new ConstantFunction(pEffectInfo->interactionName1(),
				AVERAGE_OUT_DEGREE);

		if (pEffectInfo->internalEffectParameter() == 2)
		{
			pFirstFunction = new IntSqrtFunction(pFirstFunction);
			pSecondFunction->pFunction(sqrt);
		}

		pEffect = new GenericNetworkEffect(pEffectInfo,
			new DifferenceFunction(pFirstFunction, pSecondFunction));
	}
	else if (effectName == "both")
	{
		AlterFunction * pAlterIndegreeFunction =
			new InDegreeFunction(pEffectInfo->interactionName1());
		AlterFunction * pEgoIndegreeFunction =
			new EgoInDegreeFunction(pEffectInfo->interactionName1());
		ConstantFunction * pFirstConstantFunction =
			new ConstantFunction(pEffectInfo->interactionName1(),
				AVERAGE_IN_DEGREE);
		ConstantFunction * pSecondConstantFunction =
			new ConstantFunction(pEffectInfo->interactionName1(),
				AVERAGE_IN_DEGREE);

		if (pEffectInfo->internalEffectParameter() == 2)
		{
			pAlterIndegreeFunction =
				new IntSqrtFunction(pAlterIndegreeFunction);
			pEgoIndegreeFunction =
				new IntSqrtFunction(pEgoIndegreeFunction);
			pFirstConstantFunction->pFunction(sqrt);
			pSecondConstantFunction->pFunction(sqrt);
		}

		pEffect = new GenericNetworkEffect(pEffectInfo,
			new ProductFunction(
				new DifferenceFunction(pAlterIndegreeFunction,
					pFirstConstantFunction),
				new DifferenceFunction(pEgoIndegreeFunction,
					pSecondConstantFunction)));
	}
	else if (effectName == "betweenPop")
	{
		AlterFunction * pFunction =
			new BetweennessFunction(pEffectInfo->interactionName1());

		if (pEffectInfo->internalEffectParameter() == 2)
		{
			pFunction = new IntSqrtFunction(pFunction);
		}

		pEffect = new GenericNetworkEffect(pEffectInfo, pFunction);
	}
	else if (effectName == "from")
	{
		pEffect = new GenericNetworkEffect(pEffectInfo,
			new InStarFunction(pEffectInfo->interactionName1()));
	}
	else if (effectName == "fromMutual")
	{
		pEffect = new GenericNetworkEffect(pEffectInfo,
			new ReciprocatedTwoPathFunction(pEffectInfo->interactionName1()));
	}
	else if (effectName == "to")
	{
		pEffect = new GenericNetworkEffect(pEffectInfo,
			new MixedTwoPathFunction(pEffectInfo->interactionName1(),
				pEffectInfo->variableName()));
	}
	else if (effectName == "covNetNet")
	{
		string networkName = pEffectInfo->interactionName1();
		string covariateName = pEffectInfo->interactionName2();
		pEffect = new GenericNetworkEffect(pEffectInfo,
			new ConditionalFunction(new EqualCovariatePredicate(covariateName),
				new InStarFunction(networkName),
				0),
			new ConditionalFunction(
				new MissingCovariatePredicate(covariateName),
				0,
				new ConditionalFunction(
					new EqualCovariatePredicate(covariateName),
					new InStarFunction(networkName),
					0)));
	}				
	else if (effectName == "jumpWWClosure")
	{
		string networkName = pEffectInfo->interactionName1();
		string covariateName = pEffectInfo->interactionName2();
		pEffect = new GenericNetworkEffect(pEffectInfo,
			new ConditionalFunction(new EqualCovariatePredicate(covariateName),
				0,
				new SameCovariateTwoPathFunction(networkName, 
										covariateName, false)),
			new ConditionalFunction(
				new MissingCovariatePredicate(covariateName),
				0,
				new ConditionalFunction(
					new EqualCovariatePredicate(covariateName),
					0,
					new SameCovariateTwoPathFunction(networkName, 
										covariateName, true))));
	}
	else if (effectName == "jumpWXClosure")
	{
		string networkName = pEffectInfo->interactionName1();
		string covariateName = pEffectInfo->interactionName2();
		pEffect = new GenericNetworkEffect(pEffectInfo,
			new ConditionalFunction(new EqualCovariatePredicate(covariateName),
				0,
				new SameCovariateMixedTwoPathFunction(
								pEffectInfo->variableName(),
								networkName, covariateName, false)),
			new ConditionalFunction(
				new MissingCovariatePredicate(covariateName),
				0,
				new ConditionalFunction(
					new EqualCovariatePredicate(covariateName),
					0,
					new SameCovariateMixedTwoPathFunction(
							pEffectInfo->variableName(), 
							networkName, covariateName, true))));
	}
	else if (effectName == "closure")
	{
		pEffect = new GenericNetworkEffect(pEffectInfo,
			new TwoPathFunction(pEffectInfo->interactionName1()));
	}
	else if (effectName == "cyClosure")
	{
		pEffect = new GenericNetworkEffect(pEffectInfo,
			new ReverseTwoPathFunction(pEffectInfo->interactionName1()));
	}
	else if (effectName == "sharedIn")
	{
		pEffect = new GenericNetworkEffect(pEffectInfo,
			new OutStarFunction(pEffectInfo->interactionName1()));
	}
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
		pEffect = new SimilarityEffect(pEffectInfo, true, false, false);
	}
	else if (effectName == "totSim")
	{
		pEffect = new SimilarityEffect(pEffectInfo, false, false, false);
	}
	else if (effectName == "indeg")
	{
		pEffect = new IndegreeEffect(pEffectInfo);
	}
	else if (effectName == "outdeg")
	{
		pEffect = new OutdegreeEffect(pEffectInfo);
	}
	else if (effectName == "isolate")
	{
		pEffect = new IsolateEffect(pEffectInfo);
	}
	else if (effectName == "isolateNet")
	{
		pEffect = new IsolateNetEffect(pEffectInfo);
	}
	else if (effectName == "isolatePop")
	{
		pEffect = new IsolatePopEffect(pEffectInfo);
	}
	else if (effectName == "inIsDegree")
	{
		pEffect = new InIsolateDegreeEffect(pEffectInfo);
	}
	else if (effectName == "avSimRecip")
	{
		pEffect = new ReciprocatedSimilarityEffect(pEffectInfo, true, false);
	}
	else if (effectName == "totSimRecip")
	{
		pEffect = new ReciprocatedSimilarityEffect(pEffectInfo, false, false);
	}
	else if (effectName == "avSimPopAlt")
	{
		pEffect = new SimilarityEffect(pEffectInfo, true, true, false);
	}
	else if (effectName == "totSimPopAlt")
	{
		pEffect = new SimilarityEffect(pEffectInfo, false, true, false);
	}
	else if (effectName == "popAlt")
	{
		pEffect = new PopularityAlterEffect(pEffectInfo);
	}
	else if (effectName == "avSimRecPop")
	{
		pEffect = new ReciprocatedSimilarityEffect(pEffectInfo, true, true);
	}
	else if (effectName == "totSimRecPop")
	{
		pEffect = new ReciprocatedSimilarityEffect(pEffectInfo, false, true);
	}
	else if (effectName == "avAlt")
	{
		pEffect = new AverageAlterEffect(pEffectInfo);
	}
	else if (effectName == "avRecAlt")
	{
		pEffect = new AverageReciprocatedAlterEffect(pEffectInfo);
	}
	else if (effectName == "behDenseTriads")
	{
		pEffect = new DenseTriadsBehaviorEffect(pEffectInfo);
	}
	else if (effectName == "simDenseTriads")
	{
		pEffect = new DenseTriadsSimilarityEffect(pEffectInfo);
	}
	else if (effectName == "recipDeg")
	{
		pEffect = new ReciprocalDegreeBehaviorEffect(pEffectInfo);
	}
	else if (effectName == "avSimPopEgo")
	{
		pEffect = new SimilarityEffect(pEffectInfo, true, false, true);
	}
	else if (effectName == "effFrom")
	{
		pEffect = new MainCovariateEffect(pEffectInfo);
	}
	// else if (effectName == "inflIntX")
	// {
	// 	pEffect = new InteractionCovariateEffect(pEffectInfo, false, false,
	// 		false);
	//}
	else if (effectName == "avSimEgoX")
	{
		pEffect = new InteractionCovariateEffect(pEffectInfo, true, false, false);
	}
	else if (effectName == "totSimEgoX")
	{
		pEffect = new InteractionCovariateEffect(pEffectInfo, false, true, false);
	}
	else if (effectName == "avAltEgoX")
	{
		pEffect = new InteractionCovariateEffect(pEffectInfo, false, false, true);
	}
	else if (effectName == "AltsAvAlt")
	{
		pEffect = new AltersCovariateAverageEffect(pEffectInfo);
	}
	else if (effectName == "altDist2")
	{
		//	pEffect = new CovariateDistance2AlterEffect(pEffectInfo);
		string networkName = pEffectInfo->variableName();
		string covariateName = pEffectInfo->interactionName1();
		double parameter = pEffectInfo->internalEffectParameter();
		AlterFunction * pChangeFunction =
			new CovariateDistance2AlterNetworkFunction(networkName,
				covariateName, parameter, false);
		AlterFunction * pStatisticFunction =
			new CovariateDistance2AlterNetworkFunction(networkName,
				covariateName, parameter, true);
		pEffect = new GenericNetworkEffect(pEffectInfo,
			pChangeFunction, pStatisticFunction);
	}
	else if (effectName == "simDist2")
	{
		//	pEffect = new CovariateDistance2SimilarityEffect(pEffectInfo);
		string networkName = pEffectInfo->variableName();
		string covariateName = pEffectInfo->interactionName1();
		AlterFunction * pChangeFunction =
			new CovariateDistance2SimilarityNetworkFunction(networkName,
				covariateName, false);
		AlterFunction * pStatisticFunction =
			new CovariateDistance2SimilarityNetworkFunction(networkName,
				covariateName, true);
		pEffect = new GenericNetworkEffect(pEffectInfo,
			pChangeFunction, pStatisticFunction);
	}
	else if (effectName == "altDist2W")
	{
		string networkName = pEffectInfo->interactionName1();
		string covariateName = pEffectInfo->interactionName2();
		double parameter = pEffectInfo->internalEffectParameter();
		AlterFunction * pChangeFunction =
			new CovariateDistance2AlterNetworkFunction(networkName,
				covariateName, parameter, false);
		AlterFunction * pStatisticFunction =
			new CovariateDistance2AlterNetworkFunction(networkName,
				covariateName, parameter, true);
		pEffect = new GenericNetworkEffect(pEffectInfo,
			pChangeFunction, pStatisticFunction);
	}
	else if (effectName == "simDist2W")
	{
		string networkName = pEffectInfo->interactionName1();
		string covariateName = pEffectInfo->interactionName2();
		AlterFunction * pChangeFunction =
			new CovariateDistance2SimilarityNetworkFunction(networkName,
				covariateName, false);
		AlterFunction * pStatisticFunction =
			new CovariateDistance2SimilarityNetworkFunction(networkName,
				covariateName, true);
		pEffect = new GenericNetworkEffect(pEffectInfo,
			pChangeFunction, pStatisticFunction);
	}
	else
	{
		throw domain_error("Unexpected effect name: " + effectName);
	}

	return pEffect;
}

}
