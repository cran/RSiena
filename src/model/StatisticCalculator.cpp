/******************************************************************************
 * SIENA: Simulation Investigation for Empirical Network Analysis
 *
 * Web: http://www.stats.ox.ac.uk/~snijders/siena/
 *
 * File: StatisticCalculator.cpp
 *
 * Description: This file contains the implementation of the
 * StatisticCalculator class.
 *****************************************************************************/

#include <cstdlib>
#include <stdexcept>
#include <vector>
#include <cmath>

#include "StatisticCalculator.h"
#include "data/Data.h"
#include "network/NetworkUtils.h"
#include "network/Network.h"
#include "data/NetworkLongitudinalData.h"
#include "data/BehaviorLongitudinalData.h"
#include "network/OneModeNetwork.h"
#include "network/TieIterator.h"
#include "data/ConstantCovariate.h"
#include "data/ChangingCovariate.h"
#include "model/Model.h"
#include "model/State.h"
#include "model/Function.h"
#include "model/EffectInfo.h"
#include "model/effects/EffectFactory.h"
#include "model/effects/NetworkEffect.h"
#include "model/effects/BehaviorEffect.h"
#include "model/EpochSimulation.h"
#include "model/variables/NetworkVariable.h"
#include "model/variables/BehaviorVariable.h"
#include "model/tables/Cache.h"

namespace siena
{

// ----------------------------------------------------------------------------
// Section: Construction, destruction
// ----------------------------------------------------------------------------

/**
 * Constructor.
 * @param[in] pData the observed data
 * @param[in] pModel the model whose effect statistics are to be calculated
 * @param[in] pState the current state of the dependent variables
 * @param[in] period the period under consideration
 */
StatisticCalculator::StatisticCalculator(const Data * pData,
	const Model * pModel,
	State * pState,
	int period)
{
	this->lpData = pData;
	this->lpModel = pModel;
	this->lpState = pState;
	this->lperiod = period;
	this->lpPredictorState = new State();
	this->lpStateLessMissingsEtc = new State();

	this->calculateStatistics();
}


/**
 * Deallocates this model.
 */
StatisticCalculator::~StatisticCalculator()
{
	// Delete the arrays of simulated distances

	while (!this->ldistances.empty())
	{
		int * array = this->ldistances.begin()->second;
		this->ldistances.erase(this->ldistances.begin());
		delete[] array;
	}

	// The state just stores the values, but does not own them. It means that
	// the destructor of the state won't deallocate the memory, and we must
	// request that explicitly. The basic predictor state is now
	// stored on the data so lpredictorState is just a pointer.

	//	this->lpPredictorState->deleteValues();
	delete this->lpPredictorState;
	this->lpPredictorState = 0;
	this->lpStateLessMissingsEtc->deleteValues();
	delete this->lpStateLessMissingsEtc;
	this->lpStateLessMissingsEtc = 0;
}


// ----------------------------------------------------------------------------
// Section: Accessors for statistics
// ----------------------------------------------------------------------------

/**
 * Returns the statistic for the given effect.
 */
double StatisticCalculator::statistic(EffectInfo * pEffect) const
{
	map<EffectInfo *, double>::const_iterator iter =
		this->lstatistics.find(pEffect);

	if (iter == this->lstatistics.end())
	{
		throw invalid_argument(
			"Unknown effect: The given effect is not part of the model.");
	}

	return iter->second;
}


/**
 * Returns the simulated distance for the given network and period.
 */
int StatisticCalculator::distance(LongitudinalData * pData, int period)
	const
{
	map<LongitudinalData *, int *>::const_iterator iter =
		this->ldistances.find(pData);

	if (iter == this->ldistances.end())
	{
		throw invalid_argument(
			"Unknown effect: The given basic rate is not part of the model.");
	}

	return iter->second[period];
}


/**
 * Calculates the statistics for all effects of the given model. Note that
 * this->lperiod relates to the current period when simulating, but the
 * previous when calculating targets.
 */
void StatisticCalculator::calculateStatistics()
{
	const vector<LongitudinalData *> & rVariables =
		this->lpData->rDependentVariableData();

	// set up the predictor and currentLessMissingsEtc states of these variables

	for (unsigned i = 0; i < rVariables.size(); i++)
	{
		NetworkLongitudinalData * pNetworkData =
			dynamic_cast<NetworkLongitudinalData *>(rVariables[i]);
		BehaviorLongitudinalData * pBehaviorData =
			dynamic_cast<BehaviorLongitudinalData *>(rVariables[i]);
		string name = rVariables[i]->name();

		if (pNetworkData)
		{
			const Network * pPredictor =
				pNetworkData->pNetworkLessMissing(this->lperiod);
			this->lpPredictorState->pNetwork(name, pPredictor);

			// Duplicate the current network and remove those ties that are
			// missing at either end of the period.

			Network * pNetwork =
				this->lpState->pNetwork(pNetworkData->name())->clone();

			subtractNetwork(pNetwork,
				pNetworkData->pMissingTieNetwork(this->lperiod));
			subtractNetwork(pNetwork,
				pNetworkData->pMissingTieNetwork(this->lperiod + 1));

			// for not-targets, overwrite the current network for values
			// structurally fixed for the next period. (no effect for targets)

			replaceNetwork(pNetwork,
				pNetworkData->pNetwork(this->lperiod + 1),
				pNetworkData->pStructuralTieNetwork(this->lperiod + 1));

			// for targets look backwards and mimic the simulation by carrying
			// forward structural values.

			replaceNetwork(pNetwork,
				pNetworkData->pNetwork(this->lperiod),
				pNetworkData->pStructuralTieNetwork(this->lperiod));

			this->lpStateLessMissingsEtc->pNetwork(name, pNetwork);

		}
		else if (pBehaviorData)
		{
 			// create a copy of the start of the period and zero any values
			// missing at (either end?) start of period

			const int * values =
				pBehaviorData->valuesLessMissingStarts(this->lperiod);
			this->lpPredictorState->behaviorValues(name, values);
		}
		else
		{
			throw domain_error("Unexpected class of dependent variable");
		}
	}

	for (unsigned i = 0; i < rVariables.size(); i++)
	{
		NetworkLongitudinalData * pNetworkData =
			dynamic_cast<NetworkLongitudinalData *>(rVariables[i]);
		BehaviorLongitudinalData * pBehaviorData =
			dynamic_cast<BehaviorLongitudinalData *>(rVariables[i]);

		if (pNetworkData)
		{
			this->calculateNetworkRateStatistics(pNetworkData);
			this->calculateNetworkEvaluationStatistics(pNetworkData);
			this->calculateNetworkEndowmentStatistics(pNetworkData);
			this->calculateNetworkCreationStatistics(pNetworkData);
		}
		else if (pBehaviorData)
		{
			this->calculateBehaviorRateStatistics(pBehaviorData);
			this->calculateBehaviorStatistics(pBehaviorData);
		}
		else
		{
			throw domain_error("Unexpected class of dependent variable");
		}
	}
}


/**
 * Calculates the statistics for the evaluation effects of the given
 * network variable.
 */
void StatisticCalculator::calculateNetworkEvaluationStatistics(
	NetworkLongitudinalData * pNetworkData)
{
	// We want to pass all networks to the effects in a single state,
	// hence we overwrite the network in the predictor state.
	// We do not use the predictor network with effects this network owns.

	string name = pNetworkData->name();
	const Network * pPredictorNetwork = this->lpPredictorState->pNetwork(name);
	const Network * pCurrentLessMissingsEtc =
		this->lpStateLessMissingsEtc->pNetwork(name);
	this->lpPredictorState->pNetwork(name, pCurrentLessMissingsEtc);

	// Loop through the evaluation effects, calculate the statistics,
	// and store them.

 	const vector<EffectInfo *> & rEffects =
 		this->lpModel->rEvaluationEffects(pNetworkData->name());

 	EffectFactory factory(this->lpData);
 	Cache cache;

	for (unsigned i = 0; i < rEffects.size(); i++)
	{
		EffectInfo * pInfo = rEffects[i];
		NetworkEffect * pEffect =
			(NetworkEffect *) factory.createEffect(pInfo);

		// Initialize the effect to work with our data and state of variables.

		pEffect->initialize(this->lpData,
			this->lpPredictorState,
			this->lperiod,
			&cache);

		this->lstatistics[pInfo] = pEffect->evaluationStatistic();
		delete pEffect;
	}

	// Restore the predictor network
	this->lpPredictorState->pNetwork(name, pPredictorNetwork);

}


/**
 * Calculates the statistics for the endowment effects of the given
 * network variable.
 */
void StatisticCalculator::calculateNetworkEndowmentStatistics(
	NetworkLongitudinalData * pNetworkData)
{
	// To save a lot of unnecessary effort, check first we have some
	// endowment effects
	const vector<EffectInfo *> & rEffects =
		this->lpModel->rEndowmentEffects(pNetworkData->name());

	if (rEffects.size() > 0)
	{
		// In order to calculate the statistics of network endowment effects,
		// we need the initial network of the given period, and the network
		// of lost ties, namely, the ties that are present in the initial
		// network, not missing at the start of the period, and absent in the
		// current network.

		const Network * pInitialNetwork = pNetworkData->pNetwork(this->lperiod);
		Network * pLostTieNetwork = pInitialNetwork->clone();

		// Duplicate the current network so can overwrite the structurals
		Network * pCurrentNetwork =
			this->lpState->pNetwork(pNetworkData->name())->clone();

		// Replace values for structurally determined values from previous
		// or current period
		replaceNetwork(pCurrentNetwork,
			pNetworkData->pNetwork(this->lperiod + 1),
			pNetworkData->pStructuralTieNetwork(this->lperiod + 1));

		replaceNetwork(pCurrentNetwork,
			pNetworkData->pNetwork(this->lperiod),
			pNetworkData->pStructuralTieNetwork(this->lperiod));

		// remove missings and current

		subtractNetwork(pLostTieNetwork,
			pCurrentNetwork);
		subtractNetwork(pLostTieNetwork,
			pNetworkData->pMissingTieNetwork(this->lperiod + 1));

		// overwrite the predictor network with only start missings removed
		const Network * pPredictor =
			pNetworkData->pNetworkLessMissingStart(this->lperiod);

		// We want to pass all networks to the effects in a single state,
		// hence we overwrite the network in the predictor state.

		string name = pNetworkData->name();
		const Network * pPredictorNetwork =
			this->lpPredictorState->pNetwork(name);
		this->lpPredictorState->pNetwork(name, pPredictor);

		// Loop through the endowment effects, calculate the statistics,
		// and store them.

		EffectFactory factory(this->lpData);
		Cache cache;

		for (unsigned i = 0; i < rEffects.size(); i++)
		{
			EffectInfo * pInfo = rEffects[i];
			NetworkEffect * pEffect =
				(NetworkEffect *) factory.createEffect(pInfo);

			// Initialize the effect to work with our data and state of variables.

			pEffect->initialize(this->lpData,
				this->lpPredictorState,
				this->lperiod,
				&cache);

			this->lstatistics[pInfo] =
				pEffect->endowmentStatistic(pLostTieNetwork);
			delete pEffect;
		}

		// Restore the predictor network
		this->lpPredictorState->pNetwork(name, pPredictorNetwork);

		delete pCurrentNetwork;
		delete pLostTieNetwork;
	}
}


/**
 * Calculates the statistics for the creation effects of the given
 * network variable.
 */
void StatisticCalculator::calculateNetworkCreationStatistics(
	NetworkLongitudinalData * pNetworkData)
{
	// To save a lot of unnecessary effort, check first we have some
	// endowment effects
	const vector<EffectInfo *> & rEffects =
		this->lpModel->rCreationEffects(pNetworkData->name());

	if (rEffects.size() > 0)
	{
		// We want to pass all networks to the effects in a single state,
		// hence we overwrite the network in the predictor state.
		// We do not use the predictor network with effects this network owns.

		string name = pNetworkData->name();
		const Network * pPredictorNetwork =
			this->lpPredictorState->pNetwork(name);
		const Network * pCurrentLessMissingsEtc =
			this->lpStateLessMissingsEtc->pNetwork(name);
		this->lpPredictorState->pNetwork(name, pCurrentLessMissingsEtc);

		// In order to calculate the statistics of tie creation effects, we
		// need to know the ties that are gained, namely, those that are
		// known to be absent in the beginning of the period.

		Network * pGainedTieNetwork = pCurrentLessMissingsEtc->clone();

		subtractNetwork(pGainedTieNetwork,
			pNetworkData->pNetwork(this->lperiod));
		// not sure we need this one!
		subtractNetwork(pGainedTieNetwork,
			pNetworkData->pMissingTieNetwork(this->lperiod));

		// Loop through the creation effects, calculate the statistics,
		// and store them.

		EffectFactory factory(this->lpData);
		Cache cache;

		for (unsigned i = 0; i < rEffects.size(); i++)
		{
			EffectInfo * pInfo = rEffects[i];
			NetworkEffect * pEffect =
				(NetworkEffect *) factory.createEffect(pInfo);

			// Initialize the effect to work with our data and state of variables.

			pEffect->initialize(this->lpData,
				this->lpPredictorState,
				this->lperiod,
				&cache);

			this->lstatistics[pInfo] =
				pEffect->creationStatistic(pGainedTieNetwork);

			delete pEffect;
		}

		// Restore the predictor network
		this->lpPredictorState->pNetwork(name, pPredictorNetwork);

		delete pGainedTieNetwork;
	}
}

/**
 * Calculates the statistics for effects of the given behavior variable.
 */
void StatisticCalculator::calculateBehaviorStatistics(
	BehaviorLongitudinalData * pBehaviorData)
{
	// create a copy of the current state and zero any values missing
	// at either end of period

	const int * currentState =
		this->lpState->behaviorValues(pBehaviorData->name());

	double * currentValues  = new double[pBehaviorData->n()];

	for (int i = 0; i < pBehaviorData->n(); i++)
	{
		currentValues[i] = currentState[i] - pBehaviorData->overallMean();

		if (pBehaviorData->missing(this->lperiod, i) ||
			pBehaviorData->missing(this->lperiod + 1, i))
		{
			currentValues[i] = 0;
		}
	}

	// Construct a vector of difference values relative to previous period.
    // Values for missing values are set to 0.

	int * difference = new int[pBehaviorData->n()];

	for (int i = 0; i < pBehaviorData->n(); i++)
	{
		difference[i] =
			pBehaviorData->value(this->lperiod, i) - currentState[i];

		if (pBehaviorData->missing(this->lperiod, i) ||
			pBehaviorData->missing(this->lperiod + 1, i))
		{
			difference[i] = 0;
		}
	}

	// Loop through the evaluation effects, calculate the statistics,
	// and store them.

	const vector<EffectInfo *> & rEvaluationEffects =
		this->lpModel->rEvaluationEffects(pBehaviorData->name());

 	EffectFactory factory(this->lpData);
 	Cache cache;

	for (unsigned i = 0; i < rEvaluationEffects.size(); i++)
	{
		EffectInfo * pInfo = rEvaluationEffects[i];
		BehaviorEffect * pEffect =
			(BehaviorEffect *) factory.createEffect(pInfo);

		// Initialize the effect to work with our data and state of variables.

		pEffect->initialize(this->lpData,
			this->lpPredictorState,
			this->lperiod,
			&cache);

		this->lstatistics[pInfo] =
			pEffect->evaluationStatistic(currentValues);

		delete pEffect;
	}

	// Loop through the endowment effects, calculate the statistics,
	// and store them.

	const vector<EffectInfo *> & rEndowmentEffects =
		this->lpModel->rEndowmentEffects(pBehaviorData->name());

	for (unsigned i = 0; i < rEndowmentEffects.size(); i++)
	{
		EffectInfo * pInfo = rEndowmentEffects[i];
		BehaviorEffect * pEffect =
			(BehaviorEffect *) factory.createEffect(pInfo);

		// Initialize the effect to work with our data and state of variables.

		pEffect->initialize(this->lpData,
			this->lpPredictorState,
			this->lperiod,
			&cache);

		this->lstatistics[pInfo] =
			pEffect->endowmentStatistic(difference, currentValues);

		delete pEffect;
	}

	// Loop through the creation effects, calculate the statistics,
	// and store them.

	const vector<EffectInfo *> & rCreationEffects =
		this->lpModel->rCreationEffects(pBehaviorData->name());

	for (unsigned i = 0; i < rCreationEffects.size(); i++)
	{
		EffectInfo * pInfo = rCreationEffects[i];
		BehaviorEffect * pEffect =
			(BehaviorEffect *) factory.createEffect(pInfo);

		// Initialize the effect to work with our data and state of variables.

		pEffect->initialize(this->lpData,
			this->lpPredictorState,
			this->lperiod,
			&cache);

		this->lstatistics[pInfo] =
			pEffect->creationStatistic(difference, currentValues);

		delete pEffect;
	}

	delete[] currentValues;
	delete[] difference;
}


/**
 * Calculates the statistics for the rate effects of the given
 * network variable.
 */
void StatisticCalculator::calculateNetworkRateStatistics(
	NetworkLongitudinalData * pNetworkData)
{
	Network * pNetwork = 0;
	const Network * pConstNetwork =	this->lpStateLessMissingsEtc->
			pNetwork(pNetworkData->name());
	Network * pDifference = 0;

	// if parallel running, do processing from scratch
	if (this->lpModel->parallelRun())
	{
		// Duplicate the current network and remove those ties that are
		// missing at either end of the period. TODO set leavers back.
		// (Is the TODO not done for the current network?)

		pNetwork =
			this->lpState->pNetwork(pNetworkData->name())->clone();
		subtractNetwork(pNetwork,
			pNetworkData->pMissingTieNetwork(this->lperiod));
		subtractNetwork(pNetwork,
			pNetworkData->pMissingTieNetwork(this->lperiod + 1));

		// Replace values for structurally determined values from previous
		// or current period. Siena 3 does not do this ... so we won't yet

		// 	this->replaceNetwork(pNetwork,
		// 		pNetworkData->pNetwork(this->lperiod + 1),
		// 		pNetworkData->pStructuralTieNetwork(this->lperiod + 1));

		// 	this->replaceNetwork(pNetwork,
		// 		pNetworkData->pNetwork(this->lperiod),
		// 		pNetworkData->pStructuralTieNetwork(this->lperiod));

		// instead, remove all structurally determined ties at either end
		subtractNetwork(pNetwork,
			pNetworkData->pStructuralTieNetwork(this->lperiod + 1));

		subtractNetwork(pNetwork,
			pNetworkData->pStructuralTieNetwork(this->lperiod));
	}

	// construct a network of differences between current and start
	// of period.

	Network * pStart = pNetworkData->
		pNetworkLessMissing(this->lperiod)->clone();

	// must match with above code: the net result is probably equivalent
	if (this->lpModel->parallelRun())
	{
		// remove all structurally determined ties at either end
		subtractNetwork(pStart,
			pNetworkData->pStructuralTieNetwork(this->lperiod + 1));

		subtractNetwork(pStart,
			pNetworkData->pStructuralTieNetwork(this->lperiod));

		pDifference = symmetricDifference(pStart, pNetwork);
	}
	else
	{
		replaceNetwork(pStart,
			pNetworkData->pNetwork(this->lperiod + 1),
			pNetworkData->pStructuralTieNetwork(this->lperiod + 1));

		replaceNetwork(pStart,
			pNetworkData->pNetwork(this->lperiod),
			pNetworkData->pStructuralTieNetwork(this->lperiod));

		pDifference = symmetricDifference(pStart, pConstNetwork);
	}

	// basic rate distance

	if (!this->ldistances[pNetworkData])
	{
		int * array =
			new int[pNetworkData->observationCount() - 1];

		this->ldistances[pNetworkData] = array;
	}

	this->ldistances[pNetworkData][this->lperiod] = pDifference->tieCount();

	//Rprintf("basic rate change %d\n", pDifference->tieCount());

	// Loop through the rate effects, calculate the statistics,
	// and store them.
	const vector<EffectInfo *> & rEffects =
		this->lpModel->rRateEffects(pNetworkData->name());

	for (unsigned i = 0; i < rEffects.size(); i++)
	{
		EffectInfo * pInfo = rEffects[i];
		//	double parameter = pInfo->parameter();
		string effectName = pInfo->effectName();
		string interactionName = pInfo->interactionName1();
		string rateType = pInfo->rateType();

		if (rateType == "covariate")
		{
			// Covariate-dependent rate effect

			//	if (parameter != 0)
			//	{
			ConstantCovariate * pConstantCovariate =
				this->lpData->pConstantCovariate(interactionName);
			ChangingCovariate * pChangingCovariate =
				this->lpData->pChangingCovariate(interactionName);
			BehaviorLongitudinalData * pBehavior =
				this->lpData->pBehaviorData(interactionName);

			if (pConstantCovariate)
			{
				double statistic = 0;

				for (TieIterator iter = pDifference->ties();
					 iter.valid();
					 iter.next())
				{
					statistic += pConstantCovariate->value(iter.ego()) *
						iter.value();
				}

				this->lstatistics[pInfo] = statistic;
			}
			else if (pChangingCovariate)
			{
				double statistic = 0;

				for (TieIterator iter = pDifference->ties();
					 iter.valid();
					 iter.next())
				{
					statistic +=
						pChangingCovariate->value(iter.ego(),
							this->lperiod) *
						iter.value();
				}

				this->lstatistics[pInfo] = statistic;
			}
			else if (pBehavior)
			{
				double statistic = 0;

				for (TieIterator iter = pDifference->ties();
					 iter.valid();
					 iter.next())
				{
					statistic +=
						pBehavior->value(this->lperiod, iter.ego()) *
						iter.value();
				}

				this->lstatistics[pInfo] = statistic;
			}
			else
			{
				throw logic_error(
					"No individual covariate named '" +
					interactionName +
					"'.");
			}
			//}
		}
		else
		{
			// We expect a structural (network-dependent) rate effect here.
			NetworkLongitudinalData * pExplanatoryNetwork;
			if (interactionName == "")
			{
				pExplanatoryNetwork = pNetworkData;
			}
			else
			{
				pExplanatoryNetwork =
					this->lpData->pNetworkData(interactionName);
			}

			const Network * pStructural =
				pExplanatoryNetwork->pNetworkLessMissing(this->lperiod);

			double statistic = 0;

			for (TieIterator iter = pDifference->ties();
				 iter.valid();
				 iter.next())
			{
				if (effectName == "outRate")
				{
					statistic += pStructural->outDegree(iter.ego()) *
						iter.value();
				}
				else if (effectName == "inRate")
				{
					statistic += pStructural->inDegree(iter.ego()) *
						iter.value();
				}
				else if (effectName == "recipRate")
				{
					OneModeNetwork * pOneModeNetwork =
						(OneModeNetwork *) pStructural;
					statistic +=
						pOneModeNetwork->reciprocalDegree(iter.ego()) *
						iter.value();
				}
				else if (effectName == "outRateInv")
				{
					statistic +=
						1.0 / (pStructural->outDegree(iter.ego()) + 1) *
						iter.value();
				}
				else
				{
					throw domain_error("Unexpected rate effect " + effectName);
				}
			}

			this->lstatistics[pInfo] = statistic;
		}
	}

	delete pStart;
	delete pDifference;
	if (this->lpModel->parallelRun())
	{
		delete pNetwork;
	}

}
/**
 * Calculates the statistics for the rate effects of the given
 * network variable.
 */
void StatisticCalculator::calculateBehaviorRateStatistics(
	BehaviorLongitudinalData * pBehaviorData)
{
	// create a copy of the current state and zero any values missing
	//at either end of period
	const int * currentState = this->lpState->
		behaviorValues(pBehaviorData->name());

	int * currentValues  = new int[pBehaviorData->n()];

	for (int i = 0; i < pBehaviorData->n(); i++)
	{
		currentValues[i] = currentState[i];

		if (pBehaviorData->missing(this->lperiod, i) ||
			pBehaviorData->missing(this->lperiod + 1, i))
		{
			currentValues[i] = 0;
		}
	}
	// Construct a vector of absolute differences between current and start
    // of period. Differences for missing values are set to 0.

	const int * Start = pBehaviorData->values(this->lperiod);

	int * difference  = new int[pBehaviorData->n()];

	for (int i = 0; i < pBehaviorData->n(); i++)
	{
		difference[i] = abs(currentState[i] - Start[i]);
		if (pBehaviorData->missing(this->lperiod, i) ||
			pBehaviorData->missing(this->lperiod + 1, i))
		{
			difference[i] = 0;
		}
	}

	// basic rate distance

	if (!this->ldistances[pBehaviorData])
	{
		int * array =
			new int[pBehaviorData->observationCount() - 1];

		this->ldistances[pBehaviorData] = array;
	}

	int distance = 0;
	for (int i = 0; i < pBehaviorData->n(); i++)
	{
		distance += difference[i];
	}
	this->ldistances[pBehaviorData][this->lperiod] = distance;

	// Loop through the rate effects, calculate the statistics,
	// and store them.

	const vector<EffectInfo *> & rEffects =
		this->lpModel->rRateEffects(pBehaviorData->name());

	for (unsigned i = 0; i < rEffects.size(); i++)
	{
		EffectInfo * pInfo = rEffects[i];
		//	double parameter = pInfo->parameter();
		string effectName = pInfo->effectName();
		string interactionName = pInfo->interactionName1();
		string rateType = pInfo->rateType();

		if (rateType == "covariate")
		{
			// Covariate-dependent rate effect

			//	if (parameter != 0)
			//	{
				ConstantCovariate * pConstantCovariate =
					this->lpData->pConstantCovariate(interactionName);
				ChangingCovariate * pChangingCovariate =
					this->lpData->pChangingCovariate(interactionName);
				BehaviorLongitudinalData * pBehavior =
					this->lpData->pBehaviorData(interactionName);
				if (pConstantCovariate)
				{
					double statistic = 0;

					for (int i = 0; i < pBehaviorData->n(); i++)
					{
						statistic += pConstantCovariate->value(i) *
							difference[i];
					}

					this->lstatistics[pInfo] = statistic;
				}
				else if (pChangingCovariate)
				{
					double statistic = 0;

					for (int i = 0; i < pBehaviorData->n(); i++)
					{
						statistic +=
							pChangingCovariate->value(i,
								this->lperiod) *
							difference[i];
					}

					this->lstatistics[pInfo] = statistic;
				}
				else if (pBehavior)
				{
					double statistic = 0;

					for (int i = 0; i < pBehaviorData->n(); i++)
					{
						statistic +=
							pBehavior->values(this->lperiod)[i] *
							difference[i];
					}

					this->lstatistics[pInfo] = statistic;
				}
				else
				{
					throw logic_error(
						"No individual covariate named '" +
						interactionName +
						"'.");
				}
				//}
		}
		else
		{
			// We expect a structural (network-dependent) rate effect here.

			NetworkLongitudinalData *pNetworkData = this->lpData->
				pNetworkData(interactionName);
			const Network * pStructural =
				pNetworkData->pNetworkLessMissingStart(this->lperiod);

			double statistic = 0;
			for (int i = 0; i < pBehaviorData->n(); i++)
			{
				if (effectName == "outRate")
				{
					statistic += pStructural->outDegree(i) *
						difference[i];
				}
				else if (effectName == "inRate")
				{
					statistic += pStructural->inDegree(i) *
						difference[i];
				}
				else if (effectName == "recipRate")
				{
					OneModeNetwork * pOneModeNetwork =
						(OneModeNetwork *) pStructural;
					statistic +=
						pOneModeNetwork->reciprocalDegree(i) *
						difference[i];
				}
				else if (effectName == "outRateInv")
				{
					statistic +=
						1.0 / (pStructural->outDegree(i) + 1) *
						difference[i];
				}
				else
				{
					throw domain_error("Unexpected rate effect " + effectName);
				}
			}

			this->lstatistics[pInfo] = statistic;
		}
	}

	delete[] difference;
	delete[] currentValues;
}

}
