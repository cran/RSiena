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

#include <stdexcept>
#include <vector>
#include <cmath>
#include <R.h>

#include "StatisticCalculator.h"
#include "data/Data.h"
#include "data/DataUtils.h"
#include "data/Network.h"
#include "data/NetworkLongitudinalData.h"
#include "data/BehaviorLongitudinalData.h"
#include "data/OneModeNetwork.h"
#include "data/TieIterator.h"
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

	// set up the predictor states of these variables

	for (unsigned i = 0; i < rVariables.size(); i++)
	{
		NetworkLongitudinalData * pNetworkData =
			dynamic_cast<NetworkLongitudinalData *>(rVariables[i]);
		BehaviorLongitudinalData * pBehaviorData =
			dynamic_cast<BehaviorLongitudinalData *>(rVariables[i]);
		string name = rVariables[i]->name();
		if (pNetworkData)
		{
			Network * pPredictor =
				pNetworkData->pNetwork(this->lperiod)->clone();
			this->subtractNetwork(pPredictor,
				pNetworkData->pMissingTieNetwork(this->lperiod));
			this->subtractNetwork(pPredictor,
				pNetworkData->pMissingTieNetwork(this->lperiod + 1));
			NetworkVariable * pNetworkVariable = (NetworkVariable *)
				this->lpState->pSimulation()->pVariable(name);
			pNetworkVariable->pPredictorNetwork(pPredictor);
			pNetworkVariable->changingCovariatePeriod(this->lperiod);
		}
		else if (pBehaviorData)
		{
 			// create a copy of the start of the period and zero any values
			// missing at (either end?) start of period

			// and the corresponding variable
			BehaviorVariable * pBehaviorVariable = (BehaviorVariable *)
				this->lpState->pSimulation()->pVariable(name);

			int predictorValue;

			const int * Start = pBehaviorData->values(this->lperiod);

			for (int i = 0; i < pBehaviorData->n(); i++)
			{
				predictorValue = Start[i];
				//		if (pBehaviorData->missing(this->lperiod, i) ||
				//			pBehaviorData->missing(this->lperiod + 1, i))
				if (pBehaviorData->missing(this->lperiod, i) )
				{
					predictorValue = 0;
				}
				pBehaviorVariable->predictorValue(i, predictorValue);
			}
			pBehaviorVariable->changingCovariatePeriod(this->lperiod);
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
	// Duplicate the current network and remove those ties that are
	// missing at either end of the period.

	Network * pNetwork =
		this->lpState->pNetwork(pNetworkData->name())->clone();

	this->subtractNetwork(pNetwork,
		pNetworkData->pMissingTieNetwork(this->lperiod));
	this->subtractNetwork(pNetwork,
		pNetworkData->pMissingTieNetwork(this->lperiod + 1));

	// Loop through the evaluation effects, calculate the statistics,
	// and store them.

 	const vector<EffectInfo *> & rEffects =
 		this->lpModel->rEvaluationEffects(pNetworkData->name());

 	EffectFactory factory;

	for (unsigned i = 0; i < rEffects.size(); i++)
	{
		EffectInfo * pInfo = rEffects[i];
		NetworkEffect * pEffect =
			(NetworkEffect *) factory.createEffect(pInfo);

		// Initialize the effect to work with our data and state of variables.
		pEffect->initialize(this->lpData, this->lpState, this->lperiod);

		this->lstatistics[pInfo] =
			pEffect->evaluationStatistic(pNetwork);
		delete pEffect;
	}

	delete pNetwork;
}


/**
 * Calculates the statistics for the endowment effects of the given
 * network variable.
 */
void StatisticCalculator::calculateNetworkEndowmentStatistics(
	NetworkLongitudinalData * pNetworkData)
{
	// In order to calculate the statistics of network endowment effects,
	// we need the initial network of the given period, and the network
	// of lost ties, namely, the ties that are present in the initial
	// network, not missing at the end of the period, and absent in the
	// current network.

	Network * pInitialNetwork = pNetworkData->pNetwork(this->lperiod);
	Network * pLostTieNetwork = pInitialNetwork->clone();

	this->subtractNetwork(pLostTieNetwork,
		this->lpState->pNetwork(pNetworkData->name()));
	this->subtractNetwork(pLostTieNetwork,
		pNetworkData->pMissingTieNetwork(this->lperiod + 1));

	// Loop through the endowment effects, calculate the statistics,
	// and store them.

	const vector<EffectInfo *> & rEffects =
		this->lpModel->rEndowmentEffects(pNetworkData->name());

 	EffectFactory factory;

	for (unsigned i = 0; i < rEffects.size(); i++)
	{
		EffectInfo * pInfo = rEffects[i];
		NetworkEffect * pEffect =
			(NetworkEffect *) factory.createEffect(pInfo);

		// Initialize the effect to work with our data and state of variables.
		pEffect->initialize(this->lpData, this->lpState, this->lperiod);

		this->lstatistics[pInfo] =
			pEffect->endowmentStatistic(pInitialNetwork, pLostTieNetwork);
		delete pEffect;
	}

	delete pLostTieNetwork;
}


/**
 * Removes the ties of the first network that are present in the second
 * network.
 */
void StatisticCalculator::subtractNetwork(Network * pNetwork,
	const Network * pSubtrahendNetwork) const
{
	for (TieIterator iter = pSubtrahendNetwork->ties();
		iter.valid();
		iter.next())
	{
		pNetwork->setTieValue(iter.ego(), iter.alter(), 0);
	}
}


/**
 * Calculates the statistics for effects of the given behavior variable.
 */
void StatisticCalculator::calculateBehaviorStatistics(
	BehaviorLongitudinalData * pBehaviorData)
{
	// TODO

	// create a copy of the current state and zero any values missing
	//at either end of period
	const int * currentState = this->lpState->
		behaviorValues(pBehaviorData->name());

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
		difference[i] = pBehaviorData->value(this->lperiod, i) -
			currentState[i];
		if (pBehaviorData->missing(this->lperiod, i) ||
			pBehaviorData->missing(this->lperiod + 1, i))
		{
			difference[i] = 0;
		}
	}

	// Loop through the evaluation effects, calculate the statistics,
	// and store them.

	const vector<EffectInfo *> & rEffects =
		this->lpModel->rEvaluationEffects(pBehaviorData->name());

	const BehaviorVariable * pVariable = (const BehaviorVariable *)
			this->lpState->pSimulation()->pVariable(pBehaviorData->name());

	for (unsigned i = 0; i < rEffects.size(); i++)
	{
		EffectInfo * pInfo = rEffects[i];

		const BehaviorEffect * pEffect  = (const BehaviorEffect *)
			pVariable->pEvaluationFunction()->rEffects()[i];

		this->lstatistics[pInfo] =
			pEffect->evaluationStatistic(currentValues);
	}

	// Loop through the endowment effects, calculate the statistics,
	// and store them.

	const vector<EffectInfo *> & rEndowEffects =
		this->lpModel->rEndowmentEffects(pBehaviorData->name());
	for (unsigned i = 0; i < rEndowEffects.size(); i++)
	{
		EffectInfo * pInfo = rEndowEffects[i];

		const BehaviorEffect * pEffect  = (const BehaviorEffect *)
			pVariable->pEndowmentFunction()->rEffects()[i];

		this->lstatistics[pInfo] =
			pEffect->endowmentStatistic(difference, currentValues);
	}
	delete[] currentValues;
}


/**
 * Calculates the statistics for the rate effects of the given
 * network variable.
 */
void StatisticCalculator::calculateNetworkRateStatistics(
	NetworkLongitudinalData * pNetworkData)
{
	// Duplicate the current network and remove those ties that are
	// missing at either end of the period. TODO set leavers back.

	Network * pNetwork =
		this->lpState->pNetwork(pNetworkData->name())->clone();

	this->subtractNetwork(pNetwork,
		pNetworkData->pMissingTieNetwork(this->lperiod));
	this->subtractNetwork(pNetwork,
		pNetworkData->pMissingTieNetwork(this->lperiod + 1));

	// construct a network of differences between current and start
    // of period.

	Network * pStart = pNetworkData->pNetwork(this->lperiod)->clone();

	this->subtractNetwork(pStart,
		pNetworkData->pMissingTieNetwork(this->lperiod));
	this->subtractNetwork(pStart,
		pNetworkData->pMissingTieNetwork(this->lperiod + 1));

	Network * pDifference = symmetricDifference(pStart, pNetwork);

	// basic rate distance

	if (!this->ldistances[pNetworkData])
	{
		int * array =
			new int[pNetworkData->observationCount() - 1];

		this->ldistances[pNetworkData] = array;
	}

	this->ldistances[pNetworkData][this->lperiod] = pDifference->tieCount();

//	Rprintf("basic rate change %d\n", pDifference->tieCount());

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

			Network * pStructural =
				pExplanatoryNetwork->pNetwork(this->lperiod)->clone();
			this->subtractNetwork(pStructural,
				pExplanatoryNetwork->pMissingTieNetwork(this->lperiod));

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
			delete pStructural;
		}
	}

	delete pStart;
	delete pDifference;
	delete pNetwork;
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
			Network * pStructural =
				pNetworkData->pNetwork(this->lperiod)->clone();
			this->subtractNetwork(pStructural,
				pNetworkData->pMissingTieNetwork(this->lperiod));

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
			delete pStructural;
		}
	}

	delete[] difference;
	delete[] currentValues;
}

}
