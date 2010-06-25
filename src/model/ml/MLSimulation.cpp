/******************************************************************************
 * SIENA: Simulation Investigation for Empirical Network Analysis
 *
 * Web: http://www.stats.ox.ac.uk/~snijders/siena/
 *
 * File: MLSimulation.cpp
 *
 * Description: This file contains the implementation of the class
 * MLSimulation.
 *****************************************************************************/
#include <stdexcept>
#include <string>
#include <cmath>
#include <Rinternals.h>
#include "MLSimulation.h"
#include "utils/Random.h"
#include "utils/Utils.h"
#include "data/NetworkLongitudinalData.h"
#include "data/BehaviorLongitudinalData.h"
#include "model/Model.h"
#include "model/SimulationActorSet.h"
#include "model/variables/DependentVariable.h"
#include "model/variables/BehaviorVariable.h"
#include "model/ml/Chain.h"
#include "model/ml/MiniStep.h"
#include "model/ml/NetworkChange.h"
#include "model/ml/BehaviorChange.h"
#include "model/tables/Cache.h"


namespace siena
{
SEXP getMiniStepDF(const MiniStep& miniStep);
SEXP getChainDF(const Chain& chain);


MLSimulation::MLSimulation(Data * pData, Model * pModel) :
	EpochSimulation(pData, pModel)
{
	this->lpChain = new Chain(pData);
	this->linsertDiagonalProbability = 0;
	this->lcancelDiagonalProbability = 0;
	this->lpermuteProbability = 0;
	this->linsertPermuteProbability = 0;
	this->ldeletePermuteProbability = 0;
	this->lrandomMissingProbability = 0;
	this->lmissingNetworkProbability = 0;
	this->lmissingBehaviorProbability = 0;
	this->laspect = NETWORK;
	for (int i = 0; i < 6; i++)
	{
		this->lrejections[i] = 0;
		this->lacceptances[i] = 0;
	}
}


MLSimulation::~MLSimulation()
{
	delete this->lpChain;
}


/**
 * Generates a random chain connecting the start and end observations of the
 * given data object for the given period. The chain is simple in the sense
 * that no two ministeps cancel each other out.
 */
void MLSimulation::connect(int period)
{
	this->lpChain->connect(period);
	this->initialize(period);
	this->updateProbabilities(this->lpChain,   this->lpChain->pFirst()->pNext(),
			this->lpChain->pLast()->pPrevious());
}


/**
 * Updates the probabilities for a range of ministeps of the given chain
 * (including the end-points of the range).
 */
void MLSimulation::updateProbabilities(const Chain * pChain,
	MiniStep * pFirstMiniStep,
	MiniStep * pLastMiniStep)
{
	// Initialize the variables as of the beginning of the period
    this->initialize(pChain->period());

    // Apply the ministeps before the first ministep of the required range
    // to derive the correct state before that ministep.

    this->executeMiniSteps(pChain->pFirst()->pNext(), pFirstMiniStep);

    bool done = false;
    MiniStep * pMiniStep = pFirstMiniStep;

	// set up array to store counts of structurally active ministeps by variable
	int *counts = new int[this->lvariables.size()];
	for (unsigned i = 0; i < this->lvariables.size(); i++)
	{
		counts[i] = 0;
	}
	while (!done)
    {
    	DependentVariable * pVariable =
    		this->lvariables[pMiniStep->variableId()];
		this->calculateRates();
		double rate = pVariable->rate(pMiniStep->ego());
    	double probability = pVariable->probability(pMiniStep);
    	double reciprocalTotalRate = 1 / this->totalRate();

		if (!pVariable->structural(pMiniStep))
		{
			counts[pMiniStep->variableId()] ++;
		}

    	pMiniStep->reciprocalRate(reciprocalTotalRate);
		pMiniStep->logOptionSetProbability(log(rate * reciprocalTotalRate));
		pMiniStep->logChoiceProbability(log(probability));
		pMiniStep->makeChange(pVariable);

    	if (pMiniStep == pLastMiniStep)
    	{
    		done = true;
    	}
    	else
    	{
    		pMiniStep = pMiniStep->pNext();
    	}
    }
	// calculate rate scores for the dependent variables
	for (unsigned i = 0; i < this->lvariables.size(); i++)
	{
		this->lvariables[i]->calculateMaximumLikelihoodRateScores(counts[i]);
	}

	// if necessary, calculate rate derivatives for the dependent variables
	if(this->pModel()->needDerivatives())
	{
		for (unsigned i = 0; i < this->lvariables.size(); i++)
		{
			this->lvariables[i]->
				calculateMaximumLikelihoodRateDerivatives(counts[i]);
		}

	}
	delete [] counts;
}

/**
 *  Does a few steps to increase the length of a newly constructed
 *  minimal chain.
 */

void MLSimulation::preburnin()
{
	int rejectCount = 0;
	bool accept;
	while (rejectCount < 5)
	{
		accept = this->insertDiagonalMiniStep();
		if (!accept)
		{
			rejectCount ++;
		}
	}
	rejectCount = 0;
	while (rejectCount < 5)
	{
		accept = this->insertPermute(1);
		if (!accept)
		{
			rejectCount ++;
		}
	}
}


/**
 *  Merges the probabilities, and does the required number of steps
 *  in the current chain
 *
 */

void MLSimulation::runEpoch(int period)
{
    // Initialize the rate functions of all variables and parameters for effects
	for (unsigned i = 0; i < this->lvariables.size(); i++)
	{
     	this->lvariables[i]->initializeRateFunction();
		this->lvariables[i]->updateEffectParameters();
	}
	this->setUpProbabilityArray();
	this->initialize(period);
	//	PrintValue(getMiniStepDF(*this->lpChain->pFirst()->pNext()));
	int numSteps = this->pModel()->numberMLSteps() ;

	for (int i = 0; i < numSteps; i++)
	{
		this->MLStep();
	}
}

const Chain * MLSimulation::pChain() const
{
	return this->lpChain;
}

void MLSimulation::pChain(Chain * pChain)
{
	this->lpChain = pChain;
}

/*
 * Set the chain to one uploaded from R and calculate the chain probabilities.
 *
 */

void MLSimulation::pChainProbabilities(Chain * pChain, int period)
{
	// pretty inefficient way to do this. probably need at least
	// option to store the change contributions on the ministep.
	delete this->lpChain;
	this->lpChain = pChain;
	this->runEpoch(period);

	this->updateProbabilities(this->lpChain, this->lpChain->pFirst()->pNext(),
			this->lpChain->pLast()->pPrevious());
}

/*
 * Clears the stores for MCMC.
 *
 */
void MLSimulation::initializeMCMCcycle()
{

	// clear storage for the sampled parameters.
	this->lBayesAcceptances.clear();
	this->lsampledBasicRates.clear();
	this->lcandidates.clear();

}

/*
 * Sample from the priors for the parameters.
 *
 */
void MLSimulation::MHPstep()
{
	double scaleFactor = this->pModel()->BayesianScaleFactor();
	double probabilityRatio;

	// first count how many steps for each variable and accumulate the log probs
	int * stepCount = new int[this->lvariables.size()];
	for (unsigned i = 0; i < this->lvariables.size(); i++)
	{
		stepCount[i] = 0;
	}
	double oldLogLikelihood = 0;
	double newLogLikelihood = 0;
	MiniStep * pMiniStep = this->pChain()->pFirst()->pNext();
	MiniStep * pLastMiniStep = this->pChain()->pLast();
	while (pMiniStep!= pLastMiniStep)
	{
		stepCount[pMiniStep->variableId()] ++;
		oldLogLikelihood += pMiniStep->logOptionSetProbability() +
			pMiniStep->logChoiceProbability();
		pMiniStep = pMiniStep->pNext();
	}
	double priorRatio = 0;

	for (unsigned i = 0; i < this->lvariables.size(); i++)
	{
		// generate gamma random variate for each basic rate
		this->lvariables[i]->sampleBasicRate(stepCount[i]);

		// generate normal random variate for each other effect
		priorRatio  += this->lvariables[i]->sampleParameters(scaleFactor);
	}

	//  calculate the acceptance probability for the non basic rate effects
	this->updateProbabilities(this->lpChain, this->lpChain->pFirst()->pNext(),
			this->lpChain->pLast()->pPrevious());
	pMiniStep = this->lpChain->pFirst()->pNext();
	while (pMiniStep!= pLastMiniStep)
	{
		newLogLikelihood += pMiniStep->logOptionSetProbability() +
			pMiniStep->logChoiceProbability();
		pMiniStep = pMiniStep->pNext();
	}

	probabilityRatio = exp(-oldLogLikelihood  + newLogLikelihood + priorRatio);

	if (nextDouble() < probabilityRatio)
	{
		this->lBayesAcceptances.push_back(true);
		// copy these parameters to the effectInfo's so they are the
		// default next time we reject
		for (unsigned i = 0; i < this->lvariables.size(); i++)
		{
			this->lvariables[i]->updateEffectInfoParameters();
		}
	}
	else
	{
		this->lBayesAcceptances.push_back(false);
		// re-initialize the rate functions of all variables and parameters
		// for effects from those in the effectinfo's
		for (unsigned i = 0; i < this->lvariables.size(); i++)
		{
			this->lvariables[i]->initializeRateFunction();
			this->lvariables[i]->updateEffectParameters();
		}
		// reset the chain probabilities
		this->updateProbabilities(this->lpChain,
			this->lpChain->pFirst()->pNext(),
			this->lpChain->pLast()->pPrevious());
	}
	delete[] stepCount;
}

/**
 *  sets up probability array
 *
 */

void MLSimulation::setUpProbabilityArray()
{
	this->lprobabilityArray[0] = this->linsertDiagonalProbability;
	this->lprobabilityArray[1] = this->lcancelDiagonalProbability;
	this->lprobabilityArray[2] = this->lpermuteProbability;
	this->lprobabilityArray[3] = this->linsertPermuteProbability;
	this->lprobabilityArray[4] = this->ldeletePermuteProbability;
	this->lprobabilityArray[5] = this->lrandomMissingProbability;
	for (int i = 0; i < 6; i++)
	{
		this->lrejections[i] = 0;
		this->lacceptances[i] = 0;
	}
}

/**
 *  Does one step in the current chain
 *
 */

void MLSimulation::MLStep()
{

	int stepType = nextIntWithProbabilities(6, this->lprobabilityArray);
	int c0 = 40;
	bool accept = false;
	switch (stepType)
	{
	case 0:
		accept = this->insertDiagonalMiniStep();
		break;
	case 1:
		accept = this->cancelDiagonalMiniStep();
		break;
	case 2:
		accept = this->permute(c0);
		break;
	case 3:
		accept = this->insertPermute(c0);
		break;
	case 4:
		accept = this->deletePermute(c0);
		break;
	case 5:
		//accept = this->randomMissing();
		break;
	}
	if (accept)
	{
		this->lacceptances[stepType]++;
	}
	else
	{
		this->lrejections[stepType]++;
	}
}
/**
 * Executes the given subsequence of ministeps excluding the last
 * ministep.
 */
void MLSimulation::executeMiniSteps(MiniStep * pFirstMiniStep,
	MiniStep * pLastMiniStep)
{
	MiniStep * pMiniStep = pFirstMiniStep;

	while (pMiniStep != pLastMiniStep)
	{
		DependentVariable * pVariable =
			this->lvariables[pMiniStep->variableId()];
		pMiniStep->makeChange(pVariable);

		pMiniStep = pMiniStep->pNext();
	}
}


void MLSimulation::setStateBefore(MiniStep * pMiniStep)
{
	this->resetVariables();
	this->executeMiniSteps(this->lpChain->pFirst()->pNext(), pMiniStep);
}


void MLSimulation::resetVariables()
{
	// Initialize each dependent variable

	for (unsigned i = 0; i < this->rVariables().size(); i++)
	{
		this->rVariables()[i]->initialize(this->period());
	}
}
/**
 * Returns the Bayes Acceptance for the given iteration.
 */

int MLSimulation::BayesAcceptances(unsigned iteration) const
{
	if (iteration < this->lBayesAcceptances.size())
	{
		return this->lBayesAcceptances[iteration];
	}
	else
	{
		throw std::out_of_range("The number" + toString(iteration) +
			" is not in the range [0," +
			toString(this->lBayesAcceptances.size()) + "].");
	}

}
/**
 * Returns the sampled basic rate parameter for the given iteration.
 */

double MLSimulation::sampledBasicRates(unsigned iteration) const
{
	if (iteration < this->lsampledBasicRates.size())
	{
		return this->lsampledBasicRates[iteration];
	}
	else
	{
		throw std::out_of_range("The number" + toString(iteration) +
			" is not in the range [0," +
			toString(this->lsampledBasicRates.size()) + "].");
	}

}
/**
 * Stores the sampled basic rate parameter for the next iteration.
 */

void MLSimulation::sampledBasicRates(double value)
{
	this->lsampledBasicRates.push_back(value);
}
/**
 * Returns the shape parameter used in the sampled basic rate parameter for
 * the given iteration.
 */

int MLSimulation::sampledBasicRatesDistributions(unsigned iteration) const
{
	if (iteration < this->lsampledBasicRatesDistributions.size())
	{
		return this->lsampledBasicRatesDistributions[iteration];
	}
	else
	{
		throw std::out_of_range("The number" + toString(iteration) +
			" is not in the range [0," +
			toString(this->lsampledBasicRatesDistributions.size()) + "].");
	}

}
/**
 * Stores the shape parameter used in the sampled basic rate parameter for the
 * next iteration.
 */

void MLSimulation::sampledBasicRatesDistributions(int value)
{
	this->lsampledBasicRatesDistributions.push_back(value);
}

/**
 * Returns the candidate value for the given iteration for the given effect.
 * The candidate values are updated in the MHPstep of a Bayesian simulation.
 */
double MLSimulation::candidates(const EffectInfo * pEffect,
	unsigned iteration) const
{
	map<const EffectInfo *, vector<double> >::const_iterator iter =
		this->lcandidates.find(pEffect);
	double candidate = 0;

	if (iter != this->lcandidates.end())
	{
		if (iteration < iter->second.size())
		{
			candidate = iter->second[iteration];
		}
		else
		{
			throw std::out_of_range("The number" + toString(iteration) +
				" is not in the range [0," +
				toString(iter->second.size()) + "].");
		}
	}

	return candidate;
}
/**
 * Stores the candidate value for the next iteration for the given effect.
 * The candidate values are updated in the MHPstep of a Bayesian simulation.
 */
void MLSimulation::candidates(const EffectInfo * pEffect, double value)
{
	this->lcandidates[pEffect].push_back(value);
}

// ----------------------------------------------------------------------------
// Section: Metropolis-Hastings steps
// ----------------------------------------------------------------------------

bool MLSimulation::insertDiagonalMiniStep()
{
	bool accept = false;
	MiniStep * pMiniStep = this->lpChain->randomMiniStep();
	this->setStateBefore(pMiniStep);
	this->calculateRates();
	DependentVariable * pVariable = this->chooseVariable();
	int i = this->chooseActor(pVariable);
	BehaviorVariable * pBehaviorVariable =
		dynamic_cast<BehaviorVariable *>(pVariable);
	if (!pVariable->pActorSet()->active(i) ||
		(pBehaviorVariable && pBehaviorVariable->structural(i)))
	{
		return false;
	}

	MiniStep * pNewMiniStep = 0;

	if (pBehaviorVariable)
	{
		pNewMiniStep =
			new BehaviorChange(
				dynamic_cast<BehaviorLongitudinalData *>(pVariable->pData()),
				i,
				0);
	}
	else
	{
		pNewMiniStep =
			new NetworkChange(
				dynamic_cast<NetworkLongitudinalData *>(pVariable->pData()),
				i,
				i);
	}

	double rr = 1 / this->totalRate();
	pNewMiniStep->reciprocalRate(rr);
	pNewMiniStep->logOptionSetProbability(log(pVariable->rate(i) * rr));

	double explcpr = pVariable->probability(pNewMiniStep);
	pNewMiniStep->logChoiceProbability(log(explcpr));

	double kappaFactor;

	if (this->lsimpleRates)
	{
		kappaFactor = 1 / (rr * (this->lpChain->ministepCount() + 1));
	}
	else
	{
		double sigma2 = this->lpChain->sigma2();
		double mu = this->lpChain->mu();

		kappaFactor = sqrt(sigma2 / (sigma2 + rr * rr)) *
			exp((1 - mu) * (1 - mu) / (2 * sigma2) -
				(1 - mu - rr) * (1 - mu - rr) / (2 * (sigma2 + rr * rr)));
	}

	this->lproposalProbability =
		kappaFactor *
		explcpr *
		this->lpChain->ministepCount() *
		this->lcancelDiagonalProbability /
		((this->lpChain->diagonalMinistepCount() + 1) *
			this->linsertDiagonalProbability);

	if (this->lproposalProbability > 1)
	{
		this->lproposalProbability = 1;
	}

	if (nextDouble() < this->lproposalProbability)
	{
		accept = true;
		this->lpChain->insertBefore(pNewMiniStep, pMiniStep);
	}

	return accept;
}


bool MLSimulation::cancelDiagonalMiniStep()
{
	bool accept = false;

	MiniStep * pMiniStep = this->lpChain->randomDiagonalMiniStep();
	double rr = pMiniStep->reciprocalRate();

	double kappaFactor;

	if (this->lsimpleRates)
	{
		kappaFactor = rr * this->lpChain->ministepCount();
	}
	else
	{
		double sigma2 = this->lpChain->sigma2();
		double mu = this->lpChain->mu();

		kappaFactor = sqrt(sigma2 / (sigma2 + rr * rr)) *
			exp((1 - mu) * (1 - mu) / (2 * sigma2) -
				(1 - mu + rr) * (1 - mu + rr) / (2 * (sigma2 - rr * rr)));
	}

	this->lproposalProbability =
		kappaFactor * exp(-pMiniStep->logChoiceProbability()) *
			this->lpChain->diagonalMinistepCount() *
			this->linsertDiagonalProbability /
		((this->lpChain->ministepCount() + 1) *
			this->lcancelDiagonalProbability);

	if (this->lproposalProbability > 1)
	{
		this->lproposalProbability = 1;
	}

	if (nextDouble() < this->lproposalProbability)
	{
		accept = true;
		this->lpChain->remove(pMiniStep);
	}

	return accept;
}


bool MLSimulation::permute(int c0)
{
	if (this->lpChain->ministepCount() <= 2)
	{
		return false;
	}

	bool accept = false;

	MiniStep * pMiniStepA = this->lpChain->randomMiniStep();

	while (pMiniStepA == this->lpChain->pLast())
	{
		pMiniStepA = this->lpChain->randomMiniStep();
	}

	vector<MiniStep *> interval;
	MiniStep * pMiniStep = pMiniStepA;

	while ((int) interval.size() < c0 && pMiniStep != this->lpChain->pLast())
	{
		interval.push_back(pMiniStep);
		pMiniStep = pMiniStep->pNext();
	}

	if (interval.size() <= 1)
	{
		return false;
	}

	MiniStep * pNextMiniStep = pMiniStep;

	permuteVector(interval);
	this->setStateBefore(pMiniStepA);
	bool valid = true;
	double sumlprob = 0;
	double sumlprob_new = 0;
	double mu_new = this->lpChain->mu();
	double sigma2_new = this->lpChain->sigma2();
	double * newReciprocalRate = new double[interval.size()];
	double * newOptionSetProbability = new double[interval.size()];
	double * newChoiceProbability = new double[interval.size()];

	for (unsigned i = 0; i < interval.size() && valid; i++)
	{
		pMiniStep = interval[i];
		DependentVariable * pVariable =
			this->lvariables[pMiniStep->variableId()];

		if (!pVariable->validMiniStep(pMiniStep))
		{
			valid = false;
		}
		else
		{
			sumlprob += pMiniStep->logChoiceProbability() +
				pMiniStep->logOptionSetProbability();
			double rrOld = pMiniStep->reciprocalRate();

			if (!this->simpleRates())
			{
				mu_new -= rrOld;
				sigma2_new -= rrOld * rrOld;
			}

			this->calculateRates();
			double rr = 1 / this->totalRate();
			double lospr =
				log(pVariable->rate(pMiniStep->ego()) * rr);
			double lcpr = log(pVariable->probability(pMiniStep));

			sumlprob_new += lospr + lcpr;

			if (!this->simpleRates())
			{
				mu_new += rr;
				sigma2_new += rr * rr;
			}

			pMiniStep->makeChange(pVariable);

			newReciprocalRate[i] = rr;
			newOptionSetProbability[i] = lospr;
			newChoiceProbability[i] = lcpr;
		}
	}

	if (valid)
	{
		double kappaFactor = 1;

		if (!this->simpleRates())
		{
			double sigma2 = this->lpChain->sigma2();
			double mu = this->lpChain->mu();

			kappaFactor = sqrt(sigma2 / sigma2_new) *
				exp((1 - mu) * (1 - mu) / (2 * sigma2) -
					(1 - mu_new) * (1 - mu_new) / (2 * sigma2_new));
		}

		this->lproposalProbability =
			kappaFactor * exp(sumlprob_new - sumlprob);

		if (this->lproposalProbability > 1)
		{
			this->lproposalProbability = 1;
		}

		if (nextDouble() < this->lproposalProbability)
		{
			accept = true;

			for (unsigned i = 0; i < interval.size(); i++)
			{
				pMiniStep = interval[i];

				this->lpChain->remove(pMiniStep);

				pMiniStep->reciprocalRate(newReciprocalRate[i]);
				pMiniStep->logOptionSetProbability(
					newOptionSetProbability[i]);
				pMiniStep->logChoiceProbability(newChoiceProbability[i]);
			}

			for (unsigned i = 0; i < interval.size(); i++)
			{
				this->lpChain->insertBefore(interval[i], pNextMiniStep);
			}
		}
	}

	delete[] newReciprocalRate;
	delete[] newOptionSetProbability;
	delete[] newChoiceProbability;

	return accept;
}


bool MLSimulation::insertPermute(int c0)
{
	if (this->lpChain->ministepCount() <= 1)
	{
		return false;
	}

	bool accept = false;

	MiniStep * pMiniStepA = this->lpChain->randomMiniStep();

	while (pMiniStepA == this->lpChain->pLast())
	{
		pMiniStepA = this->lpChain->randomMiniStep();
	}

	this->setStateBefore(pMiniStepA);
	this->calculateRates();
	DependentVariable * pVariable = this->chooseVariable();
	int i = this->chooseActor(pVariable);

	if (!pVariable->pActorSet()->active(i))
	{
		return false;
	}

	MiniStep * pLeftMiniStep = pVariable->randomMiniStep(i);

	if (pLeftMiniStep->diagonal() ||
		!pVariable->validMiniStep(pLeftMiniStep))
	{
		delete pLeftMiniStep;
		return false;
	}

	double rr0 = 1 / this->totalRate();
	double lospr0 = log(pVariable->rate(i) * rr0);
	double lcpr0 = pLeftMiniStep->logChoiceProbability();

	pLeftMiniStep->reciprocalRate(rr0);
	pLeftMiniStep->logOptionSetProbability(lospr0);

	bool misdat = pLeftMiniStep->missing(this->lpChain->period());
	this->lmissingData = misdat;
	MiniStep * pMiniStepB = this->lpChain->pLast();
	double choiceLength = 1;

	if (!misdat)
	{
		MiniStep * pMiniStepD =
			this->lpChain->nextMiniStepForOption(*pLeftMiniStep->pOption(),
				pMiniStepA);

		if (!pMiniStepD)
		{
			pMiniStepD = this->lpChain->pLast();
		}

		if (pMiniStepD == pMiniStepA)
		{
			delete pLeftMiniStep;
			return false;
		}

		choiceLength =
			this->lpChain->intervalLength(pMiniStepA, pMiniStepD) - 1;
		pMiniStepB =
			this->lpChain->randomMiniStep(pMiniStepA->pNext(), pMiniStepD);
	}

	vector<MiniStep *> interval;
	MiniStep * pMiniStep = pMiniStepA;

	while ((int) interval.size() < c0 && pMiniStep != pMiniStepB)
	{
		interval.push_back(pMiniStep);
		pMiniStep = pMiniStep->pNext();
	}

	permuteVector(interval);

	// Add the ministeps up to the pMiniStepB to the interval,
	// as the calculations with these ministeps are the same as with
	// the permuted ones.

	while (pMiniStep != pMiniStepB)
	{
		interval.push_back(pMiniStep);
		pMiniStep = pMiniStep->pNext();
	}

	bool valid = true;
	double sumlprob = 0;
	double sumlprob_new = 0;
	double mu_new = this->lpChain->mu();
	double sigma2_new = this->lpChain->sigma2();
	double * newReciprocalRate = new double[interval.size()];
	double * newOptionSetProbability = new double[interval.size()];
	double * newChoiceProbability = new double[interval.size()];

	// We execute the new chain, do the necessary calculations, and
	// simultaneously check if the new chain is valid.

	// Execute the new ministep that is inserted before pMiniStepA

	sumlprob_new += lcpr0 + lospr0;

	if (!this->simpleRates())
	{
		mu_new += rr0;
		sigma2_new += rr0 * rr0;
	}

	pLeftMiniStep->makeChange(pVariable);

	// Execute the existing ministeps (with some of them in a permuted order)

	for (unsigned i = 0; i < interval.size() && valid; i++)
	{
		pMiniStep = interval[i];
		DependentVariable * pVariable =
			this->lvariables[pMiniStep->variableId()];

		if (!pVariable->validMiniStep(pMiniStep))
		{
			valid = false;
		}
		else
		{
			sumlprob += pMiniStep->logChoiceProbability() +
				pMiniStep->logOptionSetProbability();
			double rrOld = pMiniStep->reciprocalRate();

			if (!this->simpleRates())
			{
				mu_new -= rrOld;
				sigma2_new -= rrOld * rrOld;
			}

			this->calculateRates();
			double rr = 1 / this->totalRate();
			double lospr =
				log(pVariable->rate(pMiniStep->ego()) * rr);
			double lcpr = log(pVariable->probability(pMiniStep));

			sumlprob_new += lospr + lcpr;

			if (!this->simpleRates())
			{
				mu_new += rr;
				sigma2_new += rr * rr;
			}

			pMiniStep->makeChange(pVariable);

			newReciprocalRate[i] = rr;
			newOptionSetProbability[i] = lospr;
			newChoiceProbability[i] = lcpr;
		}
	}

	// Execute the canceling ministep of the generated CCP.

	MiniStep * pRightMiniStep = 0;

	if (valid && !misdat)
	{
		pRightMiniStep = pLeftMiniStep->createReverseMiniStep();

		this->calculateRates();
		double rr = 1 / this->totalRate();
		double lospr =
			log(pVariable->rate(pRightMiniStep->ego()) * rr);
		double lcpr = log(pVariable->probability(pRightMiniStep));

		sumlprob_new += lospr + lcpr;

		if (!this->simpleRates())
		{
			mu_new += rr;
			sigma2_new += rr * rr;
		}

		pRightMiniStep->makeChange(pVariable);
		pRightMiniStep->reciprocalRate(rr);
		pRightMiniStep->logOptionSetProbability(lospr);
		pRightMiniStep->logChoiceProbability(lcpr);
	}

	if (valid)
	{
		double pr1 = 0;

		if (!misdat)
		{
			pr1 =
				(1 -
					this->lmissingNetworkProbability -
					this->lmissingBehaviorProbability) /
				(this->lpChain->consecutiveCancelingPairCount() + 1);
		}
		else
		{
			if (pVariable->networkVariable())
			{
				pr1 = this->lmissingNetworkProbability /
					(this->lpChain->missingNetworkMiniStepCount() + 1);
			}
			else
			{
				// pVariable is a behavior variable

				pr1 = this->lmissingBehaviorProbability /
					(this->lpChain->missingBehaviorMiniStepCount() + 1);
			}
		}

		double kappaFactor = 0;

		if (this->simpleRates())
		{
			if (!misdat)
			{
				kappaFactor = 1 /
					(rr0 * rr0 *
						(this->lpChain->ministepCount() + 1) *
						(this->lpChain->ministepCount() + 2));
			}
			else
			{
				kappaFactor = 1 / (rr0 * (this->lpChain->ministepCount() + 1));
			}
		}
		else
		{
			double sigma2 = this->lpChain->sigma2();
			double mu = this->lpChain->mu();

			kappaFactor = sqrt(sigma2 / sigma2_new) *
				exp((1 - mu) * (1 - mu) / (2 * sigma2) -
					(1 - mu_new) * (1 - mu_new) / (2 * sigma2_new));
		}

		this->lproposalProbability =
			kappaFactor *
				exp(sumlprob_new - sumlprob) *
				this->ldeletePermuteProbability *
				pr1 *
				(this->lpChain->ministepCount() + 1) *
				choiceLength /
			(this->linsertPermuteProbability * exp(lospr0 + lcpr0));

		if (this->lproposalProbability > 1)
		{
			this->lproposalProbability = 1;
		}

		if (nextDouble() < this->lproposalProbability)
		{
			// Change the chain permanently

			accept = true;

			for (unsigned i = 0; i < interval.size(); i++)
			{
				pMiniStep = interval[i];

				this->lpChain->remove(pMiniStep);

				pMiniStep->reciprocalRate(newReciprocalRate[i]);
				pMiniStep->logOptionSetProbability(
					newOptionSetProbability[i]);
				pMiniStep->logChoiceProbability(newChoiceProbability[i]);
			}

			this->lpChain->insertBefore(pLeftMiniStep, pMiniStepB);

			for (unsigned i = 0; i < interval.size(); i++)
			{
				this->lpChain->insertBefore(interval[i], pMiniStepB);
			}

			if (!misdat)
			{
				this->lpChain->insertBefore(pRightMiniStep, pMiniStepB);
			}
		}
		else
		{
			// Operation not accepted, so we don't need the generated
			// pair of ministeps.

			delete pLeftMiniStep;
			delete pRightMiniStep;
		}
	}

	delete[] newReciprocalRate;
	delete[] newOptionSetProbability;
	delete[] newChoiceProbability;

	return accept;
}


bool MLSimulation::deletePermute(int c0)
{
	this->lmissingData =
		nextDouble() <
			this->lmissingNetworkProbability +
				this->lmissingBehaviorProbability;
	MiniStep * pMiniStepA = 0;
	MiniStep * pMiniStepB = 0;

	if (!this->lmissingData)
	{
		if (this->lpChain->consecutiveCancelingPairCount() == 0)
		{
			return false;
		}

		pMiniStepA = this->lpChain->randomConsecutiveCancelingPair();
		pMiniStepB = pMiniStepA->pNextWithSameOption();
	}
	else
	{
		if (nextDouble() <
			this->lmissingNetworkProbability /
				(this->lmissingNetworkProbability +
					this->lmissingBehaviorProbability))
		{
			if (this->lpChain->missingNetworkMiniStepCount() == 0)
			{
				return false;
			}

			pMiniStepA = this->lpChain->randomMissingNetworkMiniStep();
		}
		else
		{
			if (this->lpChain->missingBehaviorMiniStepCount() == 0)
			{
				return false;
			}

			pMiniStepA = this->lpChain->randomMissingBehaviorMiniStep();
		}

		pMiniStepB = this->lpChain->pLast();
	}

	if (pMiniStepA->networkMiniStep())
	{
		this->laspect = NETWORK;
	}
	else
	{
		this->laspect = BEHAVIOR;
	}

	// Create the permuted interval of ministeps

	vector<MiniStep *> interval;
	MiniStep * pMiniStep = pMiniStepA->pNext();

	while ((int) interval.size() < c0 && pMiniStep != pMiniStepB)
	{
		interval.push_back(pMiniStep);
		pMiniStep = pMiniStep->pNext();
	}

	permuteVector(interval);

	// Add the ministeps up to pMiniStepB to the interval,
	// as the calculations with these ministeps are the same as with
	// the permuted ones.

	while (pMiniStep != pMiniStepB)
	{
		interval.push_back(pMiniStep);
		pMiniStep = pMiniStep->pNext();
	}

	this->setStateBefore(pMiniStepA);

	bool valid = true;

	double sumlprob = pMiniStepA->logChoiceProbability() +
		pMiniStepA->logOptionSetProbability();
	double sumlprob_new = 0;
	double rr = pMiniStepA->reciprocalRate();
	double mu_new = this->lpChain->mu() - rr;
	double sigma2_new = this->lpChain->sigma2() - rr * rr;
	double * newReciprocalRate = new double[interval.size()];
	double * newOptionSetProbability = new double[interval.size()];
	double * newChoiceProbability = new double[interval.size()];

	if (!this->lmissingData)
	{
		sumlprob += pMiniStepB->logChoiceProbability() +
			pMiniStepB->logOptionSetProbability();

		if (!this->simpleRates())
		{
			rr = pMiniStepB->reciprocalRate();
			mu_new -= rr;
			sigma2_new -= rr * rr;
		}
	}

	// We execute the new chain, do the necessary calculations, and
	// simultaneously check if the new chain is valid.

	for (unsigned i = 0; i < interval.size() && valid; i++)
	{
		pMiniStep = interval[i];
		DependentVariable * pVariable =
			this->lvariables[pMiniStep->variableId()];

		if (!pVariable->validMiniStep(pMiniStep))
		{
			valid = false;
		}
		else
		{
			sumlprob += pMiniStep->logChoiceProbability() +
				pMiniStep->logOptionSetProbability();
			double rrOld = pMiniStep->reciprocalRate();

			if (!this->simpleRates())
			{
				mu_new -= rrOld;
				sigma2_new -= rrOld * rrOld;
			}

			this->calculateRates();
			double rr = 1 / this->totalRate();
			double lospr =
				log(pVariable->rate(pMiniStep->ego()) * rr);
			double lcpr = log(pVariable->probability(pMiniStep));

			sumlprob_new += lospr + lcpr;

			if (!this->simpleRates())
			{
				mu_new += rr;
				sigma2_new += rr * rr;
			}

			pMiniStep->makeChange(pVariable);

			newReciprocalRate[i] = rr;
			newOptionSetProbability[i] = lospr;
			newChoiceProbability[i] = lcpr;
		}
	}

	bool accept = false;

	if (valid)
	{
		double pr1 = 0;

		if (!this->lmissingData)
		{
			pr1 =
				(1 -
					this->lmissingNetworkProbability -
					this->lmissingBehaviorProbability) /
				this->lpChain->consecutiveCancelingPairCount();
		}
		else
		{
			if (pMiniStepA->networkMiniStep())
			{
				pr1 = this->lmissingNetworkProbability /
					this->lpChain->missingNetworkMiniStepCount();
			}
			else
			{
				// pVariable is a behavior variable

				pr1 = this->lmissingBehaviorProbability /
					this->lpChain->missingBehaviorMiniStepCount();
			}
		}

		// Calculate choiceLength

		double choiceLength = 1;

		if (!this->lmissingData)
		{
			MiniStep * pMiniStepD = pMiniStepB->pNextWithSameOption();

			if (!pMiniStepD)
			{
				pMiniStepD = this->lpChain->pLast();
			}

			choiceLength =
				this->lpChain->intervalLength(pMiniStepA, pMiniStepD) - 3;
		}

		double lpr0 = pMiniStepA->logChoiceProbability() +
			pMiniStepA->logOptionSetProbability();

		double kappaFactor = 0;

		if (this->simpleRates())
		{
			if (!this->lmissingData)
			{
				kappaFactor = rr * rr *
					this->lpChain->ministepCount() *
					(this->lpChain->ministepCount() - 1);
			}
			else
			{
				kappaFactor = rr * this->lpChain->ministepCount();
			}
		}
		else
		{
			double sigma2 = this->lpChain->sigma2();
			double mu = this->lpChain->mu();

			kappaFactor = sqrt(sigma2 / sigma2_new) *
				exp((1 - mu) * (1 - mu) / (2 * sigma2) -
					(1 - mu_new) * (1 - mu_new) / (2 * sigma2_new));
		}

		this->lproposalProbability =
			kappaFactor *
				exp(sumlprob_new - sumlprob) *
				this->linsertPermuteProbability *
				exp(lpr0) /
			(this->ldeletePermuteProbability *
				pr1 *
				(this->lpChain->ministepCount() + 1) *
				choiceLength);

		if (this->lproposalProbability > 1)
		{
			this->lproposalProbability = 1;
		}

		if (nextDouble() < this->lproposalProbability)
		{
			// Change the chain permanently

			accept = true;
			this->lpChain->remove(pMiniStepA);

			for (unsigned i = 0; i < interval.size(); i++)
			{
				pMiniStep = interval[i];

				this->lpChain->remove(pMiniStep);

				pMiniStep->reciprocalRate(newReciprocalRate[i]);
				pMiniStep->logOptionSetProbability(
					newOptionSetProbability[i]);
				pMiniStep->logChoiceProbability(newChoiceProbability[i]);
			}

			for (unsigned i = 0; i < interval.size(); i++)
			{
				this->lpChain->insertBefore(interval[i], pMiniStepB);
			}

			if (!this->lmissingData)
			{
				this->lpChain->remove(pMiniStepB);
			}
		}
	}

	delete[] newReciprocalRate;
	delete[] newOptionSetProbability;
	delete[] newChoiceProbability;

	return accept;
}


// ----------------------------------------------------------------------------
// Section: Accessors
// ----------------------------------------------------------------------------

/**
 * Returns the proposal probability of the last Metropolis-Hastings step.
 */
double MLSimulation::proposalProbability() const
{
	return this->lproposalProbability;
}


/**
 * Returns if the last operation has worked with ministeps having missing
 * data at either end of the period.
 */
bool MLSimulation::missingData() const
{
	return this->lmissingData;
}


/**
 * Returns the aspect of the ministeps involved in the last
 * Metropolis-Hastings step.
 */
Aspect MLSimulation::aspect() const
{
	return this->laspect;
}


/**
 * Stores if simple rates should be used in simulations.
 */
void MLSimulation::simpleRates(bool flag)
{
	this->lsimpleRates = flag;
}


/**
 * Returns if simple rates should be used in simulations.
 */
bool MLSimulation::simpleRates() const
{
	return this->lsimpleRates;
}


/**
 * Stores the probability associated with the insertDiagonalMiniStep
 * operation.
 */
void MLSimulation::insertDiagonalProbability(double probability)
{
	this->linsertDiagonalProbability = probability;
}


/**
 * Returns the probability associated with the insertDiagonalMiniStep
 * operation.
 */
double MLSimulation::insertDiagonalProbability() const
{
	return this->linsertDiagonalProbability;
}


/**
 * Stores the probability associated with the cancelDiagonalMiniStep
 * operation.
 */
void MLSimulation::cancelDiagonalProbability(double probability)
{
	this->lcancelDiagonalProbability = probability;
}


/**
 * Returns the probability associated with the cancelDiagonalMiniStep
 * operation.
 */
double MLSimulation::cancelDiagonalProbability() const
{
	return this->lcancelDiagonalProbability;
}


/**
 * Stores the probability associated with the permute
 * operation.
 */
void MLSimulation::permuteProbability(double probability)
{
	this->lpermuteProbability = probability;
}


/**
 * Returns the probability associated with the permute
 * operation.
 */
double MLSimulation::permuteProbability() const
{
	return this->lpermuteProbability;
}

/**
 * Stores the probability associated with the insertPermute
 * operation.
 */
void MLSimulation::insertPermuteProbability(double probability)
{
	this->linsertPermuteProbability = probability;
}


/**
 * Returns the probability associated with the insertPermute
 * operation.
 */
double MLSimulation::insertPermuteProbability() const
{
	return this->linsertPermuteProbability;
}


/**
 * Stores the probability associated with the deletePermute
 * operation.
 */
void MLSimulation::deletePermuteProbability(double probability)
{
	this->ldeletePermuteProbability = probability;
}


/**
 * Returns the probability associated with the deletePermute
 * operation.
 */
double MLSimulation::deletePermuteProbability() const
{
	return this->ldeletePermuteProbability;
}

/**
 * Stores the probability associated with the randomMissing
 * operation.
 */
void MLSimulation::randomMissingProbability(double probability)
{
	this->lrandomMissingProbability = probability;
}


/**
 * Returns the probability associated with the randomMissing
 * operation.
 */
double MLSimulation::randomMissingProbability() const
{
	return this->lrandomMissingProbability;
}

/**
 * Stores the missing network probability (prmin in Tom's specification).
 */
void MLSimulation::missingNetworkProbability(double probability)
{
	this->lmissingNetworkProbability = probability;
}


/**
 * Returns the missing network probability (prmin in Tom's specification).
 */
double MLSimulation::missingNetworkProbability() const
{
	return this->lmissingNetworkProbability;
}


/**
 * Stores the missing behavior probability (prmib in Tom's specification).
 */
void MLSimulation::missingBehaviorProbability(double probability)
{
	this->lmissingBehaviorProbability = probability;
}


/**
 * Returns the missing behavior probability (prmib in Tom's specification).
 */
double MLSimulation::missingBehaviorProbability() const
{
	return this->lmissingBehaviorProbability;
}


/**
 * Returns the number of acceptances for this steptype.
 */
int MLSimulation::acceptances(int stepType) const
{
	return this->lacceptances[stepType];
}

/**
 * Returns the number of rejections for this steptype.
 */
int MLSimulation::rejections(int stepType) const
{
	return this->lrejections[stepType];
}
}
