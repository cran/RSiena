/******************************************************************************
 * SIENA: Simulation Investigation for Empirical Network Analysis
 *
 * Web: http://www.stats.ox.ac.uk/~snijders/siena/
 *
 * File: CovariateDistance2NetworkFunction.cpp
 *
 * Description: This file contains the implementation of the class
 * CovariateDistance2NetworkFunction.
 *****************************************************************************/

#include "CovariateDistance2NetworkFunction.h"
#include "network/Network.h"
#include "network/IncidentTieIterator.h"
#include "data/ConstantCovariate.h"
#include "data/ChangingCovariate.h"
#include "data/BehaviorLongitudinalData.h"

namespace siena
{

/**
 * Constructor.
 * @param[in] networkName the name of the network variable this function is
 * associated with
 * @param[in] covariateName the name of the covariate this function is
 * associated with
 */
CovariateDistance2NetworkFunction::CovariateDistance2NetworkFunction(
	string networkName, string covariateName) :
	CovariateNetworkAlterFunction(networkName, covariateName)
{
	this->laverageAlterValues = 0;
	this->ltotalAlterValues = 0;
	this->laverageAlterMissing = 0;
	this->laverageInAlterValues = 0;
	this->ltotalInAlterValues = 0;
	this->laverageInAlterMissing = 0;
}

CovariateDistance2NetworkFunction::~CovariateDistance2NetworkFunction()
{
	delete [] this->laverageAlterValues;
	delete [] this->ltotalAlterValues;
	delete [] this->laverageAlterMissing;
	delete [] this->laverageInAlterValues;
	delete [] this->ltotalInAlterValues;
	delete [] this->laverageInAlterMissing;
}

void CovariateDistance2NetworkFunction::initialize(const Data * pData,
	State * pState,
	int period,
	Cache * pCache)
{
	CovariateNetworkAlterFunction::initialize(pData, pState, period, pCache);

	if (this->laverageAlterValues)
	{
		delete [] this->laverageAlterValues;
	}
	if (this->ltotalAlterValues)
	{
		delete [] this->ltotalAlterValues;
	}
	if (this->laverageAlterMissing)
	{
		delete[] this->laverageAlterMissing;
	}
	if (this->laverageInAlterValues)
	{
		delete [] this->laverageInAlterValues;
	}
	if (this->ltotalInAlterValues)
	{
		delete [] this->ltotalInAlterValues;
	}
	if (this->laverageInAlterMissing)
	{
		delete[] this->laverageInAlterMissing;
	}
	this->laverageAlterValues = new double[this->pNetwork()->n()];
	this->ltotalAlterValues = new double[this->pNetwork()->n()];
	this->laverageAlterMissing = new bool[this->pNetwork()->n()];
	this->laverageInAlterValues = new double[this->pNetwork()->m()];
	this->ltotalInAlterValues = new double[this->pNetwork()->m()];
	this->laverageInAlterMissing = new bool[this->pNetwork()->m()];
}

/**
 * Returns the value of this function for the given alter. It is assumed
 * that the function has been initialized before and pre-processed with
 * respect to a certain ego.
 */
double CovariateDistance2NetworkFunction::averageAlterValue(int alter) const
{
	return this->laverageAlterValues[alter];
}

/**
 * Returns the value of this function for the given alter. It is assumed
 * that the function has been initialized before and pre-processed with
 * respect to a certain ego.
 */
double CovariateDistance2NetworkFunction::totalAlterValue(int alter) const
{
	return this->ltotalAlterValues[alter];
}

/**
 * Returns if the dummy covariate value for the given actor is based on
 * all missing values.
 */
bool CovariateDistance2NetworkFunction::missingDummy(int alter) const
{
	return this->laverageAlterMissing[alter];
}
/**
 * Returns the value of this function for the given alter. It is assumed
 * that the function has been initialized before and pre-processed with
 * respect to a certain ego.
 */

double CovariateDistance2NetworkFunction::averageInAlterValue(int alter) const
{
	return this->laverageInAlterValues[alter];
}

/**
 * Returns the value of this function for the given alter. It is assumed
 * that the function has been initialized before and pre-processed with
 * respect to a certain ego.
 */
double CovariateDistance2NetworkFunction::totalInAlterValue(int alter) const
{
	return this->ltotalInAlterValues[alter];
}

/**
 * Returns if the dummy covariate value for the given actor is based on
 * all missing values.
 */
bool CovariateDistance2NetworkFunction::missingInDummy(int alter) const
{
	return this->laverageInAlterMissing[alter];
}

/**
 * Does the necessary preprocessing work for calculating the
 * predicate for a specific ego. This method must be invoked before
 * calling AlterPredicate::value(...).
 */
void CovariateDistance2NetworkFunction::preprocessEgo(int ego)
{
	AlterFunction::preprocessEgo(ego);
	// set up the covariate based on current values of the network
	const Network * pNetwork = this->pNetwork();

	for (int i = 0; i < pNetwork->n(); i++)
	{
		int numberNonMissing = 0;
		this->laverageAlterMissing[i] = true;
		this->ltotalAlterValues[i] = 0;
		if (pNetwork->outDegree(i) > 0)
		{
			for (IncidentTieIterator iter = pNetwork->outTies(i);
				 iter.valid();
				 iter.next())
			{
				int j = iter.actor();
				this->ltotalAlterValues[i] += this->value(j);
				if (!this->missing(j))
				{
					numberNonMissing++;
				}
			}
			this->laverageAlterValues[i] =
					(this->ltotalAlterValues[i] / pNetwork->outDegree(i));
			if (numberNonMissing > 0)
			{
				this->laverageAlterValues[i] =
					(this->ltotalAlterValues[i] / numberNonMissing);
				this->laverageAlterMissing[i] = false;
			}
			else
			{
				this->laverageAlterValues[i] = covmean();
			}
		}
		else
		{
			this->laverageAlterValues[i] = covmean();
			this->ltotalAlterValues[i] = 0;
			this->laverageAlterMissing[i] = false;
		}
	}

	for (int i = 0; i < pNetwork->m(); i++)
	{
		int numberNonMissing = 0;
		this->laverageInAlterMissing[i] = true;
		this->ltotalInAlterValues[i] = 0;
		if (pNetwork->inDegree(i) > 0)
		{
			for (IncidentTieIterator iter = pNetwork->inTies(i);
				 iter.valid();
				 iter.next())
			{
				int j = iter.actor();
				this->ltotalInAlterValues[i] += this->value(j);
				if (!this->missing(j))
				{
					numberNonMissing++;
				}
			}
			if (numberNonMissing > 0)
			{
				this->laverageInAlterMissing[i] = false;
				this->laverageInAlterValues[i] =
					(this->ltotalInAlterValues[i] / numberNonMissing);
			}
			else
			{
				this->laverageInAlterValues[i] = covmean();
			}
		}
		else
		{
			this->laverageInAlterValues[i] = covmean();
			this->ltotalInAlterValues[i] = 0;
			this->laverageInAlterMissing[i] = false;
		}
	}
}
/**
 * Returns the centered similarity of the average alter values
   of the given actors wrt to the network with
 * which this function is associated.
 */
double CovariateDistance2NetworkFunction::similarityAvAlt(int i, int j) const
{
	double similarity = 0;

	if (this->pConstantCovariate())
	{
		similarity =
			this->pConstantCovariate()->similarity(
				this->averageAlterValue(i),
				this->averageAlterValue(j));
	}
	else if (this->pChangingCovariate())
	{
		similarity =
			this->pChangingCovariate()->similarity(
				this->averageAlterValue(i),
				this->averageAlterValue(j));
	}
	else
	{
		similarity =
			this->pBehaviorData()->similarity(
				this->averageAlterValue(i),
				this->averageAlterValue(j));
	}
	return similarity;
}
/**
 * Returns the centered similarity of the own value and the average alter value
   of the given actors wrt to the network with which this function is associated.
 */
double CovariateDistance2NetworkFunction::varOutAvSimilarity(int i, int j) const
{
	double similarity = 0;
	double outAlter = this->totalAlterValue(j);
	int tieValue = this->pNetwork()->tieValue(j, i);
	int degree = this->pNetwork()->outDegree(j);
	if (tieValue >= 1)
	{
		outAlter = outAlter - this->value(i);
		degree--;
	}

	if (degree >= 1)
	{
		outAlter /= degree;
	}
	else
	{
		outAlter = this->covmean();
	}

	if (this->pConstantCovariate())
	{
		similarity =
			this->pConstantCovariate()->similarity(this->value(i), outAlter);
	}
	else if (this->pChangingCovariate())
	{
		similarity =
			this->pChangingCovariate()->similarity(this->value(i), outAlter);
	}
	else
	{
		similarity =
			this->pBehaviorData()->similarity(this->value(i), outAlter);
	}

	return similarity;
}
/**
 * Returns the centered similarity of the own value and the average in-alter value
   of the given actors wrt to the network with
 * which this function is associated.
 */
double CovariateDistance2NetworkFunction::varInAvSimilarity(int i, int j) const
{
	double similarity = 0;
	double inAlter = this->totalInAlterValue(j);
	int tieValue = this->pNetwork()->tieValue(i, j);
	int degree = this->pNetwork()->inDegree(j);
	if (tieValue >= 1)
	{
		inAlter = inAlter - this->value(i);
		degree--;
	}

	if (degree >= 1)
	{
		inAlter /= degree;
	}
	else
	{
		inAlter = this->covmean();
	}

	if (this->pConstantCovariate())
	{
		similarity =
			this->pConstantCovariate()->similarity(this->value(i), inAlter);
	}
	else if (this->pChangingCovariate())
	{
		similarity =
			this->pChangingCovariate()->similarity(this->value(i), inAlter);
	}
	else
	{
		similarity =
			this->pBehaviorData()->similarity(this->value(i), inAlter);
	}

	return similarity;
}
}
