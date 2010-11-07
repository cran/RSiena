/******************************************************************************
 * SIENA: Simulation Investigation for Empirical Network Analysis
 *
 * Web: http://www.stats.ox.ac.uk/~snijders/siena/
 *
 * File: CovariateDistance2AlterEffect.cpp
 *
 * Description: This file contains the implementation of the
 * CovariateDistance2AlterEffect class.
 *****************************************************************************/

#include "CovariateDistance2AlterEffect.h"
#include "network/Network.h"
#include "model/variables/NetworkVariable.h"
#include "model/EffectInfo.h"

namespace siena
{

/**
 * Constructor.
 * @param[in] pEffectInfo the effect descriptor
 */
CovariateDistance2AlterEffect::CovariateDistance2AlterEffect(const EffectInfo * pEffectInfo) :
		CovariateDistance2NetworkEffect(pEffectInfo)
{
}


/**
 * Calculates the contribution of a tie flip to the given actor.
 */
double CovariateDistance2AlterEffect::calculateContribution(int alter) const
{
	double parameter = this->pEffectInfo()->internalEffectParameter();
	double change = this->averageAlterValue(alter);
	if (parameter == 2)
	{
		int tieValue =  this->pNetwork()->tieValue(alter, this->ego());
		if (tieValue == 1)
		{
			int degree = this->pNetwork()->outDegree(alter);
			//Rprintf("%d %f %d %f\n", degree, change, this->ego(),
			//	this->value(this->ego()) );
			if (degree > 1)
			{
				change -= (degree * change - this->value(this->ego()))/
					(degree - 1);
			}
			else
			{
				change = 0;
			}
			//	Rprintf("aft er %d %f %d %f\n", degree, change, this->ego(),
			//	this->value(this->ego()) );
		}
	}
	return change;
}


/**
 * The contribution of the tie from the implicit ego to the given alter
 * to the statistic. It is assumed that preprocessEgo(ego) has been
 * called before.
 */
double CovariateDistance2AlterEffect::tieStatistic(int alter)
{
	double statistic = 0;
	double parameter = this->pEffectInfo()->internalEffectParameter();

	if (!this->missingDummy(alter))
	{
		statistic = this->averageAlterValue(alter);
		if (parameter == 2)
		{
			int tieValue =  this->pNetwork()->tieValue(alter, this->ego());
			if (tieValue == 1)
			{
				int degree = this->pNetwork()->outDegree(alter);
//				Rprintf("before %d %f %d %f\n", degree, statistic,
//							this->ego(), this->value(this->ego()) );
				if (degree > 1)
				{
					statistic -=
						(degree * statistic - this->value(this->ego()))/
						(degree - 1);
				}
				else
				{
					statistic = 0;
				}
//				Rprintf("stat after %d %f %d %f\n", degree, statistic,
//					this->ego(), this->value(this->ego()) );
			}
		}

	}
//			Rprintf("%f\n", statistic );
	return statistic;
}

}
