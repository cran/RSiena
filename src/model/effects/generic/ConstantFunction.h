/******************************************************************************
 * SIENA: Simulation Investigation for Empirical Network Analysis
 *
 * Web: http://www.stats.ox.ac.uk/~snijders/siena/
 *
 * File: ConstantFunction.h
 *
 * Description: This file contains the definition of the
 * ConstantFunction class.
 *****************************************************************************/

#ifndef CONSTANTFUNCTION_H_
#define CONSTANTFUNCTION_H_

#include "AlterFunction.h"

namespace siena
{

class ConstantFunction: public AlterFunction
{
public:
	ConstantFunction(double constant);

	virtual double value(int alter);

private:
	double lconstant;
};

}

#endif /* CONSTANTFUNCTION_H_ */
