/******************************************************************************
 * SIENA: Simulation Investigation for Empirical Network Analysis
 *
 * Web: http://www.stats.ox.ac.uk/~snijders/siena/
 *
 * File: ConstantFunction.cpp
 *
 * Description: This file contains the implementation of the class
 * ConstantFunction.
 *****************************************************************************/

#include "ConstantFunction.h"

namespace siena
{

ConstantFunction::ConstantFunction(double constant)
{
	this->lconstant = constant;
}


double ConstantFunction::value(int alter)
{
	return this->lconstant;
}

}
