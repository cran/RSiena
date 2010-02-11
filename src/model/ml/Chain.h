/******************************************************************************
 * SIENA: Simulation Investigation for Empirical Network Analysis
 *
 * Web: http://www.stats.ox.ac.uk/~snijders/siena/
 *
 * File: Chain.h
 *
 * Description: This file contains the definition of the Chain class.
 *****************************************************************************/


#ifndef CHAIN_H_
#define CHAIN_H_

namespace siena
{

// ----------------------------------------------------------------------------
// Section: Forward declarations
// ----------------------------------------------------------------------------

class MiniStep;
class Data;


// ----------------------------------------------------------------------------
// Section: Class definition
// ----------------------------------------------------------------------------


/**
 * Defines a sequence of ministeps, which is the basic structure for
 * Maximum-Likelihood calculations.
 */
class Chain
{
public:
	Chain();
	virtual ~Chain();

	void emptyChain();
	void insertBefore(MiniStep * pNewMiniStep, MiniStep * pExistingMiniStep);
	void remove(MiniStep * pMiniStep);
	void connect(Data * pData, int period);

	int period() const;
	MiniStep * pFirst() const;
	MiniStep * pLast() const;

private:
	// A dummy first ministep in the chain
	MiniStep * lpFirst;

	// A dummy last ministep in the chain
	MiniStep * lpLast;

	// The period of changes represented by this chain
	int lperiod;
};

}

#endif /* CHAIN_H_ */
