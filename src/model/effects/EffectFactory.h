/******************************************************************************
 * SIENA: Simulation Investigation for Empirical Network Analysis
 * 
 * Web: http://www.stats.ox.ac.uk/~snijders/siena/
 * 
 * File: EffectFactory.h
 * 
 * Description: This file contains the definition of the
 * EffectFactory class.
 *****************************************************************************/

#ifndef EFFECTFACTORY_H_
#define EFFECTFACTORY_H_

namespace siena
{

// ----------------------------------------------------------------------------
// Section: Forward declarations
// ----------------------------------------------------------------------------

class Effect;
class EffectInfo;


// ----------------------------------------------------------------------------
// Section: EffectFactory class
// ----------------------------------------------------------------------------

/**
 * Produces concrete effects of the Effect class hierarchy from generic
 * effect descriptor objects of class EffectInfo. 
 */
class EffectFactory
{
public:
	Effect * createEffect(const EffectInfo * pEffectInfo) const;
};

}

#endif /*EFFECTFACTORY_H_*/
