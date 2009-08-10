/******************************************************************************
 * SIENA: Simulation Investigation for Empirical Network Analysis
 * 
 * Web: http://www.stats.ox.ac.uk/~snijders/siena/
 * 
 * File: ConfigurationTable.h
 * 
 * Description: This file defines the class ConfigurationTable.
 *****************************************************************************/

#ifndef CONFIGURATIONTABLE_H_
#define CONFIGURATIONTABLE_H_

namespace siena
{

// ----------------------------------------------------------------------------
// Section: Forward declarations
// ----------------------------------------------------------------------------

class NetworkVariable;


// ----------------------------------------------------------------------------
// Section: Class description
// ----------------------------------------------------------------------------

/**
 * This class defines a table storing the number of network configurations
 * depending on pairs of actors. The first actor is common to all pairs and
 * is stored indirectly in the owner of this table, which is an instance of
 * the NetworkVariable.
 * 
 * There can be various types of configurations (like the number of two-paths,
 * etc.), each implemented in a derived class. Normally, the derived classes
 * should implement only the purely virtual method vCalculate.
 * 
 * The configuration tables are used to speedup the calculations involving
 * effects.
 */
class ConfigurationTable
{
public:
	ConfigurationTable(NetworkVariable * pVariable);
	virtual ~ConfigurationTable();

	void calculate();	
	inline int get(int i) const;
	void invalidate();
	
protected:
	NetworkVariable * pVariable() const;
	
	/**
	 * An abstract method that has to be implemented by derived classes
	 * to actually calculate the number of configurations of the respective
	 * type.
	 */
	virtual void vCalculate() = 0;
	
	void set(int i, int value);
	void reset();
	
private:
	// The owner dependent variable of this table
	NetworkVariable * lpVariable;
	
	// The internal storage
	int * lpTable;
	
	// Indicates if the table is valid or should be recalculated
	bool lvalid;
};


// ----------------------------------------------------------------------------
// Section: Inline methods
// ----------------------------------------------------------------------------

/**
 * Returns the number of configurations corresponding to the given actor.
 */
inline int ConfigurationTable::get(int i) const
{
	return this->lpTable[i];
}

}

#endif /*CONFIGURATIONTABLE_H_*/
