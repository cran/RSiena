#ifndef SQRTTABLE_H_
#define SQRTTABLE_H_

namespace siena
{

class SqrtTable
{
public:
	static SqrtTable * instance();
	double sqrt(int i);
	
private:
	SqrtTable();
	virtual ~SqrtTable();
	
	static SqrtTable * lpInstance;
	double * ltable;
};

}

#endif /*SQRTTABLE_H_*/
