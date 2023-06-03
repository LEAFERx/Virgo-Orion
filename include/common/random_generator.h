#pragma once

#include <cstdlib>
#include <ctime>
#include <cstring>
#include <cassert>
/*
This is a implementation of random oracle
*/

namespace common
{

class random_oracle
{
private:
public:
	random_oracle(){}
	~random_oracle(){}
	unsigned long long query(unsigned long long x)
	{
		//most simple impl, debug purpose only
		//it should be a hhash function
		srand(x);
		return rand();
	}
};

} // namespace common

