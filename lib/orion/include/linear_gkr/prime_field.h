#ifndef __orion_prime_field
#define __orion_prime_field

//#include <boost/multiprecision/cpp_int.hpp>
//#include <boost/random.hpp>
// #include "lib/orion/include/infrastructure/constants.h"
// #include "lib/orion/include/infrastructure/utility.h"
// #include <cassert>
// #include <immintrin.h>
// #include <vector>
// #include <memory>
// #include <cstring>

#include <cassert>
#include <string>
#include <gmp.h>
#include <gmpxx.h>
#include <vector>

//using namespace boost::multiprecision;
//using namespace boost::random;

namespace orion {

namespace prime_field
{
	const __uint128_t __max_ull = ULLONG_MAX;
	
	class u256b
	{
	public:
		unsigned long long lo, mid;
		__uint128_t hi;
		u256b();
		u256b(const unsigned long long &x);
		u256b(const __uint128_t &x);
 
		u256b(const char* x, int len, int base);
 
		inline u256b operator + (const u256b &x) const {
			u256b ret;
			bool carry;
			__uint128_t midd;
			ret.lo = lo + x.lo;
			carry = (ret.lo < lo);
			midd = (__uint128_t)mid + (__uint128_t)x.mid;
			if(midd == __max_ull)
				midd += carry;
			ret.mid = mid + x.mid + carry;
			ret.hi = hi + x.hi + (midd >> 64);
			return ret;
		}
 
		inline u256b operator - (const u256b &x) const {
			u256b not_x;
			static u256b one(1ull);
			not_x.lo = ~x.lo;
			not_x.mid = ~x.mid;
			not_x.hi = ~x.hi;
			return *this + not_x + one;
		}
 
		//inline u256b operator << (const int &x) const;
		//inline u256b operator >> (const int &x) const;
 
 
		u256b operator * (const u256b &x) const;
 
		inline u256b left_128() {
			u256b ret;
			ret.lo = 0;
			ret.mid = 0;
			ret.hi = ((__uint128_t)lo | ((__uint128_t)mid << 64));
			return ret;
		}
 
		inline u256b operator & (const u256b &x) const {
			u256b ret;
			ret.lo = lo & x.lo;
			ret.mid = mid & x.mid;
			ret.hi = hi & x.hi;
			return ret;
		}

		inline int bit_at(int pos) const {
			if(pos < 64)
				return (lo >> pos) & 1;
			if(pos < 128)
				return (mid >> (pos - 64)) & 1;
			else
				return (hi >> (pos - 128)) & 1;
		}

		inline bool operator <= (const u256b &x) const {
			if(hi < x.hi)
				return true;
			if(hi > x.hi)
				return false;
			if(mid < x.mid)
				return true;
			if(mid > x.mid)
				return false;
			if(lo <= x.lo)
				return true;
			return false;
		}

		inline bool operator != (const u256b &x) const {
			return hi != x.hi || lo != x.lo || mid != x.mid;
		}
 
		inline bool operator > (const u256b &x) const {
			if(hi > x.hi)
				return true;
			if(hi < x.hi)
				return false;
			if(mid > x.mid)
				return true;
			if(mid < x.mid)
				return false;
			return lo > x.lo;
		}

		inline bool operator < (const u256b &x) const {
			if(hi < x.hi)
				return true;
			if(hi > x.hi)
				return false;
			if(mid < x.mid)
				return true;
			if(mid > x.mid)
				return false;
			return lo < x.lo;
		}

		inline u256b midhi_mul(const u256b &a, const u256b &b) {
			u256b ret;
			__uint128_t lolo = (__uint128_t)a.lo * (__uint128_t)b.lo;
	
			__uint128_t lomid1 = (__uint128_t)a.mid * (__uint128_t)b.lo;
			__uint128_t lomid2 = (__uint128_t)a.lo * (__uint128_t)b.mid;
	
			//hi * hi omitted
			ret.lo = (unsigned long long)lolo;
			__uint128_t carry = lolo >> 64; //this carry is less than 2**64
			ret.mid = ((unsigned long long)lomid1 + (unsigned long long)lomid2 + (unsigned long long)carry);				   
			return ret;
		}

 
		inline int bitLen() {
			if(hi != 0)
			{
				unsigned long long hihi, hilo;
				unsigned long long zero = 0;
				hihi = hi >> 64;
				hilo = hi & (~zero);
				if(hihi != 0)
					return 256 - __builtin_clzll(hihi);
				else
					return 256 - 64 - __builtin_clzll(hilo);
			}
			else if(mid != 0)
			{
				return 128 - __builtin_clzll(mid);
			}
			else
			{
				return 64 - __builtin_clzll(lo);
			}
		}
		
		int testBit(int i);
	};
	//extern int512_t mod;
	extern bool initialized;
	//extern independent_bits_engine<mt19937, 256, cpp_int> gen;
 
	extern u256b mod;
 
	class u512b
	{
	public:
		__uint128_t lo, mid;
		u256b hi;
 
		u512b(const u256b &x);
		u512b(const __uint128_t &x);
		u512b(const char* x, int len, int base);
 
		u512b();
		u512b operator + (const u512b &x) const;
 
		u512b operator - (const u512b &x) const;
 
		u512b operator * (const u512b &x) const;
 
		u256b operator % (const u256b &x) const;
		bool operator != (const u512b &x) const;
		bool operator > (const u512b &x) const;
 
		bool operator >= (const u512b &x) const;
 
		bool operator < (const u512b &x) const;
		void random();
		u512b barrett_reduction();
	};
	extern u512b mod_512, minus_mod_512;
 
	void init(std::string, int);
	void init_random();
	u256b convert_in(u512b x);
	

	/*
	This defines a prime field
	*/
	class field_element
	{
	private:
	public:
		u256b value;
 
		field_element();
		field_element(const int);
		field_element(const unsigned long long);
		field_element(const char*, int, int);
		field_element operator + (const field_element &b) const;
		field_element operator * (const field_element &b) const;
		field_element operator / (const field_element &b) const;
		field_element operator - (const field_element &b) const;
		bool operator == (const field_element &b) const;
		field_element mul_non_mod(const field_element &b) const;
		char* to_string();
		int bitLen();
		bool operator != (const field_element &b) const;
		mpz_class to_gmp_class();
		std::vector<bool> bit_stream();
		std::vector<unsigned char> serialize();
	};
	field_element convert_out(field_element x);
	field_element random();
    field_element fast_pow(field_element, long long);
    field_element fast_pow(field_element, u256b);
    field_element inv(field_element);
    field_element get_root_of_unity(long long);
}

} //namespace orion

#endif