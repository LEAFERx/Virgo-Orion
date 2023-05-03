#include "lib/orion/include/linear_gkr/prime_field.h"
#include <cmath>
#include <climits>
#include <cstdio>
#include <ctime>
#include <chrono>
#include <immintrin.h>

namespace orion::prime_field
{
	const __uint128_t __zero = 0LL;
	const __uint128_t __max_128 = ~(__zero);
 
	u256b::u256b(){lo = mid = 0; hi = 0;}
	u256b::u256b(const unsigned long long &x)
	{
		lo = x;
		mid = 0;
		hi = 0;
	}
	u256b::u256b(const __uint128_t &x)
	{
		lo = (unsigned long long)(x & __max_ull);
		mid = (unsigned long long)(x >> 64);
		hi = 0;
	}
 
	u256b::u256b(const char* x, int len, int base)
	{
		u256b ret = (u256b)(unsigned long long)0;
		u256b tmp = (u256b)(unsigned long long)1;
		for(int i = len - 1; i >= 0; --i)
		{
			ret = ret + tmp * ((u256b)(unsigned long long)(x[i] - '0'));
			tmp = tmp * ((u256b)(unsigned long long)base);
		}
		*this = ret;
	}
 
	u256b u256b::operator * (const u256b &x) const
	{
		u256b ret;
		__uint128_t lolo = (__uint128_t)lo * (__uint128_t)x.lo;
 
		__uint128_t lomid1 = (__uint128_t)mid * (__uint128_t)x.lo;
		__uint128_t lomid2 = (__uint128_t)lo * (__uint128_t)x.mid;
 
		//hi * hi omitted
		ret.lo = (unsigned long long)lolo;
		__uint128_t carry = lolo >> 64; //this carry is less than 2**64
		ret.mid = ((unsigned long long)lomid1 + (unsigned long long)lomid2 + (unsigned long long)carry);
		//this carry is not necessary less than 2**64
		__uint128_t carry2 = ((lomid1 >> 64) + (lomid2 >> 64)) + 
						   (((lomid1 & __max_ull) + (lomid2 & __max_ull) + carry) >> 64);
 
		ret.hi = (__uint128_t)lo * (__uint128_t)x.hi + (__uint128_t)hi * (__uint128_t)x.lo + 
				 (__uint128_t)mid * (__uint128_t)x.mid + (((__uint128_t)mid * (__uint128_t)x.hi + (__uint128_t)hi * (__uint128_t)x.mid) << 64) + carry2;
		return ret;
	}
 
	int u256b::testBit(int i)
	{
		if(i <= 64)
		{
			return (lo >> i) & 1;
		}
		else if(i <= 128)
		{
			return (mid >> (i - 64)) & 1;
		}
		else
		{
			return (hi >> (i - 128)) & 1;
		}
	}
 
	u256b one, zero;
 
	u512b::u512b(const u256b &x)
	{
		lo = ((__uint128_t)x.mid << 64) | (__uint128_t)x.lo;
		mid = x.hi;
		hi.lo = 0;
		hi.mid = 0;
		hi.hi = 0;
	}
	u512b::u512b(const __uint128_t &x)
	{
		lo = x;
		mid = 0;
		hi.lo = 0;
		hi.mid = 0;
		hi.hi = 0;
	}
	u512b::u512b(const char* x, int len, int base)
	{
		u512b ret = (u256b)(unsigned long long)0;
		u512b tmp = (u256b)(unsigned long long)1;
		for(int i = len - 1; i >= 0; --i)
		{
			ret = ret + tmp * (u512b)((u256b)(unsigned long long)(x[i] - '0'));
			tmp = tmp * (u512b)((u256b)(unsigned long long)10);
		}
		*this = ret;
	}
 
	u512b::u512b(){lo = mid = 0; hi.lo = 0; hi.mid = 0; hi.hi = 0;}
	u512b u512b::operator + (const u512b &x) const
	{
		u512b ret;
		__uint128_t carry, carry2;
		ret.lo = lo + x.lo;
		carry = ret.lo < lo;
		ret.mid = mid + x.mid + carry;
		if(carry == 0)
			carry2 = ret.mid < mid;
		else
			carry2 = ret.mid <= mid;
		ret.hi = hi + x.hi + carry2;
		return ret;
	}
 
	u512b u512b::operator - (const u512b &x) const
	{
		u512b not_x;
		not_x.hi.hi = ~x.hi.hi;
		not_x.hi.mid = ~x.hi.mid;
		not_x.hi.lo = ~x.hi.lo;
		not_x.mid = ~x.mid;
		not_x.lo = ~x.lo;
		not_x = not_x + one;
		return *this + not_x;
	}
 
	inline void mul_no_hi(const u256b &a, const u256b &b, u256b &ret)
	{
		unsigned long long lolo_lo, lolo_hi;
		lolo_lo = _mulx_u64(a.lo, b.lo, &lolo_hi);
		unsigned long long lomid1_lo, lomid1_hi;
		lomid1_lo = _mulx_u64(a.mid, b.lo, &lomid1_hi);
		unsigned long long lomid2_lo, lomid2_hi;
		lomid2_lo = _mulx_u64(a.lo, b.mid, &lomid2_hi);
		//const __uint128_t lolo = (__uint128_t)a.lo * (__uint128_t)b.lo;
 
		//const __uint128_t lomid1 = (__uint128_t)a.mid * (__uint128_t)b.lo;
		//const __uint128_t lomid2 = (__uint128_t)a.lo * (__uint128_t)b.mid;
 
		//hi * hi omitted
		ret.lo = (unsigned long long)lolo_lo;
		const __uint128_t carry = lolo_hi; //this carry is less than 2**64
		ret.mid = ((unsigned long long)lomid1_lo + (unsigned long long)lomid2_lo + (unsigned long long)carry);
		//this carry is not necessary less than 2**64
		const __uint128_t carry2 = ((__uint128_t)(lomid1_hi) + (lomid2_hi)) + 
						   (((__uint128_t)(lomid1_lo) + (lomid2_lo) + carry) >> 64);
 
		ret.hi = (__uint128_t)a.mid * (__uint128_t)b.mid + carry2;
	}
	inline u512b u512b::operator * (const u512b &x) const
	{
		if(x.mid == 0 && x.lo == 0)
			return x;
		if(mid == 0 && lo == 0)
			return *this;
		u512b ret;
		const u256b alo = (u256b)lo, xlo = (u256b)x.lo;
		const u256b amid = (u256b)mid, xmid = (u256b)x.mid;
		u256b lolo; mul_no_hi(alo, xlo, lolo);
 
		u256b lomid1; mul_no_hi(amid, xlo, lomid1);
		u256b lomid2; mul_no_hi(alo, xmid, lomid2);
 
		u256b midmid; mul_no_hi(amid, xmid, midmid);
 
 
		//hi * hi omitted
		ret.lo = (__uint128_t)lolo.lo | ((__uint128_t)lolo.mid << 64);
		const __uint128_t carry = lolo.hi; //this carry is less than 2**128
 
		const u256b tmp = (lomid1 + lomid2 + carry);
 
		ret.mid = (__uint128_t)tmp.lo | ((__uint128_t)tmp.mid << 64);
		//this carry is not necessary less than 2**128
		//double sum = (double)lomid1.hi + (double)lomid2.hi + (double)lomid1.lo + (double)((__uint128_t)lomid1.mid << 64) + (double)lomid2.lo + (double)((__uint128_t)lomid2.mid << 64) + (double)carry;
		//sum = sum / max_int128;
		const u256b carry2 = ((u256b)(lomid1.hi) + (u256b)(lomid2.hi)) + 
						   (((u256b)((__uint128_t)lomid1.lo | ((__uint128_t)lomid1.mid << 64))
						   + (u256b)((__uint128_t)lomid2.lo | ((__uint128_t)lomid2.mid << 64)) + (u256b)carry).hi);
 
		//const u256b carry2 = (unsigned long long)(sum);
 
		//ret.hi = lohi1 + lohi2 + midmid + ((midhi1 + midhi2).left_128()) + carry2;
		ret.hi = midmid + carry2;
		return ret;
	}
	u256b my_factor;
 
//https://www.nayuki.io/page/barrett-reduction-algorithm
 
	u512b u512b::barrett_reduction()
	{
		if(lo == 0 && mid == 0 && hi.lo == 0 && hi.mid == 0 && hi.hi == 0)
			return *this;
		u512b hi_factor = (u512b)hi * (u512b)my_factor;
		u512b lo256bit;
		lo256bit.hi = (unsigned long long)0;
		lo256bit.mid = mid;
		lo256bit.lo = lo;
		const u512b lo_factor = (u512b)lo256bit * (u512b)my_factor;
		const u512b res = hi_factor + (u512b)lo_factor.hi;
		u256b t;
		t.hi = res.hi.hi << 4 | (res.hi.mid >> 60);
		t.mid = res.hi.mid << 4 | (res.hi.lo >> 60);
		t.lo = (res.hi.lo << 4) | (res.mid >> 124);
		u512b t512 = (u512b)t;
		t512 = *this - t512 * mod;
 
		auto result = t512;
		if(t512 >= mod)
			result = t512 - mod;
		return result;
	}
 
 
 
	u256b mod;
	u512b minus_mod_512, mod_512;
	const int shift = 508;
	bool initialized = false;
    field_element root_of_unity;
    const int max_order = 28;
	const int reducer_bit = 256;
	u256b mask;
	u256b reciproal;
	u256b factor;
	u512b converted_one;
	u256b converted_one_256;
	u256b minus_mod;
	bool verbose = false;
 
	u256b u512b::operator % (const u256b &x) const
	{
		//montgomery reduction
		u256b self_lo_256;
		self_lo_256.lo = lo;
		self_lo_256.mid = lo >> 64;
		self_lo_256.hi = mid;
 
		u256b temp = self_lo_256 * factor;
		u256b reduced = (*this + (u512b)temp * mod).hi;
 
		// if(verbose)
		// {
		// 	prime_field::field_element low_256, tempf, reducedf, reducedf_minus_mod;
		// 	low_256.value = self_lo_256;
		// 	tempf.value = temp;
		// 	reducedf.value = reduced;
		// 	reducedf_minus_mod.value = reduced - mod;
 
		// 	std::cout << "self_lo_256: " << low_256.to_gmp_class() << std::endl;
		// 	std::cout << "temp: " << tempf.to_gmp_class() << std::endl;
		// 	std::cout << "reduced: " << reducedf.to_gmp_class() << std::endl;
		// 	std::cout << "mod - reduced " << reducedf_minus_mod.to_gmp_class() << std::endl;
		// }
 
		if(mod <= reduced)
			reduced = reduced - mod;
		return reduced;
	}
	bool u512b::operator != (const u512b &x) const
	{
		return lo != x.lo || mid != x.mid || hi != x.hi;
	}
	bool u512b::operator > (const u512b &x) const
	{
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
 
	bool u512b::operator >= (const u512b &x) const
	{
		if(hi > x.hi)
			return true;
		if(hi < x.hi)
			return false;
		if(mid > x.mid)
			return true;
		if(mid < x.mid)
			return false;
		return lo >= x.lo;
	}
 
	bool u512b::operator < (const u512b &x) const
	{
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
	void u512b::random()
	{
		lo = rand();
		mid = 0;
		hi = u256b(0ull);
	}

	u256b right_shift_1(u256b &x) {
		u256b ret;
		ret.lo = (x.lo >> 1) | (x.mid << 63);
		ret.mid = (x.mid >> 1) | ((x.hi & 1) << 63);
		ret.hi = (x.hi >> 1);
		return ret;
	}
 
	u256b convert_in(u512b x)
	{
		u512b value;
		value.lo = value.mid = 0;
		value.hi.lo = x.lo;
		value.hi.mid = x.lo >> 64;
		value.hi.hi = x.mid;
		value = value.barrett_reduction();
		u256b result;
		result.lo = value.lo;
		result.mid = value.lo >> 64;
		result.hi = value.mid;
		return result;
	}
	void init(std::string s, int base)
	{
		assert(base == 10);
		initialized = true;
		mod = u256b(s.c_str(), s.length(), base);
		one = (u256b)(unsigned long long)1;
		zero = (u256b)(unsigned long long)0;
 
		std::string factor_s = "38284845454613504619394467267190322316714506535725634610690744705837986343205";
		my_factor = u256b(factor_s.c_str(), factor_s.length(), 10);
		minus_mod_512 = prime_field::field_element(0).value - mod;
		minus_mod = prime_field::field_element(0).value - mod;
		mod_512 = mod;
 
        std::string base_root_of_unity = "6506949774839052718110406215085119770091102268120408611511048664532142289545";
        root_of_unity.value = u256b(base_root_of_unity.c_str(), base_root_of_unity.length(), 10);
		root_of_unity.value = convert_in(root_of_unity.value);
 
		//montgomery_reduction_init
 
		/*
			Python code for the montgomery reduction init
			self.mod = mod
 
			self.reducer_bits = 256
			self.reducer = 1 << self.reducer_bits
			self.mask = self.reducer - 1
 
			self.reciproal = fastpow(self.reducer, mod - 2)
			assert(self.reciproal * self.reducer % mod == 1)
			self.factor = (self.reducer * self.reciproal - 1) // mod
			self.conertedOne = self.reducer % mod
		*/
		mask.lo = ~0ull;
		mask.mid = ~0ull;
		mask.hi = 0;
		mask.hi = ~mask.hi;
		std::string reciproal_str = "9915499612839321149637521777990102151350674507940716049588462388200839649614";
		std::string factor_str = "52454480824480482120356829342366457550537710351690908576382634413609933864959";
		std::string converted_one_str = "6350874878119819312338956282401532410528162663560392320966563075034087161851";
		reciproal = u256b(reciproal_str.c_str(), reciproal_str.length(), 10);
		factor = u256b(factor_str.c_str(), factor_str.length(), 10);
		converted_one = u512b(converted_one_str.c_str(), converted_one_str.length(), 10);
		converted_one_256 = u256b(converted_one_str.c_str(), converted_one_str.length(), 10);
 
	}
	void init_random()
	{
	}
	u256b convert_in(unsigned long long x)
	{
		u512b value;
		value.lo = 0;
		value.mid = 0;
		value.hi = x;
		value = value.barrett_reduction();
		u256b result;
		result.lo = value.lo;
		result.mid = value.lo >> 64;
		result.hi = value.mid;
		return result;
	}
	field_element convert_out(field_element x)
	{
		auto product = (u512b)x.value * (u512b)reciproal;
		auto u512bvalue = product.barrett_reduction();
		x.value.hi = u512bvalue.mid;
		x.value.lo = u512bvalue.lo;
		x.value.mid = u512bvalue.lo >> 64;
		return x;
	}
 
	field_element::field_element(const char* str, int len, int base)
	{
		value = u256b(str, len, base);
		//value_original = value;
		value = convert_in(value);
	}
	field_element::field_element()
	{
		value = (u256b)(unsigned long long)0;
		//value_original = (u256b)(unsigned long long)0;
	}
	field_element::field_element(const int x)
	{
		//value_original = (u256b)(unsigned long long)x;
		switch (x)
		{
		case 0:
			value.lo = 0;
			value.mid = 0;
			value.hi = 0;
			break;
		case 1:
			value = converted_one_256;
			break;
		default:
			value = convert_in(x);
			break;
		}
	}
	field_element::field_element(const unsigned long long x)
	{
		switch (x)
		{
		case 0:
			value.lo = 0;
			value.mid = 0;
			value.hi = 0;
			break;
		case 1:
			value = converted_one_256;
			break;
		default:
			value = convert_in(x);
			break;
		}
		//value_original = (u256b)x;
	}
	field_element field_element::operator + (const field_element &b) const
	{
		field_element ret;
		ret.value = (b.value + value);
		if(mod <= ret.value)
			ret.value = ret.value + minus_mod;
		/*
		ret.value_original = b.value_original + value_original;
		if(ret.value_original >= mod_512)
			ret.value_original = ret.value_original + minus_mod_512;
 
		assert(!(convert_out(ret).value != ret.value_original));
		*/
		return ret;
	}
 
	field_element field_element::operator * (const field_element &b) const
	{
		field_element ret;
		ret.value = ((u512b)b.value * (u512b)value) % mod;
		/*
		ret.value_original = (b.value_original * value_original).barrett_reduction();
		if((convert_out(ret).value != ret.value_original))
		{
			prime_field::field_element x, y, a, bb = b, ab;
			ab.value = value * b.value;
			a.value = value;
			x.value = convert_out(ret).value;
			y.value = ret.value_original;
			std::cout << "a: " << a.to_gmp_class() << std::endl;
			std::cout << "b: " << bb.to_gmp_class() << std::endl;
			std::cout << "ab: " << ab.to_gmp_class() << std::endl;
			std::cout << "value: " << x.to_gmp_class() << std::endl;
			std::cout << "value_original: " << y.to_gmp_class() << std::endl;
 
			verbose = true;
			ret.value = (b.value * value) % mod;
		}
		assert(!(convert_out(ret).value != ret.value_original));
		*/
		return ret;
	}
	field_element field_element::operator / (const field_element &b) const
	{
		//todo
		assert(false);
		//return ret;
	}
	field_element field_element::operator - (const field_element &b) const
	{
		field_element ret;
		if(b.value <= value)
			ret.value = value - b.value;
		else
			ret.value = value + mod - b.value;
		/*
		if(value_original >= b.value_original)
			ret.value_original = value_original - b.value_original;
		else
			ret.value_original = value_original + mod_512 - b.value_original;
 
		if(convert_out(ret).value != ret.value_original)
		{
			prime_field::field_element x, y, a, bb = b, va, vb;
			a.value = value;
			x.value = convert_out(ret).value;
			y.value = ret.value_original;
			va.value = value_original;
			vb.value = b.value_original;
			std::cout << "a: " << a.to_gmp_class() << std::endl;
			std::cout << "b: " << bb.to_gmp_class() << std::endl;
			std::cout << "value: " << x.to_gmp_class() << std::endl;
			std::cout << "value_original: " << y.to_gmp_class() << std::endl;
			std::cout << "value_original a: " << va.to_gmp_class() << std::endl;
			std::cout << "value_original b: " << vb.to_gmp_class() << std::endl;
		}
		assert(!(convert_out(ret).value != ret.value_original));
		*/
		return ret;
	}
	char* field_element::to_string()
	{
		char *ret = new char[257];
		for(int i = 0; i < 64; ++i)
			ret[i] = ((value.lo >> i) & 1) + '0';
		for(int i = 0; i < 64; ++i)
			ret[i + 64] = ((value.mid >> i) & 1) + '0';
		for(int i = 0; i < 128; ++i)
			ret[i + 128] = ((value.hi >> i) & 1) + '0';
		ret[256] = '\0';
		return ret;
	}
	field_element random()
	{
		field_element ret;
		ret.value.lo = rand();
		//ret.value_original = ret.value;
		ret.value = convert_in(ret.value);
		return ret;
	}
	bool field_element::operator != (const field_element &b) const
	{
		return value != b.value;
	}
	int field_element::bitLen()
	{
		assert(false);
	}
	mpz_class field_element::to_gmp_class()
	{
		char* str_value;
		str_value = new char[257];
		for(int i = 0; i < 64; ++i)
			str_value[i] = ((value.lo >> i) & 1) + '0';
		for(int i = 0; i < 64; ++i)
			str_value[i + 64] = ((value.mid >> i) & 1) + '0';
		for(int i = 0; i < 128; ++i)
		{
			str_value[i + 128] = ((value.hi >> i) & 1) + '0';
		}
		for(int i = 0; i < 128; ++i)
		{
			char x;
			x = str_value[i];
			str_value[i] = str_value[255 - i];
			str_value[255 - i] = x;
		}
		str_value[256] = 0;
		mpz_class ret(str_value, 2);
		delete str_value;
		return ret;
	} 
	bool field_element::operator == (const field_element &b) const
	{
		return !(*this != b);
	}
	std::vector<bool> field_element::bit_stream()
	{
		std::vector<bool> out;
		for(int i = 0; i < 128; ++i)
		{
			out.push_back((value.lo >> i) & 1);
		}
 
		for(int i = 0; i < 128; ++i)
		{
			out.push_back((value.mid >> i) & 1);
		}
		return out;
	}
	std::vector<unsigned char> field_element::serialize()
	{
		std::vector<unsigned char> out;
		for(int i = 0; i < 128 / 8; ++i)
		{
			out.push_back((value.lo >> (8 * i)) & ((unsigned char)255));
		}
 
		for(int i = 0; i < 128 / 8; ++i)
		{
			out.push_back((value.mid >> (8 * i)) & ((unsigned char)255));
		}
		return out;
	}
    field_element fast_pow(field_element base, long long exp)
    {
        field_element ret;
        ret.value = converted_one_256;
        field_element tmp = base;
        while (exp)
        {
            if (exp & 1)
                ret = ret * tmp;
            exp >>= 1;
            tmp = tmp * tmp;
        }
        return ret;
 
    }
    field_element fast_pow(field_element base, u256b exp)
    {
        field_element ret;
        ret.value = converted_one_256;
        field_element tmp = base;
        while (exp != zero)
        {
            if ((exp & one) != zero)
                ret = ret * tmp;
            // exp >>= 1;
			exp = right_shift_1(exp);
            tmp = tmp * tmp;
        }
        return ret;
 
    }
    field_element inv(field_element x)
    {
        const bool mod_minus_2_bits[] = {1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 1, 1, 1, 1, 1, 0, 0, 1, 0, 0, 1, 1, 0, 1, 0, 1, 1, 1, 1, 1, 0, 0, 0, 0, 1, 1, 1, 1, 1, 0, 0, 0, 0, 1, 0, 1, 0, 0, 0, 1, 0, 0, 1, 0, 0, 0, 0, 1, 1, 1, 0, 1, 0, 0, 1, 1, 1, 0, 1, 1, 0, 0, 1, 1, 1, 1, 0, 0, 0, 0, 1, 0, 0, 1, 0, 0, 0, 0, 1, 0, 1, 1, 1, 1, 1, 0, 0, 1, 1, 0, 0, 0, 0, 0, 1, 0, 1, 0, 0, 1, 0, 1, 1, 1, 0, 1, 0, 0, 0, 0, 1, 1, 0, 1, 0, 1, 0, 0, 0, 0, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 1, 0, 1, 1, 0, 1, 1, 0, 1, 1, 0, 1, 0, 0, 0, 1, 0, 0, 0, 0, 0, 1, 0, 1, 0, 0, 0, 0, 1, 1, 1, 0, 1, 1, 0, 0, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 1, 0, 1, 1, 0, 0, 0, 1, 1, 0, 0, 1, 0, 0, 0, 0, 1, 1, 1, 0, 1, 0, 0, 1, 1, 1, 0, 0, 1, 1, 1, 0, 0, 1, 0, 0, 0, 1, 0, 0, 1, 1, 0, 0, 0, 0, 0, 1, 1};
        field_element ret;
        ret.value = converted_one_256;
		//ret.value_original = 1;
        field_element tmp = x;
        for(int i = 0; i < 254; ++i)
        {
			//printf("%d\n", i);
            if(mod_minus_2_bits[i])
            {
                ret = ret * tmp;
            }
			//std::cout << ret.to_gmp_class() << std::endl;
            tmp = tmp * tmp;
			//std::cout << tmp.to_gmp_class() << std::endl;
		}
 
        return ret;
    }
    field_element get_root_of_unity(long long log_order)
    {
        field_element ret = root_of_unity;
 
        for(int i = 0; i < (max_order - log_order); ++i)
        {
            ret = ret * ret;
        }
		field_element tmp = ret;
		for(int i = 0; i < log_order; ++i)
		{
			tmp = tmp * tmp;
		}
		assert(tmp == prime_field::field_element(1));
        return ret;
    }
}