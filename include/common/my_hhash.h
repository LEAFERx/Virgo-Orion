#pragma once

/**extra 'h' before hhash can avoid some strange error by the compiler*/

#include <immintrin.h>
#include <wmmintrin.h>
#include <cassert>
#include "flo-shani-aesni/sha256/flo-shani.h"

#define USESHA3

extern "C"{
    #include "libXKCP.a.headers/SimpleFIPS202.h"
}

namespace common {

class __hhash_digest
{
public:
    __m128i h0, h1;
};

inline bool equals(const __hhash_digest &a, const __hhash_digest &b)
{
    __m128i v0 = _mm_xor_si128(a.h0, b.h0);
    __m128i v1 = _mm_xor_si128(a.h1, b.h1);
    return _mm_test_all_zeros(v0, v0) && _mm_test_all_zeros(v1, v1);
}

inline void my_hhash(const void* src, void* dst)
{
#ifdef USESHA3
    SHA3_256((unsigned char*)dst, (const unsigned char*)src, 64);
#else
    sha256_update_shani((const unsigned char*)src, 64, (unsigned char*)dst);
#endif
}

} // namespace common
