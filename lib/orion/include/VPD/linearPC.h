#pragma once
#include "lib/orion/include/linear_gkr/prime_field.h"
#include "lib/orion/include/infrastructure/merkle_tree.h"

namespace orion {

//commit
__hhash_digest* commit(const prime_field::field_element *src, long long N);
//open

//verify
::std::pair<prime_field::field_element, bool> open_and_verify(prime_field::field_element x, long long N, __hhash_digest *com_mt);
::std::pair<prime_field::field_element, bool> open_and_verify(prime_field::field_element *r, int size_r, int N, __hhash_digest *com_mt);

} // namespace orion
