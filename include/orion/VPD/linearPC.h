#pragma once

#include "common/prime_field.h"
#include "orion/infrastructure/merkle_tree.h"

namespace orion {

//commit
common::__hhash_digest* commit(const common::prime_field::field_element *src, long long N);
//open

//verify
::std::pair<common::prime_field::field_element, bool> open_and_verify(common::prime_field::field_element x, long long N, common::__hhash_digest *com_mt);
::std::pair<common::prime_field::field_element, bool> open_and_verify(common::prime_field::field_element *r, int size_r, int N, common::__hhash_digest *com_mt);

} // namespace orion
