#pragma once

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <immintrin.h>
#include <wmmintrin.h>
#include "common/my_hhash.h"
#include "common/prime_field.h"


namespace orion {

namespace merkle_tree
{
extern int size_after_padding;

common::__hhash_digest hash_single_field_element_zero(common::prime_field::field_element x);

common::__hhash_digest hash_single_field_element(common::prime_field::field_element x, common::__hhash_digest prev_hash);

namespace merkle_tree_prover
{
    //Merkle tree functions used by the prover
    //void create_tree(void* data_source_array, int lengh, void* &target_tree_array, const int single_element_size = 256/8)
    void create_tree(void* data, int ele_num, common::__hhash_digest* &dst, const int element_size, bool alloc_required);
}

namespace merkle_tree_verifier
{
    //Merkle tree functions used by the verifier
    bool verify_claim(common::__hhash_digest root_hhash, const common::__hhash_digest* tree, common::__hhash_digest hhash_element, int pos_element, int N, bool *visited, long long &proof_size);
}
}

} // namespace orion
