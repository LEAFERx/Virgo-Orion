#pragma once

#include <vector>
#include "common/prime_field.h"
#include "common/my_hhash.h"

namespace orion {

common::__hhash_digest vpd_prover_init(common::prime_field::field_element *l_eval, common::prime_field::field_element *&l_coef, int log_input_length, int slice_size, int slice_count);

} // namespace orion
