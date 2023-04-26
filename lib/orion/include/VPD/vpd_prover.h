#ifndef __orion_vpd_prover
#define __orion_vpd_prover
#include <vector>
#include "lib/orion/include/linear_gkr/prime_field.h"
#include "lib/orion/include/infrastructure/my_hhash.h"

namespace orion {

__hhash_digest vpd_prover_init(prime_field::field_element *l_eval, prime_field::field_element *&l_coef, int log_input_length, int slice_size, int slice_count);

} // namespace orion
#endif