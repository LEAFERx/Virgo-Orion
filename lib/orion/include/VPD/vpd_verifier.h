#ifndef __orion_vpd_verifier
#define __orion_vpd_verifier
#include <vector>
#include "lib/orion/include/linear_gkr/prime_field.h"

namespace orion {

bool vpd_verify(prime_field::field_element all_mask_sum, double &v_time);

} // namespace orion

#endif