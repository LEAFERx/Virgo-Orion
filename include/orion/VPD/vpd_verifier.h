#pragma once

#include <vector>
#include "common/prime_field.h"

namespace orion {

bool vpd_verify(common::prime_field::field_element all_mask_sum, double &v_time);

} // namespace orion
