#include "orion/linear_code/linear_code_encode.h"

namespace orion {

common::prime_field::field_element *scratch[2][100];
bool __encode_initialized = false;

} // namespace orion