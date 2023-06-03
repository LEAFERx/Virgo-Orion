#pragma once

#include <cassert>
#include "common/constants.h"
#include "common/my_hhash.h"
#include "common/prime_field.h"
#include "orion/infrastructure/merkle_tree.h"

namespace orion {

void init_scratch_pad(long long order);

extern common::prime_field::field_element* twiddle_factor;

//Given coefficient of polynomials, output the evaluations on the group generated by root of unity(RoT) [RoT^i for i in range(order)]
void fast_fourier_transform(const common::prime_field::field_element *coefficients, int coef_len, int order, common::prime_field::field_element root_of_unity, common::prime_field::field_element *result, common::prime_field::field_element *twiddle_fac = twiddle_factor);


//see https://codeforces.com/blog/entry/48798 for tutorial
void inverse_fast_fourier_transform(common::prime_field::field_element *evaluations, int coef_len, int order, common::prime_field::field_element root_of_unity, common::prime_field::field_element *result);

}
