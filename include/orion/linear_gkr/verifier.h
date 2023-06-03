#pragma once

#include <utility>
#include "common/polynomial.h"
#include "orion/linear_gkr/circuit_fast_track.h"
#include "orion/linear_gkr/prover.h"
#include "orion/poly_commitment/poly_commit.h"

namespace orion {

enum gate_types
{
	add = 0,
	mult = 1,
	dummy = 2,
	sum = 5,
	exp_sum = 12,
	direct_relay = 4,
	not_gate = 6,
	minus = 7,
	xor_gate = 8,
	bit_test = 13,
	relay = 10,
	custom_linear_comb = 14,
	input = 3
};

class zk_verifier
{
public:
	int proof_size;
	double v_time;
	poly_commit::poly_commit_verifier poly_ver;
	/** @name Randomness&Const 
	* Storing randomness or constant for simplifying computation/
	*/
	///@{
	common::prime_field::field_element *beta_g_r0_first_half, *beta_g_r0_second_half;
	common::prime_field::field_element *beta_g_r1_first_half, *beta_g_r1_second_half;
	common::prime_field::field_element *beta_u_first_half, *beta_u_second_half;
	common::prime_field::field_element *beta_v_first_half, *beta_v_second_half;

	common::prime_field::field_element *beta_g_r0_block_first_half, *beta_g_r0_block_second_half;
	common::prime_field::field_element *beta_g_r1_block_first_half, *beta_g_r1_block_second_half;
	common::prime_field::field_element *beta_u_block_first_half, *beta_u_block_second_half;
	common::prime_field::field_element *beta_v_block_first_half, *beta_v_block_second_half;
	///@}
	layered_circuit C; //!< The circuit
	zk_prover *p; //!< The prover
	void beta_init(int depth, common::prime_field::field_element alpha, common::prime_field::field_element beta,
	const common::prime_field::field_element* r_0, const common::prime_field::field_element* r_1, 
	const common::prime_field::field_element* r_u, const common::prime_field::field_element* r_v, 
	const common::prime_field::field_element* one_minus_r_0, const common::prime_field::field_element* one_minus_r_1, 
	const common::prime_field::field_element* one_minus_r_u, const common::prime_field::field_element* one_minus_r_v);
	void read_circuit(const char *, const char*);
	void read_r1cs(const char *, const char*, const char*, const char*, const char*);
	bool verify(const char*);
	void get_prover(zk_prover*);
	void delete_self();
	void init_array(int);
	::std::vector<common::prime_field::field_element> predicates(int depth, common::prime_field::field_element *r_0, common::prime_field::field_element *r_1, common::prime_field::field_element *r_u, common::prime_field::field_element *r_v, common::prime_field::field_element alpha, common::prime_field::field_element beta);
	common::prime_field::field_element direct_relay(int depth, common::prime_field::field_element *r_g, common::prime_field::field_element *r_u);
	common::prime_field::field_element V_in(const common::prime_field::field_element*, const common::prime_field::field_element*, common::prime_field::field_element*, int, int);
	common::prime_field::field_element *VPD_randomness, *one_minus_VPD_randomness;
	void self_inner_product_test(common::prime_field::field_element alpha_beta_sum, common::prime_field::field_element v_in);
	/**Test the evaluation of all mask polys after doing random linear combination for them. */
	bool verify_poly_commitment(common::prime_field::field_element* all_sum, double &v_time, int &proof_size, double &p_time, common::__hhash_digest merkle_tree_l, common::__hhash_digest merkle_tree_h);
};

} // namespace orion
