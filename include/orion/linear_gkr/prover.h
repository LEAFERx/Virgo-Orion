#pragma once

#include <cstring>
#include <utility>
#include <vector>
#include <chrono>
#include "common/prime_field.h"
#include "common/polynomial.h"
#include "common/my_hhash.h"
#include "orion/linear_gkr/circuit_fast_track.h"
#include "orion/VPD/vpd_prover.h"
#include "orion/VPD/fri.h"
#include "orion/poly_commitment/poly_commit.h"

namespace orion {

class zk_prover
{
public:
	poly_commit::poly_commit_prover poly_prover;
	/** @name Basic
	* Basic information and variables about the arithmetic circuit C.*/
	///@{
	common::prime_field::field_element v_u, v_v; //!< two random gates v_u and v_v queried by V in each layer
	int total_uv;
	layered_circuit *C;
	common::prime_field::field_element* circuit_value[1000000];

	int sumcheck_layer_id, length_g, length_u, length_v;
	///@}

	/** @name Randomness
	* Some randomness or values during the proof phase. */
	///@{
	common::prime_field::field_element alpha, beta;
	const common::prime_field::field_element *r_0, *r_1;
	common::prime_field::field_element *one_minus_r_0, *one_minus_r_1;

	common::linear_poly *addV_array;
	common::linear_poly *V_mult_add;
	common::prime_field::field_element *beta_g_r0_fhalf, *beta_g_r0_shalf, *beta_g_r1_fhalf, *beta_g_r1_shalf, *beta_u_fhalf, *beta_u_shalf;
	common::prime_field::field_element *beta_u, *beta_v, *beta_g;
	common::linear_poly *add_mult_sum;
	///@}


	double total_time;
	/**Initialize some arrays used in the protocol.*/
	void init_array(int);
	/**Read the circuit from the input file.*/
	void get_circuit(layered_circuit &from_verifier);
	/**Evaluate the output of the circuit.*/
	common::prime_field::field_element* evaluate();
	void get_witness(common::prime_field::field_element*, int);
	/**Generate random mask polynomials and initialize parameters. */
	void proof_init();
	
	/** @name Group function for interior sumcheck protocol. 
	*Run linear GKR protocol: for each layer of the circuit, our protocol is based on sumcheck protocol. And there are total three phases of the sumcheck:
	*1. sumcheck for gate u
	*2. sumcheck for gate v
	*3. final round for the extra mask bit
 	*/ 
 	///@{
	void sumcheck_init(int layer_id, int bit_length_g, int bit_length_u, int bit_length_v, const common::prime_field::field_element &,
		const common::prime_field::field_element &, const common::prime_field::field_element*, const common::prime_field::field_element*, common::prime_field::field_element*, common::prime_field::field_element*);
	void sumcheck_phase1_init();
	void sumcheck_phase2_init(common::prime_field::field_element, const common::prime_field::field_element*, const common::prime_field::field_element*);
	common::quadratic_poly sumcheck_phase1_update(common::prime_field::field_element, int);
	common::quadratic_poly sumcheck_phase2_update(common::prime_field::field_element, int);

	///@}
	/**I do not know what it is*/
	common::prime_field::field_element V_res(const common::prime_field::field_element*, const common::prime_field::field_element*, const common::prime_field::field_element*, int, int);
	::std::pair<common::prime_field::field_element, common::prime_field::field_element> sumcheck_finalize(common::prime_field::field_element);
	void delete_self();
	zk_prover()
	{
		memset(circuit_value, 0, sizeof circuit_value);
	}
	~zk_prover(); 
	common::__hhash_digest prover_vpd_prepare();
	common::__hhash_digest prover_vpd_prepare_post_gkr(common::prime_field::field_element *r_0, common::prime_field::field_element *one_minus_r_0, int r_0_len, common::prime_field::field_element target_sum, common::prime_field::field_element *all_sum);
	///@}
};

} // namespace orion
