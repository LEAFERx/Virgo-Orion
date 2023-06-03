#pragma once

#include <cstring>
#include <utility>
#include <vector>
#include <chrono>
#include "common/prime_field.h"
#include "common/polynomial.h"
#include "common/my_hhash.h"
#include "frontend/circuit_fast_track.h"

namespace frontend
{

class zk_prover
{
public:
	/** @name Basic
	* Basic information and variables about the arithmetic circuit C.*/
	///@{
	common::prime_field::field_element v_u, v_v; //!< two random gates v_u and v_v queried by V in each layer
	int total_uv;
	layered_circuit C;
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
	void get_circuit(const layered_circuit &from_verifier);
	/**Evaluate the output of the circuit.*/
	common::prime_field::field_element* evaluate();
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
	common::quintuple_poly sumcheck_phase1_updatelastbit(common::prime_field::field_element, int);
	common::quadratic_poly sumcheck_phase2_update(common::prime_field::field_element, int);
	common::quintuple_poly sumcheck_phase2_updatelastbit(common::prime_field::field_element, int);

	common::quadratic_poly sumcheck_finalround(common::prime_field::field_element, int, common::prime_field::field_element);
	///@}
	/**I do not know what it is*/
	common::prime_field::field_element V_res(const common::prime_field::field_element*, const common::prime_field::field_element*, const common::prime_field::field_element*, int, int);
	std::pair<common::prime_field::field_element, common::prime_field::field_element> sumcheck_finalize(common::prime_field::field_element);
	void delete_self();
	zk_prover()
	{
		memset(circuit_value, 0, sizeof circuit_value);
	}
	~zk_prover(); 

	/** @name Masks
	* generate mask polynomials of each layer and input layer
	*/
	///@{
	// void generate_maskR(int layer_id);
	// std::vector<common::prime_field::field_element> all_pri_mask;
	// common::prime_field::field_element *maskpoly; 
	// common::prime_field::field_element maskpoly_sumc;
	// common::prime_field::field_element maskpoly_sumr;
	// common::prime_field::field_element rho;
	// void generate_maskpoly_pre_rho(int, int);
	// void generate_maskpoly_after_rho(int, int);
	// common::prime_field::field_element maskR[6], preR[6];

	// common::prime_field::field_element maskR_sumcu, maskR_sumcv, preZu, preZv, Zu, Zv, preu1, prev1, Iuv, prepreu1, preprev1;
	// common::quadratic_poly Rg1, Rg2, sumRc;
	///@}

	/** @name Query
	* Answer the queries of V for mask polynomials
	*/
	///@{
	// common::prime_field::field_element query(common::prime_field::field_element*, common::prime_field::field_element*, common::prime_field::field_element);
	// common::prime_field::field_element queryRg1(common::prime_field::field_element);
	// common::prime_field::field_element queryRg2(common::prime_field::field_element);

	///@}

	/** @name Commit
	* Some function for committing the mask polys before the protocol. */
	///@{
	void sumcheck_maskpoly_init();
	std::vector<common::prime_field::field_element> maskr;
	std::vector<common::prime_field::field_element> commit_input();
	std::vector<common::prime_field::field_element> commit_maskpoly();
	common::__hhash_digest prover_vpd_prepare();
	common::__hhash_digest prover_vpd_prepare_post_gkr(std::vector<common::prime_field::field_element> &all_pub_msk, common::prime_field::field_element *r_0, common::prime_field::field_element *one_minus_r_0, int r_0_len, common::prime_field::field_element target_sum, common::prime_field::field_element *all_sum);
	///@}
};

} // namespace frontend

