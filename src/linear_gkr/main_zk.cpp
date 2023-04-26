#define __debug__
#define __timer__
#include "linear_gkr/zk_verifier.h"
#include "linear_gkr/zk_prover.h"
#include "linear_gkr/prime_field.h"
#include "lib/orion/include/VPD/linearPC.h"
#include "lib/orion/include/linear_gkr/prime_field.h"
#include "lib/orion/include/linear_code/linear_code_encode.h"
#include <iostream>
#include <cassert>
zk_verifier v;
zk_prover p;

int main(int argc, char** argv)
{
	//std::cout << "hello world" << std::endl;
	// prime_field::init();

	// Orion test
	orion::prime_field::init();
	
	const int N = 1 << 8;
    orion::expander_init(N / orion::column_size);

    orion::prime_field::field_element *coefs = new orion::prime_field::field_element[N];
    for(int i = 0; i < N; ++i)
        coefs[i] = orion::prime_field::random();
    
    auto h = orion::commit(coefs, N);
	auto result = orion::open_and_verify(orion::prime_field::random(), N, h);
    printf("%s\n", result.second ? "succ" : "fail");

	// p.total_time = 0;
	// v.get_prover(&p);
	// //std::cout << "come in" << std::endl;
	// assert(argc == 4);
	// v.read_circuit(argv[1], argv[2]);
	// //std::cout << "after readfile" << std::endl;
	// p.get_circuit(v.C);
	// bool result = v.verify(argv[3]);
	// printf("%s\n", result ? "Pass" : "Fail");
	return 0;
}