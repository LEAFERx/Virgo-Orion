#define __debug__
#define __timer__
#include <iostream>
#include <cassert>
#include "common/prime_field.h"
#include "frontend/zk_verifier.h"
#include "frontend/zk_prover.h"
#include "orion/linear_code/linear_code_encode.h"

frontend::zk_verifier v;
frontend::zk_prover p;

int main(int argc, char** argv)
{
	std::cout << std::endl << ">>>>>>>>>>>>>>>>>>>>>>>>>>>>> New Session <<<<<<<<<<<<<<<<<<<<<<<<<<<<<" << std::endl << std::endl;
	common::prime_field::init("21888242871839275222246405745257275088548364400416034343698204186575808495617", 10);

	p.total_time = 0;
	v.get_prover(&p);
	//std::cout << "come in" << std::endl;
	assert(argc == 4);
	v.read_circuit(argv[1], argv[2]);
	//std::cout << "after readfile" << std::endl;
	p.get_circuit(v.C);
	bool result = v.verify(argv[3]);
	printf("%s\n", result ? "Pass" : "Fail");
	return 0;
}