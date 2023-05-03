#define __debug__
#define __timer__
#include "linear_gkr/verifier_fast_track.h"
#include "linear_gkr/prover_fast_track.h"
#include "linear_gkr/prime_field.h"

verifier v;
prover p;

int main(int argc, char** argv)
{
	prime_field::init("21888242871839275222246405745257275088548364400416034343698204186575808495617", 10);
	p.total_time = 0;
	v.get_prover(&p);
	assert(argc == 3);
	v.read_circuit(argv[1], argv[2]);
	p.get_circuit(v.C);
	bool result = v.verify();
	printf("%s\n", result ? "Pass" : "Fail");
	return 0;
}
