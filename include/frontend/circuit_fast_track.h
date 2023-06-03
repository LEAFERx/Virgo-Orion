#pragma once

#include <utility>
#include <unordered_map>
#include <vector>
#include "common/prime_field.h"

namespace frontend
{

class gate
{
public:
	int ty;
	long long u, v;

	// For weight summation gate (ty=14) only
	long long wsum_k; // number of gates to be summed
	std::vector<long long> wsum_gates;
	std::vector<common::prime_field::field_element> wsum_weights;

	gate(){}
	gate(int t, long long U, long long V)
	{
		ty = t, u = U, v = V;
	}
	gate(int t, long long K, std::vector<long long> _wsum_gates, std::vector<common::prime_field::field_element> _wsum_weights)
	{
		assert(ty == 15);
		assert(_wsum_gates.size() == wsum_k);
		assert(_wsum_weights.size() == wsum_k);
		ty = t;
		u = v = 0;
		wsum_k = K;
		wsum_gates = _wsum_gates;
		wsum_weights = _wsum_weights;
	}
};

class layer
{
public:
	gate *gates;
	int bit_length;
	std::unordered_map<int, std::vector<std::pair<int, std::pair<int, int> > > > u_gates;
	std::unordered_map<int, std::vector<std::pair<int, std::pair<int, int> > > > v_gates;
	bool is_parallel;
	int block_size;
	int log_block_size;
	int repeat_num;
	int log_repeat_num;
	layer()
	{
		is_parallel = false;
		gates = NULL;
		bit_length = 0;
	}
	~layer()
	{
	}
};

class layered_circuit
{
public:
	layer *circuit;
	int total_depth;
	layered_circuit()
	{
		circuit = NULL;
		total_depth = 0;
	}
	~layered_circuit()
	{
		
	}
};

} // namespace frontend

