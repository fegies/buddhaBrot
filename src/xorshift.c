//
// Created by felix on 1/21/18.
//

#include "xorshift.h"

uint64_t xorshift128plus(uint64_t xorshift_state[2])
{
	uint64_t x = xorshift_state[0];
	uint64_t const y = xorshift_state[1];
	xorshift_state[0] = y;
	x ^= x << 23; // a
	xorshift_state[1] = x ^ y ^ (x >> 17) ^ (y >> 26); // b, c
	return xorshift_state[1] + y;
}
