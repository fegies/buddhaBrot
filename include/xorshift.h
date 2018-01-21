#pragma once

#include <stdint.h>

#define XORSHIFT_MAX 0xFFFFFFFFFFFFFFFF


uint64_t xorshift128plus(uint64_t xorshift_state[2]);
