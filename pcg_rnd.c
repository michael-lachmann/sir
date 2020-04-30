// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
// code for temporal network SIR by Petter Holme (2018)

// this file contains the random number generator, derived from the PCG
// RNG v0.94 http://www.pcg-random.org under the Apache License 2.0
// http://www.apache.org/licenses/LICENSE-2.0

// 32-bit Output, 64-bit State: PCG-XSH-RS

#include "sir.h"

extern GLOBALS g;

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

void pcg_init ( uint64_t initstate ) {
	int fd;
	unsigned long sz;

	g.state = (initstate + 1442695040888963407ULL) * 6364136223846793005ULL + 1442695040888963407ULL;
}

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

uint32_t pcg_32 () {
	uint64_t state = g.state;
	uint32_t value, rot;

	g.state = g.state * 6364136223846793005ULL + 1442695040888963407ULL;
	value = (uint32_t)(((state >> 18u) ^ state) >> 27u);
	rot = state >> 59u;
	return (value >> rot) | (value << ((- rot) & 31));
}

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

uint32_t pcg_32_bounded (uint32_t bound) {
	uint32_t threshold = -bound % bound, r;

	for ( ; ; ) {
		r = pcg_32();
		if (r >= threshold) return r % bound;
	}
}

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

uint16_t pcg_16 () {
	static unsigned int exist;

	if (exist) {
		exist = 0;
		return g.rmem >> 16;
	}
	exist = 1;
	g.rmem = pcg_32();
	return (uint16_t) g.rmem;
}

uint8_t pcg_8() {
	static unsigned int n=0 ;
	static uint32_t r ;

	if (n == 0) {
		r = pcg_32() ;
		n = 4 ;
	}
	else {
		r /= 0x100 ;
		n-- ;
	}
	return (uint8_t) (r & 0xff);

}

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
