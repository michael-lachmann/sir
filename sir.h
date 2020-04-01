// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
// code for SIR on networks by Petter Holme (2018)

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#define _USE_MATH_DEFINES
#include <math.h>
#include <float.h>
#include <fcntl.h>
#ifdef __linux__
#include <unistd.h>
#elif _WIN32
#include <io.h>
#pragma warning(disable : 4996)
#else
#endif

#include <limits.h>
#include <stdint.h>
#ifdef TIME
#include <time.h>
#endif

typedef double real ;

#define NAVG 10 // number of runs for averages

#define I_OR_R (UINT_MAX - 1)
#define NONE UINT_MAX

#define S(x) (n[(x)].time > now)

// auxiliary macro
#define SQ(x) ((x) * (x))

typedef struct GLOBALS {
	// INPUT PARAMETERS
	real beta, INV_beta, gamma, INV_gamma; // infection rate
	real Em, Ev ; // mean and variance of time in E

	real *w_time, *w_val ;
	int w_n ;
	// NETWORK SPECS
	unsigned int n;
	// OTHER GLOBALS
	unsigned int nheap, *heap;
	real *weight  ; // weight for edge indexes - so th
	unsigned int nweight ;
	// OUTBREAK STATS
	unsigned int s[4]; // number of S E I R.
	real t, now;
	// FOR RNG
	uint64_t state;
	uint32_t rmem;
	real rexp[0x10000];
	real rnorm_E[0x10000]; // distribution of times in E.
} GLOBALS;

enum infection_state
{
	s_S, s_E, s_I, s_R 
};


typedef struct NODE {
	unsigned int deg; // degree (num of nb and w)
	unsigned int *nb ; // neighbors
	unsigned int *w;  // edge weight index

	enum infection_state state ; 

	unsigned int heap; // place on heap
	real time;        // time of next event
} NODE;


// heap.c
extern void up_heap (unsigned int);
extern void down_heap( unsigned int) ;
extern void del_root ();

// misc.c
extern void init_rng ();
extern void read_data (FILE *, unsigned int w);
extern void read_data3 (FILE *);
void gen_full_nwk( unsigned int n, unsigned int clust, unsigned int w) ;
void gen_nwk( unsigned int n, unsigned int deg, unsigned int w) ;

// pcg_rnd.c
extern uint16_t pcg_16 ();
extern uint32_t pcg_32 ();
extern uint32_t pcg_32_bounded ();
extern void pcg_init ();

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
