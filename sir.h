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
typedef char   bool ;
typedef unsigned char state ;
#define TRUE (1==1)
#define FALSE (1!=1)




#define NAVG 10 // number of runs for averages

#define NSTATES (256)
#define NQUART (0x10000)
#define RANDOM_QUART pcg_16
#define EPSILON (1e-10)
#define NEVER (DBL_MAX)
#define MAXRAND (4294967295.0)

#define I_OR_R (UINT_MAX - 1)
#define NONE UINT_MAX

#define S(x) (n[(x)].time > now)

// auxiliary macro
#define SQ(x) ((x) * (x))

typedef struct GLOBALS {
	// INPUT PARAMETERS

	real *w_time, *w_val ;
	int w_n ;
	// NETWORK SPECS
	unsigned int n;
	// OTHER GLOBALS
	unsigned int nheap, *heap;

	real *weight  ; // weight for edge indexes - so th
	unsigned int nweight ;

	unsigned int ngroup ; // number of groups read from net. Can be age groups, for example.
	// OUTBREAK STATS
	unsigned int s[NSTATES]; // number of S E I R.
	real t, now;
	real av_deg ;
	// FOR RNG
	uint64_t state;
	uint32_t rmem;
	real          *time_dist [NSTATES] ;
	state         state_trans[NSTATES][3] ;
	uint32_t	  state_trans_p[NSTATES][3] ;
	real          state_infect_rate[NSTATES] ;
	bool	      self_change[NSTATES] ; // can change without infection?
	bool          infectable[NSTATES] ;


	real beta, gamma_h, gamma_y, gamma_a, sigma, eta, tau, omega_y, omega_h, mu;

	real nu[50],pi[50],omega_e[50], omega_a[50] ;

	real rexp[0x10000];
	real rnorm_E[0x10000]; // distribution of times in E.
} GLOBALS;

enum infection_state
{
	s_S, s_E, s_I, s_R 
};


typedef struct NODE {
	unsigned int deg; // degree (num of nb and w)
	real w_deg_sum ;
	unsigned int *nb ; // neighbors
	unsigned int *w;  // edge weight index

	unsigned int state, next_state ; // state also includes things like age risk etc.
	unsigned int who_infected ;      // who is responsible for next infection. Important so that it can be taken back when needed.
	real time, when_infected;        // time of next event, when infection happened.

	unsigned int heap; // place on heap
} NODE;


// heap.c
extern void up_heap (unsigned int);
extern void down_heap( unsigned int) ;
extern void del_root ();

// misc.c
extern void init_rng ();
extern void read_data (FILE *, unsigned int w);
extern void read_data3 (FILE *);
extern void read_data5 (FILE *);
void gen_full_nwk( unsigned int n, unsigned int clust, unsigned int w) ;
void gen_nwk( unsigned int n, unsigned int deg, unsigned int w) ;

// pcg_rnd.c
extern uint16_t pcg_16 ();
extern uint8_t pcg_8() ;
extern uint32_t pcg_32 ();
extern uint32_t pcg_32_bounded ();
extern void pcg_init ( uint64_t seed);

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
