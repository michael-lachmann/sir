// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
// code for SIR on networks by Petter Holme (2018)
//#define TIME
#include "sir.h"
//test


GLOBALS g;
NODE *nd;

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
// this routine does the bookkeeping for an infection event


//Quick & dirty, tolerance under +-6e-3. 
// Work based on "A handy approximation for the error function and its inverse" by Sergei Winitzki.
// https://stackoverflow.com/questions/631629/erfx-and-math-h
real myErfInv2(real x){
   real tt1, tt2, lnx, sgn;
   sgn = (x < 0) ? -1.0 : 1.0;

   x = -x*x + 1.0;        // x = 1 - x*x;
   lnx = log(x);

   tt1 = 2.0/(M_PI*0.147) + 0.5 * lnx;
   tt2 = lnx/(0.147) ;

   return(sgn*sqrt(-tt1 + sqrt(tt1*tt1 - tt2)));
}

real myQnorm( real p) {
	return( M_SQRT2 * myErfInv2(2.0 * p -1.0) );
}

void infect () {
	unsigned int i, you, me = g.heap[1];
	real t, now = nd[me].time;

	del_root();
	nd[me].heap = I_OR_R;
	// get the recovery time
	nd[me].time += g.rexp[pcg_16()] * g.beta; // bcoz g.rexpr has a / g.beta factor
	// if (nd[me].time > g.t) g.t = nd[me].time;

	// go through the neighbors of the infected node . .
	for (i = 0; i < nd[me].deg; i++) {
		you = nd[me].nb[i];
		if (nd[you].heap != I_OR_R) { // if you is S, you can be infected
			t = now + g.rexp[pcg_16()] * nd[me].w[i] ; // get the infection time

			if ((t < nd[me].time) && (t < nd[you].time)) {
				nd[you].time = t;
				if (nd[you].heap == NONE) { // if not listed before, then extend the heap
					g.heap[++g.nheap] = you;
					nd[you].heap = g.nheap;
				}
				up_heap(nd[you].heap); // this works bcoz the only heap relationship 
									  // that can be violated is the one between you and its parent
			}
		}
	}
}



void conditional_infect_neighbours( unsigned int me) {
	unsigned int i, you ; 
	real t;

	// go through the neighbors of the infected node . .
	for ( i = 0; i < nd[me].deg; i++) {
		you = nd[me].nb[i];
		if (nd[you].state == s_S) { // if you is S, you can be infected
			t = g.now + g.rexp[pcg_16()] * g.weight[ nd[me].w[i] ] ; // get the infection time
			if ((t < nd[me].time) &&     // When me becomes resistant
				(t < nd[you].time)) {    // when you become exposed, if you do
				nd[you].time = t;
				if (nd[you].heap == NONE) { // if not listed before, then extend the heap
					g.heap[++g.nheap] = you;
					nd[you].heap = g.nheap;
				}
				up_heap(nd[you].heap); // this works bcoz the only heap relationship 
									   // that can be violated is the one between you and its parent
			}
		}
	}
}


void unconditional_infect_neighbours( unsigned int me) {
	unsigned int i, you ; 
	real t, old_t;

	// go through the neighbors of the infected node . .
	for ( i = 0; i < nd[me].deg; i++) {
		you = nd[me].nb[i];
		if (nd[you].state == s_S) { // if you is S, you can be infected
			t = g.now + g.rexp[pcg_16()] * g.weight[ nd[me].w[i] ] ; // get the infection time. (weight already includes 1/beta)
			if ( t < nd[me].time) {    // when you become exposed, if you do
				old_t = nd[you].time ;
				nd[you].time = t;
				if (nd[you].heap == NONE) { // if not listed before, then extend the heap
					g.heap[++g.nheap] = you;
					nd[you].heap = g.nheap;
				}
				if( t < old_t)
					up_heap(nd[you].heap); // this works bcoz the only heap relationship 
										   // that can be violated is the one between you and its parent
				else
					down_heap(nd[you].heap);
			}
		}
	}
}


// basically we get to a node, and change it's state, make it's time earlier
// or later.
// in special cases, we infect other nodes, and make their time earlier or later
// (later only in the special case of changing parameters)

void infect_SEIR () {
	unsigned int me ; 

	me = g.heap[1] ;
	g.now = nd[me].time ;


	switch( nd[me].state ) { 
		case s_S: // Became exposed but not infectious
			nd[me].time += g.rnorm_E[pcg_16()] ; // Normal with variance 1, mean 4
			down_heap( nd[me].heap ) ;

			nd[me].state = s_E ;
			g.s[s_S]-- ;
			g.s[s_E]++ ;
			break ;
		case s_E: // Finished exposed period, start spreading, switch to I.
			nd[me].time += g.rexp[pcg_16()] * g.INV_gamma; // Time to become resistant
			down_heap( nd[me].heap ) ;

			nd[me].state = s_I ;
			g.s[s_E]--;
			g.s[s_I]++;

			conditional_infect_neighbours( me) ;
			break ;	
		case s_I: // finished infectious period, resistant
			del_root(); // This is only because R doesn't do anything any more. Totally taken off heap.
			nd[me].state = s_R ;
			g.s[s_I]--;
			g.s[s_R]++;
		break ;
		case s_R:
			perror("We shouldn't be here, tried to do something to resistant node\n") ;
			break ;
	} 

}

void reinfect_s_I_neighbours()
{
	unsigned int i  ;
	for (i = 0; i < g.n; i++) {
		if (nd[i].state == s_I) {
			unconditional_infect_neighbours( i) ;

		}
	}
}


// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
// this routine runs one SIR outbreak from a random seed node

void seir_init( unsigned int nS) {
	unsigned int i, source;

	// g.t = 0.0;
	g.s[(unsigned int)s_S] = g.n ;
	g.s[(unsigned int)s_E] = 0 ;
	g.s[(unsigned int)s_I] = 0 ;
	g.s[(unsigned int)s_R] = 0 ;

	// initialize
	for (i = 0; i < g.n; i++) {
		nd[i].heap = NONE;
		nd[i].time = DBL_MAX; // to a large value
		nd[i].state = s_S ;
	}

	for( i=1; i <= nS; i++) {
	// get & infect the source
		do {
			source = pcg_32_bounded(g.n);
		} while( nd[source].time==0.0 ) ;
		nd[source].time = 0.0;
		nd[source].state = s_S ;
		nd[source].heap = i;
		g.heap[ g.nheap = i ] = source;
	}

}

void seir ( real T) {

	// run the outbreak until either not enough s_S or the heap is empty.
	while ( ( g.now < T) &&  g.nheap ) {
		infect_SEIR();
		printf("%f %d %d %d %d\n",g.now,g.s[0],g.s[1],g.s[2],g.s[3]) ;
//		printf("%d\n",g.nheap) ;
	}
}



// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
// main function handling input


int main (int argc, char *argv[]) {
	unsigned int i,d,uk;
	int pos,npos,k, T_i,r ;
	unsigned int w ;
	char arg_type ;

	FILE *fp;
#ifdef TIME
	struct timespec t0, t1;
#endif
	pcg_init();
	g.w_n = 0 ;

	for( k=1; k<argc; k++) {
	sscanf(argv[k],"%c",&arg_type) ;
	switch( arg_type) {
	case 'N': // total network
		sscanf(argv[k],"%c,%d,%d,%u",&arg_type,&i,&d,&w) ;
		gen_nwk( i,d, w ) ;
		break ;
	case 'H': // home networks 
		sscanf(argv[k],"%c,%d,%d,%u",&arg_type,&i,&d,&w) ;
		if( i % d != 0) {
			printf("for full network, size (%d) needs to be a multiple of clust (%d)\n",i,d) ;
			return 1 ;
		}
		gen_full_nwk( i, d, w) ;
		break ;
	case 'S': // SEIR model
		sscanf(argv[k],"%c,%lg,%lg,%lg,%lg",&arg_type,&g.beta, &g.gamma,&g.Em, &g.Ev) ;
		g.INV_beta = 1.0 / g.beta ;
		g.INV_gamma = 1.0 / g.gamma ;
		for (i = 0; i < 0x10000; i++)
			g.rnorm_E[i] = g.Ev * myQnorm( (1.0 + i) / 0x10000 ) + g.Em ;
		break ;
	case 'R': // random seed
		sscanf( argv[k],"%c,%u",&arg_type,&i) ;
		g.state = (uint64_t)i ;
		break ;
	case 'F': // read network from file
		fp = fopen(argv[k]+2, "r"); // skip over 'F,'
		if (!fp) {
			fprintf(stderr, "can't open '%s'\n", argv[1]);
			return 1;
		}
		read_data3(fp);
		fclose(fp);
		break ;
	case 'w':
		if( g.nweight == 0 ) {
			perror("number of weights is 0, read file first\n") ;
			exit(1) ;
		}
		if( g.w_n==0 ) {
			g.w_time = malloc( (1+100) * sizeof(real)) ;
			g.w_val = malloc( (1+ 100) * g.nweight * sizeof(real)) ;
		} 
		if( g.w_n >= 100) {
			perror("Can't have more than 100 periods\n") ;
			exit(1) ;
		}
		pos = 2 ;
		r = sscanf( argv[k]+pos,"%lg%n",&g.w_time[g.w_n],&npos) ;
		if (r == 1) {
			pos += npos + 1 ;
		}
		for( uk=1; uk< g.nweight; uk++) {
			if (argv[k][pos] == '\0') {
				r = 0 ; npos = 0 ;
			} else {
				r = sscanf( argv[k]+pos,"%lf%n", &g.w_val[ g.w_n*g.nweight + uk],&npos )  ;
				if( r==1) pos += npos ;
				if( argv[k][pos] != '\0') pos++ ;
			}
			if( r != 1) {
				perror("Didn't find weight, assuming 1\n") ;
				g.w_val[ g.w_n*g.nweight + uk] = 1.0 ;
			}
		}
		g.w_n++ ;
		break ;
	}
	}


	g.weight = malloc( g.nweight * sizeof(real)) ;
	// allocating the heap (N + 1) because it's indices are 1,...,N
	g.heap = malloc((g.n + 1) * sizeof(unsigned int));

	for (i = 0; i < 0x10000; i++)
		g.rexp[i] = -log((i + 1.0) / 0x10000) ;

	for( T_i=0; T_i < g.w_n; T_i++) {
		for( i=1; i< g.nweight; i++) 
			g.weight[i] = 1.0 / (g.w_val[  T_i*(g.nweight) +i ] * g.beta) ;
		if( T_i==0)
			seir_init( 10) ;
		else
			reinfect_s_I_neighbours() ;
		seir( g.w_time[T_i] ) ;
	}


	// print result
	// printf("avg. outbreak size: %lg (%lg)\n", ss1, sqrt((ss2 - SQ(ss1)) / (NAVG - 1)));
	// printf("avg. time to extinction: %lg (%lg)\n", st1, sqrt((st2 - SQ(st1)) / (NAVG - 1)));

	// cleaning up
	for (i = 0; i < g.n; i++) {
		free(nd[i].nb);
		free(nd[i].w);
	}
	free(nd); 
	free(g.heap);
	free( g.w_time) ;
	free( g.w_val) ;


	return 0;
}

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
