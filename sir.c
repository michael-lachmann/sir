// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
// code for SIR on networks by Petter Holme (2018)
//#define TIME
#include "sir.h"

GLOBALS g;
NODE *nd;

#define DEBUG(...)


void set_parameters()
{


/*	Symp_H_Ratio_dict = df_US_dict['YHR']
		Symp_H_Ratio_L_dict = df_US_dict['YHR_low']
		Symp_H_Ratio_H_dict = df_US_dict['YHR_high']
		Hosp_F_Ratio_L_dict = df_US_dict['HFR_low']
		Hosp_F_Ratio_H_dict = df_US_dict['HFR_high']
		Symp_H_Ratio = np.array([Symp_H_Ratio_dict[i] for i in age_group_dict[n_age]])
		Symp_H_Ratio_w_risk = np.array([[Symp_H_Ratio_L_dict[i] for i in age_group_dict[n_age]], \
			[Symp_H_Ratio_H_dict[i] for i in age_group_dict[n_age]]])
		Hosp_F_Ratio_w_risk = np.array([[Hosp_F_Ratio_L_dict[i] for i in age_group_dict[n_age]], \
			[Hosp_F_Ratio_H_dict[i] for i in age_group_dict[n_age]]])

	R0 = 2.2
		DOUBLE_TIME = {'high': 4., 'low': 7.2, 'fit': 2.797335}
		T_EXPOSED_PARA = [5.6, 7, 8.2]
		T_EXPOSED_DIST = lambda x: np.random.triangular(*x)
		T_Y_TO_R_PARA = np.array([21.1, 22.6, 24.4])
		T_Y_TO_R_DIST = lambda x: np.random.triangular(*x)
		T_H_TO_R = 14.
		ASYMP_RATE = 0.179
		PROP_TRANS_IN_E = 0.126
		T_ONSET_TO_H = 5.9
		T_H_TO_D = 14.
		H_RELATIVE_RISK_IN_HIGH = 10
		D_RELATIVE_RISK_IN_HIGH = 10
		HIGH_RISK_RATIO = {'0-4': 8.2825, '5-17': 14.1121, '18-49': 16.5298, '50-64': 32.9912, '65+': 47.0568}  # in %
		H_FATALITY_RATIO = {'0-9': 0., '10-19': 0.2, '20-29': 0.2, '30-39': 0.2, '40-49': 0.4, \
		'50-59': 1.3, '60-69': 3.6, '70-79': 8, '80+': 14.8}  # in %
		INFECTION_FATALITY_RATIO = {'0-9': 0.0016, '10-19': 0.007, '20-29': 0.031, '30-39': 0.084, '40-49': 0.16, \
		'50-59': 0.6, '60-69': 1.9, '70-79': 4.3, '80+': 7.8}  # in %
		OVERALL_H_RATIO = {'0-9': 0.04, '10-19': 0.04, '20-29': 1.1, '30-39': 3.4, '40-49': 4.3, \
		'50-59': 8.2, '60-69': 11.8, '70-79': 16.6, '80+': 18.4}  # in %


		symp_h_ratio_overall
		[0.00048721, 0.00048721, 0.03287572, 0.11337395, 0.17733063]
		hosp_f_ratio
		[[0.04      , 0.12365475, 0.03122403, 0.10744644, 0.23157691],
		 [0.04      , 0.12365475, 0.03122403, 0.10744644, 0.23157691]]
		symp_h_ratio
		[[2.79135866e-04, 2.14621858e-04, 1.32154040e-02, 2.85633688e-02, 3.38733218e-02],
		 [2.79135866e-03, 2.14621858e-03, 1.32154040e-01, 2.85633688e-01, 3.38733218e-01]]

		      ,Hosp       , Death      ,  Death | Hosp
		 0-4  ,0.000487211, 1.94884E-05, 0.04
		 5-17 ,0.000487211, 6.02459E-05, 0.123654749
		 18-49,0.032875723, 0.001026513, 0.031224031
		 50-64,0.113373952, 0.012181628, 0.107446442
		 65+  ,0.177330634, 0.04106568,  0.231576911


*/
	real 
		T_Y_TO_R_PARA = 22.6 , // median of  [21.1, 22.6, 24.4]
		T_EXPOSED_PARA = 7.0, // median of  [5.6, 7, 8.2]
		ASYMP_RATE = 0.179,
		PROP_TRANS_IN_E = 0.126,
		T_ONSET_TO_H = 5.9,
		T_H_TO_R = 14.0 ,
		T_H_TO_D = 14.0 ;
	real symp_h_ratio_overall[] = {0.00048721, 0.00048721, 0.03287572, 0.11337395, 0.17733063} ;

	real hosp_f_ratio[] = {0.04 , 0.12365475, 0.03122403, 0.10744644, 0.23157691} ;

	real symp_h_ratio[] = 
		{2.79135866e-04, 2.14621858e-04, 1.32154040e-02, 2.85633688e-02, 3.38733218e-02, 
		2.79135866e-03, 2.14621858e-03, 1.32154040e-01, 2.85633688e-01, 3.38733218e-01  } ; // for now only first row is used because there is only one risk group.

	real symp_h_ratio_corrrect[] = {
		1.610326595443765, 1.924464960134284, 2.31133016137442, 3.724051596082457, 4.95257504190157
	} ; // hospitalization is 10 times higher for high risk individuals. I made an average of 1*[num low]+10*[num high] for each age class.

	real gamma_y_c, sigma_c ;
	int i ;

	g.gamma_h = 1.0 / T_H_TO_R ;
	gamma_y_c = 1.0 / T_Y_TO_R_PARA ;
	g.gamma_y = gamma_y_c ;
	g.gamma_a = g.gamma_y ;

	sigma_c = 1.0 / T_EXPOSED_PARA ;
	g.sigma = sigma_c ;

	g.eta = 1.0 / T_ONSET_TO_H ;

	g.mu = 1.0 / T_H_TO_D ;

	g.tau = 1.0 - ASYMP_RATE ;

	g.omega_y = 1.0 ;
	g.omega_h = 0.0 ;

	for( i=0; i<g.ngroup; i++)
		symp_h_ratio[i] *= symp_h_ratio_corrrect[i] ;

	for( i=0; i<g.ngroup; i++) {
		g.nu[i] = hosp_f_ratio[i] * g.gamma_h / (g.mu + (g.gamma_h - g.mu  ) * hosp_f_ratio[i] ) ; // hosp_f_ratio is an array of size age.
		g.pi[i] = symp_h_ratio[i] * g.gamma_y / (g.eta + (g.gamma_y - g.eta) * symp_h_ratio[i] ) ; // symp_h_ratio is an array
		// omega_e - relative infectiousness in E, IY, IA // symp_h_ratio_overall length age. 
		g.omega_e[i] = ((symp_h_ratio_overall[i] / g.eta) + ((1 - symp_h_ratio_overall[i] ) / g.gamma_y)) * 
						g.omega_y * g.sigma * PROP_TRANS_IN_E / (1 - PROP_TRANS_IN_E) ;
		g.omega_a[i] = ((symp_h_ratio_overall[i] / g.eta) + ((1 - symp_h_ratio_overall[i]) / gamma_y_c)) * 
						g.omega_y * sigma_c * PROP_TRANS_IN_E / (1 - PROP_TRANS_IN_E) ;
	}
	

}


char *SN="SEAYHRD" ;

print_state(state s) {
	int a = s % 5 ;
	int i = s / 5 ;
	printf("%c%d",SN[i],a) ;
}



// One question is if to go through the whole disease progression when we pick an individual up.
// Is the heap just for infections or for any events?
// Events other than infections are independent. No one cares how it progresses.
// Therefore, putting them on a heap is a waste. You don't care if they occured before or after someone else, they do not cause any further events in other individuals.

// However, in the Meyer lab mode, you are still infectious when you do transitions. So mayeb then it does make sense to put them on the heap?
// So maybe the rule should be do as little ass possible so you don't have to do things twice?

// If we knew when the next change in the network would be, we could compare to it. But is it worth it to compare all the time, when it is so rare?



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


real myQexp(real q) {
	return -log(q) ;
}

state next_state(state s) {
	int i=0;
	while (g.state_trans_p[s][i] < MAXRAND) {
		if( pcg_32() < g.state_trans_p[s][i] ) return g.state_trans[s][i] ;
		i++ ;
	}
	return g.state_trans[s][i] ;
}

void conditional_infect_neighbours( real t_me, unsigned int me, real rate) {
	unsigned int i, you ; 
	real t ;

	// go through my neighbors 
	for ( i = 0; i < nd[me].deg; i++) {
		you = nd[me].nb[i];
		if (g.infectable[ nd[you].state ]) { // if you can be infected
			t = g.now + 
				g.time_dist[ nd[you].state ][RANDOM_QUART()] / (g.weight[ nd[me].w[i]  ] *rate ) ; // get the infection time
			if ((t < t_me ) &&     // When me becomes resistant
				(t < nd[you].time)) {    // when you become exposed, if you do
				DEBUG( printf("%g %g inf:%d\n",t_me,t,you) ) ;
				nd[you].time = t;
				nd[you].next_state = next_state( nd[you].state ) ; // ***** NOT RANDOM YET, because we always switch to E.
				nd[you].who_infected = me ;
				nd[you].when_infected = g.now ;
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

void conditional_infect_neighbours_down( real t_me, unsigned int me, real rate) {
	unsigned int i, you ; 
	real t, old_t ;


	// ****************************************** This is the problem: if I infected you in the past, but now I chose a time that is larger than when I stay infectious, then 
	// ****************************************** I don't take things back. Also need to figure out what to do if multiple infections. 
	// What really needs to be done is to take back all infections and reroll.
	// go through my neighbors 
	for ( i = 0; i < nd[me].deg; i++) {
		you = nd[me].nb[i];
		if (g.infectable[ nd[you].state ]) { // if you can be infected
			t = g.now + 
				g.time_dist[ nd[you].state ][RANDOM_QUART()] / (g.weight[ nd[me].w[i]  ] *rate ) ; // get the infection time
			if (t < t_me )  {  // When me becomes resistant
				if( (nd[you].when_infected < g.now-EPSILON) || // infection is old, can ignore
					(t < nd[you].time) // infection is before your next event
				) {    // when you become exposed, if you do
					DEBUG( printf("%g %g inf:%d\n",t_me,t,you) ) ;
					old_t = nd[you].time ;
					nd[you].time = t;
					nd[you].next_state = next_state( nd[you].state ) ; // ***** NOT RANDOM YET, because we always switch to E.
					nd[you].who_infected = me ;
					nd[you].when_infected = g.now ;
					if (nd[you].heap == NONE) { // if not listed before, then extend the heap
						g.heap[++g.nheap] = you;
						nd[you].heap = g.nheap;
					}
					if( t < old_t)
						up_heap(nd[you].heap); // this works bcoz the only heap relationship 
											   // that can be violated is the one between you and its parent
					else
						down_heap(nd[you].heap);

					continue ;
				}
			}
			if (nd[you].who_infected == me) { // I didn't infect you this time, but it was my fault before
				nd[you].time = NEVER ;
				down_heap( nd[you].heap) ;
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
		if (g.infectable[ nd[you].state ]) { // if you is S, you can be infected
			t = g.now + g.rexp[pcg_16()] / (g.weight[ nd[me].w[i] ]*g.beta)  * nd[me].w_deg_sum ; // get the infection time. (weight already includes 1/beta)
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



// basically we get to a node, and change its state, make its time earlier
// or later.
// in special cases, we infect other nodes, and make their time earlier or later
// (later only in the special case of changing parameters)




void epi_timestep( ) {
	unsigned int me ;
	state next_s, s ;
	real  next_t ;



		me = g.heap[1] ;
		g.now = nd[me].time ;

	
		// Advance node to its next state
//		nd[me].time  = nd[me].next_time ;

		g.s[ nd[me].state      ]-- ;
		g.s[ nd[me].next_state ]++ ;
		DEBUG( printf("%d ",me) ) ;
		DEBUG( print_state( nd[me].state) ) ;
		DEBUG( printf("->") ) ;
		DEBUG( print_state( nd[me].next_state) ) ;
		DEBUG( printf("\n") ) ;
		nd[me].state = nd[me].next_state ;


		// state detemines what happens next
		s = nd[me].state ;


		// Node just switched to state s. 
		// What to do next with it?
		// If it is in a state that will change in the future on its own, then schedule that.
		// If it is a state that can infect others, schedule those infections.
		if( g.self_change[s] ) { // State is one that will change on its own in the future
			// time of next own change
			next_t = g.now + g.time_dist [s][RANDOM_QUART()] ; 
			// choose next state
			next_s =    next_state(s) ;

			// If infection rate is bigger than 0, test for infections.
			if( g.state_infect_rate[s] > EPSILON ) 
				conditional_infect_neighbours( next_t, me, g.state_infect_rate[s]) ; // only infect neighbours that get infected before switch to next state.

			// schedule next event for node
			nd[me].next_state = next_s ;
			nd[me].time       = next_t ;
			down_heap( nd[me].heap ) ;
		} else { // State doesn't change on its own, schedule a change for NEVER.
			// We switched to a state that should never change
			nd[me].time = NEVER ;
			// next_state stays what it is
			down_heap( nd[me].heap ) ;
		}
}



void reinfect_s_I_neighbours()
{
	unsigned int me, i  ;
	state s;

	for (i = 1; i <= g.nheap; i++) {
		me = g.heap[i] ;
		s = nd[me].state ;
		if ( g.state_infect_rate[ s ] > EPSILON) {
			conditional_infect_neighbours_down(  nd[me].time, me,  g.state_infect_rate[ s ]) ;
		}
	}
}


// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
// this routine runs one SIR outbreak from a random seed node

void seir_init( unsigned int nS) {
	unsigned int i, source;
	const state S=0 ;
	// g.t = 0.0;

	for(i=0; i<NSTATES; i++)
		g.s[i] = 0 ;


	// initialize
	for (i = 0; i < g.n; i++) {
		nd[i].heap = NONE;
		nd[i].time = NEVER; // to a large value
		nd[i].state = nd[i].state % g.ngroup ; // reset to modulo ngroup
		g.s[ nd[i].state ]++ ;
	}

	for( i=1; i <= nS; i++) {
	// get & infect the source
		do {
			source = pcg_32_bounded(g.n);
		} while( nd[source].time==0.0 ) ;
		nd[source].time = 0.0;
		nd[source].next_state = nd[source].state+ 1* g.ngroup ; // assuming for each group, 0 is S, 1 is E.
		nd[source].heap = i;
		g.heap[ g.nheap = i ] = source;
	}

}



struct _named_seir_args {
	state S, del;
	real beta, gammaY, gammaA, gammaH, 
		omegaY, omegaA, omegaE, omegaH,
		sigma, tau, eta, mu, nu, pi ;
};

void _init_mats_seir_full( struct _named_seir_args *args) {
	int i ;
	state E = args->S+  args->del,
		  IA= args->S+2*args->del, 
		  IY= args->S+3*args->del, 
		  IH= args->S+4*args->del, 
		  R = args->S+5*args->del, 
		  D = args->S+6*args->del,
		  S = args->S;

	real rateIY, rateIH, frac1 ;


	for (i = S; i <= D;  i+= args->del) {
		g.time_dist[i] = malloc( NQUART * sizeof(real) ) ;
	}


	// S
	for( i=0; i< NQUART; i++) g.time_dist [S][i] = myQexp( (1.0+i)/NQUART ) / args->beta ;
	g.state_trans_p[S][0] = MAXRAND ;  g.state_trans[S][0]   = E ;	

	g.state_infect_rate[S] = 0.0 ;
	g.self_change[S] = FALSE ;
	g.infectable[S] = TRUE ;

	// E
	for( i=0; i< NQUART; i++) g.time_dist [E][i] = myQexp( (1.0+i)/NQUART ) / args->sigma;
	g.state_trans_p[E][0] = args->tau * MAXRAND ; g.state_trans[E][0]   = IY ;	
	g.state_trans_p[E][1] = MAXRAND ;             g.state_trans[E][1]   = IA ;	

	g.state_infect_rate[E] = args->omegaE ;
	g.self_change[E] = TRUE ;
	g.infectable[E] = FALSE ;

	// IY
	rateIY = (1- args->pi) * args->gammaY + args->pi * args->eta;
	frac1 =  (1- args->pi) * args->gammaY                        / rateIY ;

	for( i=0; i< NQUART ; i++) g.time_dist [IY][i] = myQexp( (1.0+i)/NQUART ) / rateIY;
	g.state_trans_p[IY][0] = frac1 * MAXRAND;  g.state_trans[IY][0]   = R ;
	g.state_trans_p[IY][1] = MAXRAND;          g.state_trans[IY][1]   = IH ;

	g.state_infect_rate[IY] = args->omegaY ;
	g.self_change[IY] = TRUE ;
	g.infectable[IY] = FALSE ;

	// IA
	for( i=0; i<NQUART; i++) g.time_dist [IA][i] = myQexp( (1.0+i)/NQUART ) / args->gammaA;
	g.state_trans_p[IA][0] = MAXRAND;     g.state_trans[IA][0]   = R ;

	g.state_infect_rate[IA] = args->omegaA ;
	g.self_change[IA] = TRUE ;
	g.infectable[IA] = FALSE ;

	// IH
	rateIH = (1- args->nu) * args->gammaH  +  args->nu * args->mu ;
	frac1  = (1- args->nu) * args->gammaH                          / rateIH ;

	for( i=0; i< NQUART ; i++) g.time_dist [IH][i] = myQexp( (1.0+i)/NQUART ) / rateIH;
	g.state_trans_p[IH][0] = frac1 * MAXRAND; g.state_trans[IH][0]   = R ;
	g.state_trans_p[IH][1] = MAXRAND;         g.state_trans[IH][1]   = D ;
	
	g.state_infect_rate[IH] = args->omegaH ;
	g.self_change[IH] = TRUE ;
	g.infectable[IH] = FALSE ;


	// R
	for( i=0; i<NQUART; i++) g.time_dist [R][i] = NEVER ;
	g.state_trans_p[R][0] = MAXRAND; g.state_trans[R][0]   = R ;

	g.state_infect_rate[R] = 0.0 ;
	g.self_change[R] = FALSE ;
	g.infectable[R] = FALSE ;

	// D
	for( i=0; i<NQUART; i++) g.time_dist [R][i] = NEVER ;
	g.state_trans_p[D][0] = MAXRAND ;  g.state_trans[D][0]   = D ;

	g.state_infect_rate[R] = 0.0 ;
	g.self_change[R] = FALSE ;
	g.infectable[R] = FALSE ;


}

#define init_mats_seir_full(...) _init_mats_seir_full( &(struct _named_seir_args){__VA_ARGS__})

void free_mats_seir_full( state S, state del) {
	int i ;
	state D= S+6*del ;
	for (i = S; i < D; i+=del) {
		free( g.time_dist[i]) ;
	}
}


// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
// main function handling input


int main (int argc, char *argv[]) {
	unsigned int i,d,uk;
	real par, R0 , last_now;
	int pos,npos,k, T_i,r ;
	unsigned int w ;
	char arg_type, par_type ;

	FILE *fp;
#ifdef TIME
	struct timespec t0, t1;
#endif
	pcg_init( 42) ;
	g.w_n = 0 ;
	g.av_deg = 0.0 ;
	g.ngroup = 1 ;

	for( k=1; k<argc; k++) {
	sscanf(argv[k],"%c",&arg_type) ;
	switch( arg_type) {
	// Generate networks
	//
	// 
	case 'N': // global contact network network
		sscanf(argv[k],"%c,%d,%d,%u",&arg_type,&i,&d,&w) ;
		gen_nwk( i,d, w ) ;
		break ;
	case 'H': // home networks 
		sscanf(argv[k],"%c,%d,%d,%u",&arg_type,&i,&d,&w) ;
		if( i % d != 0) {
			fprintf(stderr,"for full network, size (%d) needs to be a multiple of clust (%d)\n",i,d) ;
			return 1 ;
		}
		gen_full_nwk( i, d, w) ;
		break ;
	// specify SEIR model
	case 'S': // SEIR model S,R[R0] or S,B[beta]
		sscanf(argv[k],"%c,%c%lg",&arg_type,&par_type,&par) ;
		if (g.av_deg == 0) {
			perror("Please read network file first in parameters\n") ;
			exit(1) ;
		}
		if (par_type == 'B') {
			g.beta = par ;
		}
		else if (par_type == 'R') {
			R0 = par ;
			g.beta = R0 * g.gamma_y ;  /// we / ( g.av_deg - R0) ;
			fprintf( stderr, "calculated beta is %lg\n",g.beta) ;
		} else {
			perror("Please enter type of 2nd parameter: B or R\n") ;
			exit(1) ;
		}
		break ;
	case 'R': // random seed
		sscanf( argv[k],"%c,%u",&arg_type,&i) ;
		pcg_init( (uint64_t)i );
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
	case 'f': // read network with groups from file
		fp = fopen(argv[k]+2, "r"); // skip over 'F,'
		if (!fp) {
			fprintf(stderr, "can't open '%s'\n", argv[1]);
			return 1;
		}
		read_data5(fp);
		fclose(fp);
		break ;
	case 'w': // read periods of weight.
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

	fprintf(stderr,"done\n") ;
	g.weight = malloc( g.nweight * sizeof(real)) ;
	// allocating the heap (N + 1) because its indices are 1,...,N
	g.heap = malloc((g.n + 1) * sizeof(unsigned int));

	//init_mats_seir( g.beta, g.gamma, g.sigma) ;

	set_parameters() ;
	for( i=0; i<g.ngroup; i++) {
		init_mats_seir_full( .S=i, .del=g.ngroup,   
			.beta = g.beta, .gammaA = g.gamma_a, .gammaH= g.gamma_h, .gammaY = g.gamma_y, .sigma = g.sigma,
			.omegaY = g.omega_y, .omegaA = g.omega_a[i], .omegaE = g.omega_e[i], .omegaH = g.omega_h ,
			.sigma = g.sigma, .mu= g.mu, .nu= g.nu[i], .pi = g.pi[i], .tau = g.tau, .eta= g.eta       ) ;
	}


	for( T_i=0; T_i < g.w_n; T_i++) {
		for( i=1; i< g.nweight; i++) {
			g.weight[i] = (g.w_val[  T_i*(g.nweight) +i ] ) ;
		}
		if( T_i==0)
			seir_init( 10) ;
		else
			reinfect_s_I_neighbours() ;

		last_now = g.now ;
		while( g.now < g.w_time[T_i] ) {
			epi_timestep(  ) ;
			if (g.now - last_now > 1) {
				printf("%g ",g.now) ;
				for( i=0; i< g.ngroup*7; i++)
					printf("%d ",g.s[i]) ;
				printf("\n") ;
				last_now = g.now ;
			}
		}
	}


	// print result
	// printf("avg. outbreak size: %lg (%lg)\n", ss1, sqrt((ss2 - SQ(ss1)) / (NAVG - 1)));
	// printf("avg. time to extinction: %lg (%lg)\n", st1, sqrt((st2 - SQ(st1)) / (NAVG - 1)));

	// cleaning up


	for (i = 0; i < g.n; i++) {
		free(nd[i].nb);
		free(nd[i].w);
	}
	free( nd); 
	free( g.heap);
	free( g.w_time) ;
	free( g.w_val) ;
	for( i=0; i<g.ngroup; i++)
		free_mats_seir_full( i,  g.ngroup) ;


	return 0;
}

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
