// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
// code for SIR on networks by Petter Holme (2018)


#include "sir.h"

extern NODE *nd;
extern GLOBALS g;

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
// reads the network, assumes an edge list with vertex label 0,N-1
// if your network has nodes with degree zero, make sure that none of them is
// the node with largest index

void read_data (FILE *fp, unsigned int w) {
	unsigned int i, me, you;

	g.n = 0;

	// scan the system size
	while (2 == fscanf(fp, "%u %u\n", &me, &you)) {
		if (g.n < me) g.n = me;
		if (g.n < you) g.n = you;
	}

	g.n++;

	nd = calloc(g.n, sizeof(NODE));

	rewind(fp);

	// scan the degrees
	while (2 == fscanf(fp, "%u %u\n", &me, &you)) {
		nd[me].deg++;
		nd[you].deg++;
	}

	// allocate adjacency lists
	for (i = 0; i < g.n; i++) {
		nd[i].nb = malloc(nd[i].deg * sizeof(unsigned int));
		nd[i].w =  malloc(nd[i].deg * sizeof(real));
		nd[i].deg = 0;
		nd[i].state = s_S ;
	}

	rewind(fp);

	// fill adjacency lists
	while (2 == fscanf(fp, "%u %u\n", &me, &you)) {
		nd[me].nb[ nd[me].deg] = you;
		nd[me].w [ nd[me].deg  ] = w ;
		nd[me].deg++ ;

		nd[you].nb[ nd[you].deg] = me;
		nd[you].w [ nd[you].deg  ] = w ;
		nd[you].deg++ ;
	}
}

void read_data3 (FILE *fp) {
	unsigned int i, me, you, w, max_w;

	g.n = 0;

	// scan the system size
	max_w = g.nweight ;
	while (3 == fscanf(fp, "%u %u %u\n", &me, &you, &w)) {
		if (g.n < me) g.n = me;
		if (g.n < you) g.n = you;
		if( w > max_w) max_w = w ;
	}
	g.nweight = max_w+1 ;

	g.n++;

	nd = calloc(g.n, sizeof(NODE));

	rewind(fp);

	// scan the degrees
	while (3 == fscanf(fp, "%u %u %u\n", &me, &you, &w)) {
		nd[me].deg++;
		nd[you].deg++;
	}

	// allocate adjacency lists
	for (i = 0; i < g.n; i++) {
		nd[i].nb = malloc(nd[i].deg * sizeof(unsigned int));
		nd[i].w =  malloc(nd[i].deg * sizeof(unsigned int));
		nd[i].deg = 0;
		nd[i].state = s_S ;
	}

	rewind(fp);


	// fill adjacency lists
	while (3 == fscanf(fp, "%u %u %u\n", &me, &you, &w)) {
		nd[me].nb[ nd[me].deg] = you;
		nd[me].w [ nd[me].deg  ] = w ;
		nd[me].deg++ ;

		nd[you].nb[ nd[you].deg] = me;
		nd[you].w [ nd[you].deg  ] = w ;
		nd[you].deg++ ;
	}

}


void gen_nwk( unsigned int n, unsigned int deg, unsigned int w) {
	unsigned int i,j ;
	for( i=0; i<n; i++) {
		for( j=0; j<deg; j++) {
			printf("%u %u %u\n",i,pcg_32_bounded(n),w) ;
		}
	}
}

void gen_full_nwk( unsigned int n, unsigned int clust, unsigned int w) {
	unsigned int *m ;
	unsigned int i,j,k, tmp, rep_i ;
	m = malloc( n * sizeof( unsigned int) ) ;
	if( m == NULL) {
		perror("malloc failed in gen_full_nwk\n") ;
		exit(1) ;
	}

	g.nweight = 0 ;

	for( i=0; i<n; i++)  
		m[i] = i ;

	for( i=0; i<n; i++) {
		tmp = m[i] ;
		rep_i = pcg_32_bounded( n) ;  
		m[i] = m[ rep_i] ;
		m[rep_i] = tmp ; 
	}


	for( i=0; i<n; i+= clust) {
		for( j=0; j< clust-1; j++ ) {
			for( k=j+1; k<clust; k++ ) {
				printf("%u %u %u\n", m[i+j], m[i+k],w) ;
			}
		}
	}
	free(m) ;
}


