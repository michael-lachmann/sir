// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
// code for SIR on networks by Petter Holme (2018)

// routines for maintaining the binary heap (for the priority queue)
// the root of the heap is 1 (although it is allocated from 0) for simplicity

#include "sir.h"

extern NODE *nd;
extern GLOBALS g;

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
// performing down_heap (a.k.a. percolate down)
// it restores the heap property if there is an inconsistency between 'here'
// and its children (and no other inconsistencies)

// Heap - g.heap is an array, key nd[].time
// for element g.heap[i] on the heap, g.heap[i*2] and g.heap[i+2+1]
// are the roots of the left and right trees.
// The el

void down_heap (unsigned int here) {
	unsigned int utmp, smallest, left, right ; 
		
	smallest = here;
	left = here  *2; // = here * 2 (I know this is silly and saves no time)
	right = left +1; // = left + 1

	// if the heap property is violated vs the children, find the smallest child 
	if ((left <= g.nheap) && (nd[g.heap[left]].time < nd[g.heap[smallest]].time))
		smallest = left;
	if ((right <= g.nheap) && (nd[g.heap[right]].time < nd[g.heap[smallest]].time))
		smallest = right;

	if (smallest != here) {

		// swap smallest and here
		utmp = g.heap[smallest];
		g.heap[smallest] = g.heap[here];
		g.heap[here] = utmp;

		nd[g.heap[smallest]].heap = smallest;
		nd[g.heap[here]].heap = here;

		// continue checking below
		down_heap(smallest);
	}
}

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
// deleting the root of the heap

void del_root () {

	g.heap[1] = g.heap[g.nheap--];
	nd[g.heap[1]].heap = 1;
	down_heap(1);
}

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
// performing up_heap (a.k.a. percolate up)
// for adding an element to the heap
//
// up the heap means earlier and earlier times.
// bubble the start element up till we find an ancestor that is earlier.
// since all ancestors are earlier than their decendents, heap property is kept
// when they move down.

void up_heap (unsigned int start) {
	unsigned int above, here, mem_start ;
	real start_time ;
	mem_start = g.heap[start] ;
		
	here = start ;
	start_time = nd[g.heap[start]].time;

	while (here > 1) {
		above = here >>1 ; // = here / 2 // direct ancestor

		
		if ( start_time >= nd[g.heap[above]].time) break; // found a younger one
		nd[ g.heap[above] ].heap = here; // above node points to here now 
		g.heap[here] = g.heap[above]; // place [here] on heap points to above node 
		// nd[ g.heap[here] ].heap = here; // used to be like this, but harder to understand for me.
		
		here = above; // now that above node moved down, we can continue with that place
	}
	
	g.heap[here] = mem_start;
	nd[ mem_start ].heap = here;
}

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
