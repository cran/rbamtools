/*
 * align_counts.h
 *
 *  Created on: 25.11.2014
 *      Author: kaisers
 */

#ifndef ALIGN_COUNTS_H_
#define ALIGN_COUNTS_H_

#include <Rconfig.h>

///////////////////////////////////////////////////////////////////////////////
// basic definitions


typedef struct seg_align_counts
{
	int * position;	// Genomic positions of segment boundaries
	int * count;	// Number of aligns in each segment
	int n;			// Length of position - 1 = (number of segments) - 2
					// Used as criterion in seg_align_count function
	int index;		// Index of current counting segment (0 <= index < n)

} seg_align_counts;


///////////////////////////////////////////////////////////////////////////////
// basic functions

seg_align_counts * init_align_counts()
{
	seg_align_counts * c = (seg_align_counts*) calloc(1, sizeof(seg_align_counts));
	return c;
}

static R_INLINE void seg_align_counts_destroy(seg_align_counts *c)
{
	// position and count shall be pointers from SEXP
	// so they are not free'd here
	free(c);
}

void seg_align_count(seg_align_counts *c, const bam1_t *align)
{
	while( (align->core.pos >= c->position[(c->index) + 1]) && (c->index < c->n) )
		++c->index;
	++(c->count[c->index]);
}


#endif /* ALIGN_COUNTS_H_ */
