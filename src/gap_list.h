/*
 *  File      :  gap_list.h
 *
 *  Created on:  12.09.2012
 *  Author    :  Wolfgang Kaisers
 *
 *  Change log:
 *
 *  29.Okt.12 : Corrected Error in list_gaps: Multiple gapped Aligns are now correct handled
 */

#ifndef GAP_LIST_H_
#define GAP_LIST_H_

#include "samtools/sam.h"
#include "samtools/bam.h"

///////////////////////////////////////////////////////////////////////////////////////////////////
// basic definitions


// Number of unsigned int's in gap_data
const unsigned data_size=9;

typedef struct gap_data
{
	unsigned refid;
	unsigned position;
	unsigned left_cigar_len;
	unsigned left_cigar_type;
	unsigned left_stop;
	unsigned gap_len;
	unsigned right_start;
	unsigned right_cigar_len;
	unsigned right_cigar_type;
} gap_data;


typedef struct gap_element
{
	unsigned *pos_data;
	struct gap_element *last_el;
	struct gap_element *next_el;
} gap_element;

typedef struct {
	gap_element *first_el;
	gap_element *last_el;
	gap_element *curr_el;
	unsigned long size;
	unsigned long nAlings;
	unsigned long nGapAlings;
} gap_list;

inline gap_list * init_gap_list()
{
	return (gap_list*) calloc(1,sizeof(gap_list));
}

inline gap_element* init_gap_elem(const gap_data data)
{
	gap_element *e=calloc(1,sizeof(gap_element));
	e->pos_data=calloc(data_size,sizeof(unsigned));
	(e->pos_data)[0]=data.refid;
	(e->pos_data)[1]=data.position;
	(e->pos_data)[2]=data.left_cigar_len;
	(e->pos_data)[3]=data.left_cigar_type;
	(e->pos_data)[4]=data.left_stop;
	(e->pos_data)[5]=data.gap_len;
	(e->pos_data)[6]=data.right_start;
	(e->pos_data)[7]=data.right_cigar_len;
	(e->pos_data)[8]=data.right_cigar_type;
	return e;
}

inline void destroy_gap_elem(gap_element *e)
{
	free(e->pos_data);
	free(e);
}

///////////////////////////////////////////////////////////////////////////////////////////////////
// list generic accessor functions

inline void gap_list_push_back_elem(gap_list *l,gap_element *e)
{
	if(l->size==0)
	{
		l->first_el=e;
		l->last_el=e;
		l->size=1;
	}
	else
	{
		e->last_el=l->last_el;
		e->last_el->next_el=e;
		l->last_el=e;
		++(l->size);
	}
}

inline void gap_list_push_back(gap_list *l,const gap_data data) { gap_list_push_back_elem(l,init_gap_elem(data)); }

inline void gap_list_push_front_elem(gap_list *l,gap_element *e)
{
	if(l->size==0)
	{
		l->first_el=e;
		l->last_el=e;
		l->size=1;
	}
	else
	{
		e->next_el=l->first_el;
		e->next_el->last_el=e;
		l->first_el=e;
		++(l->size);
	}
}

inline void gap_list_push_front(gap_list *l,const gap_data data) { gap_list_push_front_elem(l,init_gap_elem(data)); }

inline void pop_back_gap_list(gap_list *l)
{
	if(l->first_el!=l->last_el)
	{
		gap_element *e=l->last_el;
		e->last_el->next_el=0;
		l->last_el=e->last_el;
		destroy_gap_elem(e);
		--(l->size);
	}
	else if(l->last_el>0)
	{
		destroy_gap_elem(l->last_el);
		l->first_el=0;
		l->last_el=0;
		l->size=0;
	}
}

inline void pop_front_gap_list(gap_list *l)
{
	//printf("[pop_front_gap_list] list size %lu\n",l->size);
	if(l->first_el!=l->last_el)
	{
		gap_element *e=l->first_el;
		e->next_el->last_el=0;
		l->first_el=e->next_el;
		destroy_gap_elem(e);
		--(l->size);
	}
	else if(l->first_el>0)
	{
		destroy_gap_elem(l->first_el);
		l->first_el=0;
		l->last_el=0;
		l->size=0;
	}
}

inline void wind_back_gap_list(gap_list *l) { l->curr_el=NULL; }

///////////////////////////////////////////////////////////////////////////////
// higher level functions

inline void destroy_gap_list(gap_list *l)
{
	while(l->size>0)
		pop_front_gap_list(l);
	free(l);
}

inline void list_gaps(gap_list *l,const bam1_t* align)
{
	uint32_t position,n_cigar,i;
	uint32_t *cigar;
	gap_data g;
	int op;

	cigar=bam1_cigar(align);
	position=align->core.pos;
	n_cigar=align->core.n_cigar;

	g.refid=align->core.tid;
	g.position=position;
	// Count total Aligns
	++(l->nAlings);
	//printf("[list_gaps] refid: %i\tposition: %i\tn_cigar: %i\n",g.refid,g.position,n_cigar);


	if(n_cigar>2)
	{
		// Assumes that first and last cigar type cannot be N!
		// In this loop there's always a valid left and right cigar

		// Add size of first cigar to position to get *1-based* left stop
		// BAM is *0-based*
		position+= (cigar[0] >> BAM_CIGAR_SHIFT);

		for(i=1;i<(n_cigar-1);++i)
		{
			// BAM_CMATCH		= 0 = M match or mismatch
			// BAM_CREF_SKIP	= 3 = N skip on the reference (e.g. spliced alignment) *
			// BAM_CINS		= 1 = S clip on the read with clipped sequence

			// ERROR: First CIGAR item cannot be I or D or N.
			op = cigar[0] & BAM_CIGAR_MASK;
			if((op!= BAM_CMATCH) && (op!=BAM_CINS))
				return;

			op = cigar[i] & BAM_CIGAR_MASK;
			//printf("op %i\t",op);
			if(op==BAM_CREF_SKIP)
			{
				// Count gaps
				++(l->nGapAlings);
				// N -> Add gap to list
				g.left_stop=position;
				g.gap_len=cigar[i] >> BAM_CIGAR_SHIFT;
				position += g.gap_len;
				g.right_start=position+1;

				// Add left cigar data
				g.left_cigar_len=cigar[i-1]>>BAM_CIGAR_SHIFT;
				// Adding one -> always positive values -> output as factor
				g.left_cigar_type=(cigar[i-1] & BAM_CIGAR_MASK)+1;

				// Add right cigar data
				g.right_cigar_len=cigar[i+1]>>BAM_CIGAR_SHIFT;
				// Adding one -> always positive values -> output as factor
				g.right_cigar_type=(cigar[i+1] & BAM_CIGAR_MASK)+1;

				gap_list_push_back(l,g);
			}
			else if(op == BAM_CMATCH || op == BAM_CINS)
			{
				// shift position
				position +=(cigar[i] >> BAM_CIGAR_SHIFT);
			}
		}
		//printf("\n");
	}
}

#endif /* GAP_LIST_H_ */
