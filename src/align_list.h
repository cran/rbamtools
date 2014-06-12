/*
 *	File		: align_list.h
 *	Content		: Double linked list which contains bam1_t align structs
 *
 * 	Created on	: 25.01.2012
 *      Author		: Wolfgang Kaisers
 *
 *	Changelog	:
 *			01.Nov.12 [get_const_next_align] Function added (returns align without copying).
 */

#ifndef ALIGN_LIST_H_
#define ALIGN_LIST_H_
#include "samtools/sam.h"
#include "samtools/bam.h"
#include "samtools/rdef.h"

/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 * basic definitions
 */


typedef struct align_element
{
	bam1_t *align;
	struct align_element *last_el;
	struct align_element *next_el;
} align_element;

typedef struct {
	align_element *first_el;
	align_element *last_el;
	align_element *curr_el;
	unsigned long size;
	int32_t min_seqlen;		/* bam1_core_t->l_qseq */
	int32_t max_seqlen;

	/* 0-based BAM coordinates */
	unsigned      seqid;
	unsigned long range_begin;
	unsigned long range_end;

	/* BamHeader values */
	char * refname;
	unsigned long seq_LN;

	unsigned complex;
} align_list;

/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 * basic functions
 */

align_list * init_align_list()
{
	align_list * l=(align_list*) calloc(1,sizeof(align_list));
	/* Create large number because min_seqlen must decrease	*/
	l->min_seqlen=10000;
	return l;
}

// Copy of bam_copy1 in bam.h
static R_INLINE void copy_align(bam1_t *target,const bam1_t * const source)
{
	/* see bam.h duplicate_align	*/
	if(target==NULL)
		return;
	*target=*source;
	target->m_data=source->data_len;
	free(target->data);
	target->data=(uint8_t*)calloc((size_t)(target->data_len),1);
	memcpy(target->data,source->data,(size_t)target->data_len);

#ifdef BAM1_ADD_CIGAR
	// + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + //
	// Fill additional cigar field with copy of cigar data
	// + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + //
	// double free!
	target->cigar=0;
	COPY_CIGAR_VALUES(target);
	// + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + //
#endif
}

// Copy of bam_dup1 in bam.h
static R_INLINE bam1_t *duplicate_align(const bam1_t *src)
{
	bam1_t *b;
	b = bam_init1();

	// Static copy (replaces *b=*src)
	b->core=src->core;
	b->l_aux=src->l_aux;
	b->data_len=src->data_len;
	b->m_data = b->data_len;	// max data len=data_len (new allocated)
	// b->cigar must stay 0!

	b->data = (uint8_t*)calloc((size_t)(b->data_len), 1);
	memcpy(b->data, src->data, (size_t)(b->data_len));

#ifdef BAM1_ADD_CIGAR
	// + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + //
	// Fill additional cigar field with copy of cigar data
	// + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + //
	COPY_CIGAR_VALUES(b);
	// + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + //
#endif

	return b;
}

static R_INLINE align_element* align_list_init_elem(const bam1_t *align)
{
	align_element *e=calloc(1,sizeof(align_element));
	e->align=duplicate_align(align);
	return e;
}

static R_INLINE void align_list_destroy_elem(align_element *e)
{
	//Rprintf("[align_list_destroy_elem] e=%u\n",e);
	bam_destroy1(e->align);
	free(e);
}


/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 * list generic accessor functions
 */

void align_list_push_back(align_list *l, const bam1_t *align)
{
	align_element *e;

	if((l->max_seqlen)<(align->core.l_qseq))
			l->max_seqlen=align->core.l_qseq;

	if((l->min_seqlen)>(align->core.l_qseq) && l->size>1)
		l->min_seqlen=align->core.l_qseq;

	e=align_list_init_elem(align);
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
	//Rprintf("[align_list_push_back] e=%x\n",e);
}

void align_list_push_front(align_list *l,const bam1_t *align)
{
	align_element *e;

	if((l->max_seqlen)<(align->core.l_qseq))
			l->max_seqlen=align->core.l_qseq;

	if((l->min_seqlen)>(align->core.l_qseq) && l->size>1)
		l->min_seqlen=align->core.l_qseq;

	e=align_list_init_elem(align);
	if(l->first_el==0)
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
	//Rprintf("[align_list_push_front] e=%x\n",e);
}

void align_list_pop_back(align_list *l)
{
	//Rprintf("[align_list_pop_back] e=%x\n",l->last_el);
	if(l->first_el!=l->last_el)
	{
		align_element *e=l->last_el;
		e->last_el->next_el=0;
		l->last_el=e->last_el;
		--(l->size);
		align_list_destroy_elem(e);
	}
	else if(l->last_el!=0)
	{
		align_list_destroy_elem(l->last_el);
		l->first_el=0;
		l->last_el=0;
		l->size=0;
	}

}

void align_list_pop_front(align_list *l)
{
	//Rprintf("[align_list_pop_front] e=%x\n",l->first_el);
	align_element *e;
	if(l->first_el!=l->last_el)
	{
		e=l->first_el;
		e->next_el->last_el=0;
		l->first_el=e->next_el;
		--(l->size);
		align_list_destroy_elem(e);
	}
	else if(l->first_el!=0)
	{
		align_list_destroy_elem(l->first_el);
		l->first_el=0;
		l->last_el=0;
		l->size=0;
	}
}

void wind_back(align_list *l)
{
	l->curr_el=NULL;
	return;
}


void align_list_mv_curr_elem(align_list *src,align_list *target)
{
	/*
	 * Moves current element from src to end of target list
	 * and moves curr_el pointer to next align
	 *
	 */

	align_element *e;

	if(src->curr_el==NULL)
		return;

	e=src->curr_el;
	src->curr_el=e->next_el;

	/*
	 * Remove e from src list
	 */
	if(e->next_el!=NULL)
		e->next_el->last_el=e->last_el;
	if(e->last_el!=NULL)
		e->last_el->next_el=e->next_el;

	/*
	 * Insert e into end of target list
	 */
	if(target->size==0)
	{
		target->first_el=e;
		target->last_el=e;
		target->size=1;
	}
	else
	{
		target->last_el->next_el=e;
		e->last_el=target->last_el;
		target->last_el=e;
		++(target->size);
	}
}

/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 * higher level convenience functions
 */

void destroy_align_list(align_list *l)
{
	while(l->size>0)
		align_list_pop_front(l);

	free(l->refname);
	free(l);
}

/*
 * Returns a *COPY* of current align
 */
bam1_t * get_next_align(align_list *l)
{

	if(l->first_el==NULL)
	{
		return (bam1_t*) NULL;
	}

	if(l->curr_el==NULL)
	{
		l->curr_el=l->first_el;
		return duplicate_align(l->curr_el->align);	/* Copy */
	}
	if(l->curr_el->next_el==NULL)
	{
		l->curr_el=NULL;
		return (bam1_t*) NULL;
	}

	l->curr_el=l->curr_el->next_el;
	return duplicate_align(l->curr_el->align);		/* Copy */
}

const bam1_t * get_const_next_align(align_list *l)		/* Returns a *CONSTANT REFERENCE* to current align */
{

	if(l->first_el==NULL)
	{
		return (bam1_t*) NULL;
	}

	if(l->curr_el==NULL)
	{
		l->curr_el=l->first_el;
		return l->curr_el->align;	/* No Copy! */
	}
	if(l->curr_el->next_el==NULL)
	{
		l->curr_el=NULL;
		return (bam1_t*) NULL;
	}

	l->curr_el=l->curr_el->next_el;
	return l->curr_el->align;		/* No copy! */
}


bam1_t * get_prev_align(align_list *l)
{
	if((l->last_el)==NULL)
		return (bam1_t*) NULL;
	if((l->curr_el)==NULL)
	{
		l->curr_el=l->last_el;
		return duplicate_align(l->curr_el->align);	/* Copy! */
	}
	if((l->curr_el->last_el)==NULL)
	{
		l->curr_el=NULL;
		return (bam1_t*) NULL;
	}
	l->curr_el=l->curr_el->last_el;
	return duplicate_align(l->curr_el->align);
}

void write_current_align(align_list *l,bam1_t *align)
{
	if((l->curr_el)!=NULL)
		copy_align(l->curr_el->align,align);
	return;
}

void pp_curr_align(align_list *l)
{
	if((l->curr_el)==NULL)
	{
		l->curr_el=l->first_el;
		return;
	}
	l->curr_el=(l->curr_el->next_el);
}

void mm_curr_align(align_list *l)
{
	if((l->curr_el)==NULL)
	{
		l->curr_el=l->last_el;
		return;
	}
	l->curr_el=(l->curr_el->last_el);
}

void insert_past_curr_align(align_list *l,bam1_t *align)
{
	align_element *e;

	if((l->max_seqlen)<(align->core.l_qseq))
			l->max_seqlen=align->core.l_qseq;

	if((l->min_seqlen)>(align->core.l_qseq) && l->size>1)
		l->min_seqlen=align->core.l_qseq;

	e=align_list_init_elem(align);
	if(l->first_el==NULL)
	{
		l->first_el=e;
		l->last_el=e;
		(l->size)=1;
		return;
	}
	if(l->curr_el==NULL)
	{
		l->first_el->last_el=e;
		e->next_el=l->first_el;
		l->first_el=e;
		++(l->size);
		return;
	}
	if(l->curr_el->next_el==NULL)
	{
		l->curr_el->next_el=e;
		e->last_el=l->curr_el;
		l->last_el=e;
		++(l->size);
		return;
	}
	e->next_el=l->curr_el->next_el;
	e->last_el=l->curr_el;
	l->curr_el->next_el->last_el=e;
	l->curr_el->next_el=e;
	++(l->size);

}

void insert_pre_curr_align(align_list *l,bam1_t *align)
{
	align_element *e;
	if((l->max_seqlen)<(align->core.l_qseq))
			l->max_seqlen=align->core.l_qseq;

	if((l->min_seqlen)>(align->core.l_qseq) && l->size>1)
		l->min_seqlen=align->core.l_qseq;

	e=align_list_init_elem(align);
	if(l->first_el==NULL)
	{
		l->first_el=e;
		l->last_el=e;
		l->size=1;
		return;
	}
	if(l->curr_el==NULL)
	{
		l->last_el->next_el=e;
		e->last_el=l->last_el;
		l->last_el=e;
		++(l->size);
		return;
	}
	if(l->curr_el->last_el==NULL)
	{
		l->curr_el->last_el=e;
		e->next_el=l->curr_el;
		l->first_el=e;
		++(l->size);
		return;
	}
	e->last_el=l->curr_el->last_el;
	e->next_el=l->curr_el;
	l->curr_el->last_el->next_el=e;
	l->curr_el->last_el=e;
	++(l->size);
}

static R_INLINE void add_match_depth(unsigned long  *ald, long begin, long end, long position, uint32_t cigar_len)
{
	/*
	 * 0-based index of last count value
	 * nPos=range_end+1 (size of count)
	 */
	long range_end, w_start, w_end,i;
	long align_end;

	range_end=((long)end)-((long)begin);
	if(range_end<1)
		return;

	/*
	 * position = 0-based align_begin
	 * align_end = 0-based
	 */
	align_end=position+cigar_len-1;

	/*
	 * first and last writing index (0-based). <0 allowed
	 */
	w_start= position - begin;
	w_end  = align_end- begin;

	if((w_start>=range_end) | (w_end<0))
		return;

	/*
	 * secure array boundaries (also all trimming)
	 */
	w_start=(w_start<0)       ? 0         : w_start;
	w_end  =(w_end>range_end) ? range_end : w_end;

	/* Do counting */
	for(i=w_start;i<=w_end;++i)
		++(ald[i]);

	return;
}

void count_align_depth (unsigned long *ald,unsigned long begin,unsigned long end,const bam1_t * align)
{
	uint32_t n_cigar,i;
	int32_t position;
	int op;
	const uint32_t* cigar;


	/*
	 * All positions are 0-based handled.
	 * position: 0-based position of first cigar-op nuc
	 * pos += cigar_shift: shifts to next (first cigar-op nuc)
	 */
	if(!align)
		return;

	cigar=bam1_cigar(align);
	position=align->core.pos;
	n_cigar=align->core.n_cigar;

	/*
	 * Add first cigar (must be match)
	 *
	 */
	add_match_depth(ald,(long) begin,(long) end,position,cigar[0]>>BAM_CIGAR_SHIFT);
	position += (int32_t) (cigar[0] >> BAM_CIGAR_SHIFT);

	if(n_cigar>1)
	{
		/* n_cigar>2 */
		for(i=1;i<(n_cigar-1);++i)
		{
			op = cigar[i] & BAM_CIGAR_MASK;
			if((op==BAM_CREF_SKIP) | (op == BAM_CDEL))
			{
				/*
				 * N or D -> shift position
				 * Rprintf("[count_align_depth] + + NorD + + pos: %lu\tlen: %u\n",position,cigar[i]>>BAM_CIGAR_SHIFT);
				 */
				position += (int32_t) (cigar[i] >> BAM_CIGAR_SHIFT);

			}
			else if(op == BAM_CMATCH)
			{
				/*
				 * M en=cigar[i+1]>>BAM_CIGAR_SHIFT;
				 * Rprintf("[count_align_depth] pos: %lu\tlen: %u\n",position,cigar[i]>>BAM_CIGAR_SHIFT);
				 */
				add_match_depth(ald,(long) begin,(long) end,position,cigar[i]>>BAM_CIGAR_SHIFT);
				/* position then points on rightmost nuc of exon */
				position += (int32_t) (cigar[i] >> BAM_CIGAR_SHIFT);
			}
			/* I: Do nothing */
		}
		/*
		 * Add last cigar (must be match)
		 * Rprintf("[count_align_depth] pos: %lu\tlen: %u\n",position,cigar[i]>>BAM_CIGAR_SHIFT);
		 */
		add_match_depth(ald,(long) begin,(long) end,position,cigar[i]>>BAM_CIGAR_SHIFT);
	}
}

void count_align_gap_depth (unsigned long *ald,unsigned long begin, unsigned long end,const bam1_t * align)
{
	uint32_t n_cigar,i;
	int32_t position, right_cigar_pos;
	int op;
	const uint32_t *cigar;
	uint32_t right_cigar_len;

	/*
	 * All positions are 0-based handled.
	 */
	if(!align)
		return;


	cigar=bam1_cigar(align);
	position=align->core.pos;
	n_cigar=align->core.n_cigar;
	/*
	 * Store count coords for right adjacent match
	 * right_cigar_len==0 says that no cigar op is
	 * still to be counted
	 */
	right_cigar_len=0;


	/*
	 * Always count Matches on left side of N
	 */
	if(n_cigar>2)
	{

		/* Shift position for first cigar op (must be M) */
		position += (int32_t)(cigar[0] >> BAM_CIGAR_SHIFT);
		for(i=1;i<(n_cigar-1);++i)
		{
			/* There always is a left and right cigar */
			op = cigar[i] & BAM_CIGAR_MASK;
			if(op==BAM_CREF_SKIP)
			{
				/* Count left adjacent Match */
				add_match_depth(ald,(long)begin,(long) end,position,cigar[i-1]>>BAM_CIGAR_SHIFT);
				/* shift position */
				position += (int32_t) (cigar[i] >> BAM_CIGAR_SHIFT);
				/* Save position and length for right match */
				right_cigar_pos=position;
				right_cigar_len=(cigar[i+1] >> BAM_CIGAR_SHIFT);
			}
			else if(op == BAM_CDEL)
			{
				position += (int32_t) (cigar[i] >> BAM_CIGAR_SHIFT);

				/*
				 * When D lies behind right adjacent
				 * cigar match of a skip (N):
				 * count saved region and reset.
				 */
				if(right_cigar_len)
				{
					add_match_depth(ald,(long)begin,(long)end,right_cigar_pos,right_cigar_len);
					right_cigar_len=0;
				}
			}
			else if(op == BAM_CINS)
			{
				/*
				 * When I lies behind right adjacent
				 * cigar match of a skip (N):
				 * count saved region and reset.
				 */
				if(right_cigar_len)
				{
					add_match_depth(ald,(long)begin,(long)end,right_cigar_pos,right_cigar_len);
					right_cigar_len=0;
				}
			}
			else if(op == BAM_CMATCH) /* M */
			{
				position += (int32_t)(cigar[i] >> BAM_CIGAR_SHIFT);
			}
		}
		/* Add last unsaved match region */
		if(right_cigar_len)
			add_match_depth(ald,(long)begin,(long) end,right_cigar_pos,right_cigar_len);
	}
	return;
}

#endif /* ALIGN_LIST_H_ */
