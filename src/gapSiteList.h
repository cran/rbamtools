/*
 * gapSiteList.h
 *
 *  Created on: 16.12.2012
 *      Author: wolfgang kaisers
 */

#ifndef GAPSITELIST_H_
#define GAPSITELIST_H_

#include <stdlib.h>
#include <R.h>
#include "bitmask.h"
#include "samtools/sam.h"
#include "samtools/bam.h"

typedef unsigned long pos_type;
typedef index_type    idx_type;
typedef bitmap_type   sle_type;

// + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + +
//
// site_list_element
//
// + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + +


typedef struct site_list_element
{
	pos_type   lend;
	pos_type   rstart;
	index_type gap_len;

	// contains maximum of right-cigar-size
	index_type r_cigar_size;

	idx_type nAligns;
	idx_type nProbes;
	sle_type lcl; 		// left    cigar length
	sle_type mcl;           // minimal cigar length

	// List
	struct site_list_element *next;
	struct site_list_element *last;
} site_list_element;


static R_INLINE int site_list_el_compare(site_list_element *lhs,site_list_element *rhs)
{
	if(lhs->lend>rhs->lend || lhs->rstart>rhs->rstart)
		return 1;
	if(lhs->lend<rhs->lend || lhs->rstart<rhs->rstart)
		return -1;

	return 0;
}


static R_INLINE site_list_element* site_list_el_init(pos_type lend, pos_type rstart,index_type gap_len,unsigned char lcl, index_type rcs, unsigned char mcl)
{
	site_list_element *sle=calloc(sizeof(site_list_element),1);
	sle->lend=lend;
	sle->rstart=rstart;
	sle->gap_len=gap_len;
	sle->r_cigar_size=rcs;
	sle->nAligns=1;
	sle->nProbes=1;

	sle->lcl=(sle_type) lcl;
	sle->mcl=(sle_type) mcl;
	return sle;
}

static R_INLINE site_list_element* copy_site_list_element(const site_list_element *in)
{
	site_list_element *el=calloc(sizeof(site_list_element),1);
	el->lend=in->lend;
	el->rstart=in->rstart;
	el->gap_len=in->gap_len;
	el->r_cigar_size=in->r_cigar_size;
	el->nAligns=in->nAligns;
	el->nProbes=in->nProbes;
	el->lcl=in->lcl;
        el->mcl=in->mcl;
	return el;
}

void static R_INLINE site_list_add_cs(struct site_list_element* el,sle_type lcl,index_type rcs,sle_type mcl)
{
	++(el->nAligns);
	// Maximum of right cigar size
	if(el->r_cigar_size<rcs)
		el->r_cigar_size=rcs;
	r_addVal(&(el->lcl),lcl);
	r_addVal_f(&(el->mcl),mcl);
}

// static _inline_ void site_list_el_destroy(site_list_element* el) { free(el); }

// + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + +
//
// site_list
//
// + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + +


// + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + +
// Basic definitions

typedef struct site_list
{
	struct site_list_element *first;
	struct site_list_element *last;
	struct site_list_element *curr;
	pos_type size;

	unsigned refid;
	unsigned long int nAligns;
	unsigned long int nAlignGaps;
} site_list;

site_list* site_list_init(){ return(calloc(sizeof(site_list),1));}

void site_list_destroy(site_list *l)
{
	//printf("[site_list_destroy] refid: %u\n",l->refid);
	if(l==0)
		return;
	site_list_element *el;
	while(l->size>0)
	{
		//printf("[sl_destroy] size %lu\n",l->size);
		el=l->first;
		if(l->size>1)
			l->first=l->first->next;
		--(l->size);
		// Destroy element
		free(el);
	}
	free(l);
}

// + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + +
// Basic element operations

void site_list_set_curr_last (site_list *l) { l->curr=l->last;  }
void site_list_set_curr_first(site_list *l) { l->curr=l->first; }

static R_INLINE void site_list_insert_at_first(site_list *l,pos_type lend,pos_type rstart,index_type gap_len,unsigned char lcl,index_type rcs,unsigned char mcl)
{
	site_list_element *el=site_list_el_init(lend,rstart,gap_len,lcl,rcs,mcl);
	if(l->size==0)
	{
		l->first=el;
		l->last=el;
		l->size=1;
	}
	else
	{
		el->next=l->first;
		l->first->last=el;
		l->first=el;
		++(l->size);
	}
}

static R_INLINE void site_list_insert_at_last(site_list *l,pos_type lend,pos_type rstart,index_type gap_len,unsigned char lcl,index_type rcs,unsigned char mcl)
{
	site_list_element *el=site_list_el_init(lend,rstart,gap_len,lcl,rcs,mcl);
	if(l->size==0)
	{
		l->first=el;
		l->last=el;
		l->size=1;
	}
	else
	{
		el->last=l->last;
		l->last->next=el;
		l->last=el;
		++(l->size);
	}
}

static R_INLINE void site_list_insert_el_at_last(site_list *l, site_list_element *el)
{
	if(l->size==0)
	{
		l->first=el;
		l->last=el;
		l->size=1;
	}
	else
	{
		el->last=l->last;
		l->last->next=el;
		l->last=el;
		++(l->size);
	}
}

static R_INLINE void site_list_copy_to_last(site_list *l, const site_list_element *in)
{
	site_list_element *el=copy_site_list_element(in);
	if(l->size==0)
	{
		l->first=el;
		l->last=el;
		l->size=1;
	}
	else
	{
		el->last=l->last;
		l->last->next=el;
		l->last=el;
		++(l->size);
	}
}

static R_INLINE void site_list_insert_pre_current(site_list *l,pos_type lend,pos_type rstart,index_type gap_len,unsigned char lcl,index_type rcs,unsigned char mcl)
{
	if(l->curr==l->first)
	{
		site_list_insert_at_first(l,lend,rstart,gap_len,lcl,rcs,mcl);
		return;
	}

	site_list_element *el=site_list_el_init(lend,rstart,gap_len,lcl,rcs,mcl);
	el->last=l->curr->last;
	el->next=l->curr;
	l->curr->last->next=el;
	l->curr->last=el;
	++(l->size);
}

static R_INLINE void site_list_insert_post_current(site_list *l,pos_type lend,pos_type rstart,index_type gap_len,unsigned char lcl,index_type rcs,unsigned char mcl)
{
	if(l->curr==l->last)
	{
		site_list_insert_at_last(l,lend,rstart,gap_len,lcl,rcs,mcl);
		return;
	}
	site_list_element *el=site_list_el_init(lend,rstart,gap_len,lcl,rcs,mcl);
	el->last=l->curr;
	el->next=l->curr->next;
	l->curr->next=el;
	el->next->last=el;
	++(l->size);
}

site_list_element* site_list_get_curr_mm(site_list *l)
{
	site_list_element *el=l->curr;
	if(l->curr!=NULL)
		l->curr=l->curr->last;
	return el;
}

site_list_element* site_list_get_curr_pp(site_list *l)
{
	site_list_element *el=l->curr;
	if(l->curr!=NULL)
		l->curr=l->curr->next;
	return el;
}

// + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + +
// Higher level functions

void site_list_add_element(site_list *l,pos_type lend,pos_type rstart,index_type gap_len,unsigned char lcl,index_type rcs,unsigned char mcl)
{
	//printf("[site_list_add_element] lcl: %u\trcs: %u\n",lcl,rcs);
	site_list_element *el;
	// + + + Empty List + + +
	if(l->size==0)
	{
		el=site_list_el_init(lend,rstart,gap_len,lcl,rcs,mcl);
		l->first=el;
		l->last=el;
		l->size=1;
		return;
	}

	// + + + New Element is last element
	if(lend>l->last->lend)
	{
		site_list_insert_at_last(l,lend,rstart,gap_len,lcl,rcs,mcl);
		return;
	}

	// + + + One Element in list + + +
	if(l->size==1)
	{
		if(lend==l->first->lend)
		{
			if(rstart>l->first->rstart)
			{
				site_list_insert_at_last(l,lend,rstart,gap_len,lcl,rcs,mcl);
				return;
			}
			if(rstart==l->first->rstart)
			{
				site_list_add_cs(l->first,lcl,rcs,mcl);
				return;
			}
		}
		if(lend>l->first->lend)
		{
			site_list_insert_at_last(l,lend,rstart,gap_len,lcl,rcs,mcl);
			return;
		}
		site_list_insert_at_first(l,lend,rstart,gap_len,lcl,rcs,mcl);
		return;
	}

	// Start searching right position in list:
	l->curr=l->last;

	// walk back until lend >= curr->lend and rstart>=l->curr->rstart
	while(lend<l->curr->lend && l->curr!=l->first)
		l->curr=l->curr->last;
	while(lend==l->curr->lend && rstart<l->curr->rstart && l->curr!=l->first)
		l->curr=l->curr->last;

	el=l->curr;
	//printf("[site_list_add_element] Start search. curr: lend=%lu\trstart%lu\n",el->lend,el->rstart);

	if(lend>el->lend)
	{
		site_list_insert_post_current(l,lend,rstart,gap_len,lcl,rcs,mcl);
		return;
	}
	if(lend==el->lend)
	{

		if(rstart>el->rstart)
		{
			site_list_insert_post_current(l,lend,rstart,gap_len,lcl,rcs,mcl);
			return;
		}
		if(rstart==el->rstart)
		{
			// No new site is inserted but instead information
			// added to existing sitesiteList_add_llcl
			site_list_add_cs(l->curr,lcl,rcs,mcl);
			return;
		}
	}
	// lend < el->lend or rstart < el->rstart
	site_list_insert_pre_current(l,lend,rstart,gap_len,lcl,rcs,mcl);
}

// args can't be const because *curr* values are intermediately changed
site_list* site_list_merge(site_list *lhs, site_list *rhs,unsigned refid)
{
	site_list_element *lcurr,*rcurr,*ins;
	lcurr=lhs->curr;
	rcurr=rhs->curr;
	lhs->curr=lhs->first;
	rhs->curr=rhs->first;

	site_list_element *le=site_list_get_curr_pp(lhs);
	site_list_element *re=site_list_get_curr_pp(rhs);

	site_list *res=site_list_init();
	res->refid=refid;

	while(le!=NULL || re !=NULL)
	{
		if(le==NULL)
		{
			while(re!=NULL)
			{
				site_list_copy_to_last(res,re);
				re=site_list_get_curr_pp(rhs);
			}
		}
		else if(re==NULL)
		{
			while(le!=NULL)
			{
				site_list_copy_to_last(res,le);
				le=site_list_get_curr_pp(lhs);
			}
		}
		else
		{
			if(le->lend<re->lend)
			{
				//printf("[site_list_merge] ++le-le<re-le\n");
				site_list_copy_to_last(res,le);
				le=site_list_get_curr_pp(lhs);
			}
			else if(re->lend<le->lend)
			{
				site_list_copy_to_last(res,re);
				re=site_list_get_curr_pp(rhs);
			}
			// le->lend==re->lend
			else if(le->rstart<re->rstart)
			{
				site_list_copy_to_last(res,le);
				le=site_list_get_curr_pp(lhs);
			}
			else if(re->rstart<le->rstart)
			{
				site_list_copy_to_last(res,re);
				re=site_list_get_curr_pp(rhs);
			}
			else
			{
				// le->lend == re->lend && le->rstart == re->rstart
				// This is the important part
				// where the code departs from standard merging:
				ins=copy_site_list_element(le);
				r_zip(le->lcl,re->lcl,&(ins->lcl));
				r_zip_f(le->mcl,re->mcl,&(ins->mcl));
				ins->nProbes+=re->nProbes;
				ins->nAligns+=re->nAligns;
				site_list_insert_el_at_last(res,ins);
				le=site_list_get_curr_pp(lhs);
				re=site_list_get_curr_pp(rhs);
			}
		}
	}
	lhs->curr=lcurr;
	rhs->curr=rcurr;
	return res;
}

idx_type* site_list_get_max_match(site_list *l)
{
	idx_type *res=calloc(sizeof(idx_type),l->size);
	unsigned i=0;
	site_list_element *el,*curr;

	curr=l->curr;
	site_list_set_curr_first(l);

	for(i=0;i<l->size;++i)
	{
		el=site_list_get_curr_pp(l);
		res[i]=(idx_type) ((el->mcl)>>idx[0])&0xFF; //getByte(el->mcl,0);
	}
	// Reset curr
	l->curr=curr;
	return res;
}

idx_type* site_list_get_n_lpos(site_list *l)
{
	idx_type *res=calloc(sizeof(idx_type),l->size);
	unsigned i=0;
	site_list_element *el,*curr;

	curr=l->curr;
	site_list_set_curr_first(l);

	for(i=0;i<l->size;++i)
	{
		el=site_list_get_curr_pp(l);
		res[i]=bitmask_nPos(el->lcl);
	}
	// Reset curr
	l->curr=curr;
	return res;

}

idx_type* site_list_get_n_match(site_list *l)
{
	idx_type *res=calloc(sizeof(idx_type),l->size);
	unsigned i=0;
	site_list_element *el,*curr;

	curr=l->curr;
	site_list_set_curr_first(l);
	for(i=0;i<l->size;++i)
	{
		el=site_list_get_curr_pp(l);
		res[i]=bitmask_nPos(el->mcl);
	}
	// Reset curr
	l->curr=curr;
	return res;
}

idx_type* site_list_get_sum_match(site_list *l)
{
	idx_type *res=calloc(sizeof(idx_type),l->size);
	unsigned i=0;
	site_list_element *el,*curr;

	curr=l->curr;
	site_list_set_curr_first(l);

	for(i=0;i<l->size;++i)
	{
		el=site_list_get_curr_pp(l);
		res[i]=bitmask_sumPos(el->mcl);
	}
	// Reset curr
	l->curr=curr;
	return res;
}


site_list* site_list_copy(site_list *l)
{
	site_list* res=site_list_init();
	res->nAligns=l->nAligns;
	res->nAlignGaps=l->nAlignGaps;
	res->refid=l->refid;

	site_list_element *et,*el=l->first;
	unsigned i;
	for(i=0;i<l->size;++i)
	{
		et=copy_site_list_element(el);
		site_list_insert_el_at_last(res,et);
		el=el->next;
	}
	return res;
}


// + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + +
//
// list_ext_gaps: Function for listing BAM-align-gaps
//
// + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + +

void list_gap_sites(site_list *l,const bam1_t* align)
{
	// Count total number of Aligns
	++(l->nAligns);

	// Ungapped align -> nothing to do
	if(align->core.n_cigar==1)
		return;

	// + + + + + + + + + + + + + + + + + + + + + + + + + + + + + +
	// Assumes that first and last cigar type cannot be N
	// Assumes that left  and right of a N there's always a M.
	// + + + + + + + + + + + + + + + + + + + + + + + + + + + + + +
	uint32_t position,n_cigar,lend,rstart,gap_len,i,left_len, right_len, min_mtc;
	uint32_t *cigar;
	int op;

	cigar    = bam1_cigar(align);
	// Add size of first cigar to position to get 1-based left stop (BAM is 0-based)
	position = align->core.pos + (cigar[0] >> BAM_CIGAR_SHIFT);
	n_cigar  = align->core.n_cigar;
	op = cigar[0] & BAM_CIGAR_MASK;
	if((op!= BAM_CMATCH))
		return;	   // ERROR: First CIGAR item cannot be I or D or N.

	//Rprintf("[list_gap_sites] refid: %i\tposition: %i\tn_cigar: %i\n",align->core.tid,align->core.pos,n_cigar);


	if(n_cigar>2)
	{
		// In this loop there's always a valid left and right cigar

		// BAM_CMATCH		= 0 = M match or mismatch
		// BAM_CREF_SKIP	= 3 = N skip on the reference (e.g. spliced alignment) *
		// BAM_CINS		    = 1 = S clip on the read with clipped sequence

		// First cigar item must be M! Otherwise drop align.
		if((cigar[0] & BAM_CIGAR_MASK)!=BAM_CMATCH)
			return;

		// Add size of first cigar to position to get *1-based* left stop
		// BAM is *0-based*
		// position always points to *1-based* last nucleotide of active cigar item
		position = align->core.pos + (cigar[0] >> BAM_CIGAR_SHIFT);

		for(i=1;i<(n_cigar-1);++i)
		{

			op = cigar[i] & BAM_CIGAR_MASK;
			//printf("op %i\t",op);
			if(op==BAM_CREF_SKIP)
			{
				// Count gaps
				++(l->nAlignGaps);
				// N -> Add gap to list
				lend=position;
				gap_len=cigar[i] >> BAM_CIGAR_SHIFT;

				position += gap_len;
				// position now points to *1-based* last nuc of gap

				// Right start (first nuc) is one after last intron position
				rstart=position+1;

				// left side of N must be M otherwise skip align
				if((cigar[i-1] & BAM_CIGAR_MASK)!=BAM_CMATCH)
					return;
				// Add left cigar data
				left_len=cigar[i-1]>>BAM_CIGAR_SHIFT;
				// Adding one -> always positive values -> output as factor
				//g.left_cigar_type=BAM_CMATCH+1;

				// right side of N must be M otherwise skip align
				if((cigar[i+1] & BAM_CIGAR_MASK)!=BAM_CMATCH)
					return;
				// Add right cigar data
				right_len=cigar[i+1]>>BAM_CIGAR_SHIFT;
				// Adding one -> always positive values -> output as factor
				// g.right_cigar_type=BAM_CMATCH+1;

				min_mtc=left_len>right_len ? right_len : left_len;
				site_list_add_element(
						l,								// site_list
						lend,							// lend
						rstart,							// rstart
						gap_len,						// gap_len
						left_len,						// lcl
						right_len,						// rcs
						min_mtc                         // mcl
						);

				// next cigar can also be processed because we know it's an M
				++i;
				position+=(cigar[i] >> BAM_CIGAR_SHIFT);
				// position now points to *1-based* last nuc of match
			}
			else if(op == BAM_CMATCH)
			{
				// position then points on rightmost nuc of exon
				position +=(cigar[i] >> BAM_CIGAR_SHIFT);
			}
			else if(op == BAM_CDEL)
			{
				// shift position
				position +=(cigar[i] >> BAM_CIGAR_SHIFT);
			}
		}
		//printf("\n");
	}
}


#endif /* GAPSITELIST_H_ */
