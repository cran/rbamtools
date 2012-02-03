/*
 * align_list.h
 *
 *  Created on: 25.01.2012
 *      Author: wolfgang
 */

#ifndef ALIGN_LIST_H_
#define ALIGN_LIST_H_
#include "samtools/sam.h"
#include "samtools/bam.h"

///////////////////////////////////////////////////////////////////////////////////////////////////
// basic definitions

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
} align_list;

///////////////////////////////////////////////////////////////////////////////////////////////////
// basic definitions

inline align_list * init_align_list()
{
	return (align_list*) calloc(1,sizeof(align_list));
}

inline void copy_align(bam1_t *target,const bam1_t * const source)
{
	// see bam.h duplicate_align
	if(target==NULL)
		return;
	*target=*source;
	target->m_data=source->data_len;
	free(target->data);
	target->data=(uint8_t*)calloc(target->data_len,1);
	memcpy(target->data,source->data,target->data_len);
}

inline bam1_t *duplicate_align(const bam1_t *src)
{
	bam1_t *b;
	b = bam_init1();
	*b = *src;
	b->m_data = b->data_len;
	b->data = (uint8_t*)calloc(b->data_len, 1);
	memcpy(b->data, src->data, b->data_len);
	return b;
}

inline align_element* init_elem(const bam1_t *align)
{
	align_element *e=calloc(1,sizeof(align_element));
	e->align=duplicate_align(align);
	return e;
}


///////////////////////////////////////////////////////////////////////////////////////////////////
// list generic accessor functions

inline void push_back(align_list *l, const bam1_t *align)
{
	align_element *e=init_elem(align);
	if(l->size==0)	// list is empty
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
	//printf("push_back pos: %i\n",align->core.pos);
}

inline void push_front(align_list *l,const bam1_t *align)
{
	align_element *e=init_elem(align);
	if(l->first_el==0)	// list is empty
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

inline void pop_back(align_list *l)
{
	if(l->first_el!=l->last_el)
	{
		align_element *e=l->last_el;
		e->last_el->next_el=0;
		l->last_el=e->last_el;
		bam_destroy1(e->align);
		free(e);
		--(l->size);
	}
	else if(l->last_el>0)
	{
		//printf("pop_back pos: %i\n",l->first_el->align->core.pos);
		bam_destroy1(l->first_el->align);
		free(l->first_el);
		l->first_el=0;
		l->last_el=0;
		l->size=0;
	}

}

inline void pop_front(align_list *l)
{
	//printf("pop_front pos %i\n",l->first_el->align->core.pos);
	if(l->first_el!=l->last_el)
	{
		align_element *e=l->first_el;
		e->next_el->last_el=0;
		l->first_el=e->next_el;
		bam_destroy1(e->align);
		free(e);
		--(l->size);
	}
	else if(l->first_el>0)
	{
		bam_destroy1(l->first_el->align);
		free(l->first_el);
		l->first_el=0;
		l->last_el=0;
		l->size=0;
	}
}

inline void wind_back(align_list *l)
{
	l->curr_el=NULL;
	return;
}


///////////////////////////////////////////////////////////////////////////////////////////////////
// higher level convenience functions

inline void al_destroy(align_list *l)
{
	while(l->size>0)
		pop_front(l);
	free(l);
}
inline bam1_t * get_next_align(align_list *l)		// Returns A COPY of current align
{
	if(l->first_el==NULL)
		return (bam1_t*) NULL;

	if(l->curr_el==NULL)
	{
		l->curr_el=l->first_el;
		return duplicate_align(l->curr_el->align);	// Copy!
	}
	if(l->curr_el->next_el==NULL)
		return (bam1_t*) NULL;

	l->curr_el=l->curr_el->next_el;
	return duplicate_align(l->curr_el->align);		// Copy!
}

inline bam1_t * get_prev_align(align_list *l)
{
	if((l->last_el)==NULL)
		return (bam1_t*) NULL;
	if((l->curr_el)==NULL)
	{
		l->curr_el=l->last_el;
		return duplicate_align(l->curr_el->align);	// Copy!
	}
	if((l->curr_el->last_el)==NULL)
		return (bam1_t*) NULL;
	l->curr_el=l->curr_el->last_el;
	return duplicate_align(l->curr_el->align);
}

inline void write_current_align(align_list *l,bam1_t *align)
{
	if((l->curr_el)!=NULL)
		copy_align(l->curr_el->align,align);
	return;
}

inline void pp_curr_align(align_list *l)
{
	if((l->curr_el)==NULL)
	{
		l->curr_el=l->first_el;
		return;
	}
	l->curr_el=(l->curr_el->next_el);
}

inline void mm_curr_align(align_list *l)
{
	if((l->curr_el)==NULL)
	{
		l->curr_el=l->last_el;
		return;
	}
	l->curr_el=(l->curr_el->last_el);
}

inline void insert_past_curr_align(align_list *l,bam1_t *align)
{
	align_element *e=init_elem(align);
	if(l->first_el==NULL)
	{
		l->first_el=e;
		l->last_el=e;
		l->size=1;
		return;
	}
	if(l->curr_el==NULL)
	{
		l->first_el->last_el=e;
		e->next_el=e;
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

inline void insert_pre_curr_align(align_list *l,bam1_t *align)
{
	align_element *e=init_elem(align);
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


#endif /* ALIGN_LIST_H_ */
