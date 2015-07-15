/*
 * gapSiteListList.h
 *
 *  Created on: 11.01.2013
 *      Author: kaisers
 */

#ifndef GAPSITELISTLIST_H_
#define GAPSITELISTLIST_H_

#include "gapSiteList.h"

typedef struct site_ll_element
{
	site_list *l;
	struct site_ll_element *next;
	struct site_ll_element *last;
}site_ll_element;

static R_INLINE site_ll_element * site_ll_element_init(site_list *l)
{
	site_ll_element* el=Calloc(1,site_ll_element);
	el->l=l;
	return el;
}

static R_INLINE void site_ll_element_destroy(site_ll_element *el)
{
	site_list_destroy(el->l);
	Free(el);
}


typedef struct site_ll
{
	struct site_ll_element *first;
	struct site_ll_element *last;
	struct site_ll_element *curr;
	unsigned size;
} site_ll;

site_ll * site_ll_init()
{
	site_ll * sll = Calloc(1, site_ll);
	return sll;
}

void site_ll_add_site_list(site_ll *l,site_list *sl)
{
	site_ll_element *el=site_ll_element_init(sl);
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

void site_ll_destroy(site_ll *l)
{
	if(l==0)
		return;
	site_ll_element *el;
	while(l->size>0)
	{
		el=l->first;
		//printf("[site_ll_destroy] refid %u\n",el->l->refid);
		if(l->size>1)
			l->first=l->first->next;
		--(l->size);
		site_ll_element_destroy(el);
	}
	Free(l);
}

site_list * site_ll_get_curr_site_list_pp(site_ll *l)
{
	if(l->curr==0)
		return 0;
	site_list* sl=l->curr->l;
	l->curr=l->curr->next;
	return sl;
}

void site_ll_set_curr_first(site_ll *l) { l->curr=l->first; }

pos_type sum_ll_sizes(const site_ll *l)
{
	if(l->size==0)
		return 0;

	pos_type size=0;
	unsigned i;
	site_ll_element *el=l->first;
	size=el->l->size;
	for(i=1;i<l->size;++i)
	{
		el=el->next;
		size+=el->l->size;
	}
	return size;
}

sle_type get_nAligns(const site_ll *l)
{
	if(l->size == 0)
		return 0;

	sle_type nAligns = 0LLU;
	unsigned i;
	site_ll_element *el = l->first;
	nAligns = el->l->nAligns;
	for(i=1;i<l->size;++i)
	{
		el = el->next;
		nAligns += el->l->nAligns;
	}
	return nAligns;
}

sle_type get_nAlignGaps(const site_ll *l)
{
	sle_type nAlignGaps = 0LLU;
	unsigned i;
	site_ll_element *el = l->first;
	nAlignGaps = el->l->nAlignGaps;
	for(i=1; i<l->size; ++i)
	{
		el = el->next;
		nAlignGaps += el->l->nAlignGaps;
	}
	return nAlignGaps;
}


#endif /* GAPSITELISTLIST_H_ */
