/*
 * rdef.h
 *
 *  Created on: 30.05.2014
 *      Author: wolfgang
 */

#ifndef RDEF_H_
#define RDEF_H_


/* * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 *
 * Turn R definition on and off
 *
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
#define R_CRAN

#ifdef R_CRAN
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 * Library is compiled under R
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
#include <R.h>
#else
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 *  Library is compiled without R
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
#define R_INLINE inline
#define Rprintf printf
/*
 * Variadic macros
 * https://gcc.gnu.org/onlinedocs/cpp/Variadic-Macros.html
 */
#define error(...)	\
	printf(__VA_ARGS__);	\
	exit (EXIT_FAILURE);
#endif /* R_CRAN */



/* * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 * Change bam1_t related code in order to
 * correct misaligns
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
#define BAM1_ADD_CIGAR


#ifdef BAM1_ADD_CIGAR
// bam.h: bam_copy1 (796), bam_dup1 (821)
// align_list.h: copy_align, duplicate_align
#define COPY_CIGAR_VALUES(b) 																				\
		do																									\
		{																									\
				free((b)->cigar);																			\
				(b)->cigar=calloc((b)->core.n_cigar,sizeof(uint32_t));										\
				memcpy((b)->cigar,((b)->data + (b)->core.l_qname),(b)->core.n_cigar*sizeof(uint32_t));		\
		}																									\
		while(0)

#endif /*BAM1_ADD_CIGAR		*/


#endif /* RDEF_H_           */
