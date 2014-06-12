/*
 *  File      	: bitmask.h
 *  Created on	: 22.12.2012
 *  Author		: Wolfgang Kaisers
 *  Content		: Usage of 32(64) bit Variables as container for 4(8) sorted Byte-Values
 *              l-version: left  adjusted, descending ordered values
 *              r-version: right adjusted, ascending  ordered values
 *
 *	Changelog 	: 05.02.2013 Corrected bitmask_nPos (max_bitmap_index to bitmap_size) so that all fields are counted.
 *				: 06.02.2013 Added detection of 64bit size and getQmean function.
 *				: 25.11.2013 Corrected r_addVal and r_addVal_f (CRAN Address sanitizer complained)
 */

#ifndef BITMASK_H_
#define BITMASK_H_

#include "samtools/rdef.h"
#ifdef R_CRAN
#include <R.h>
#else
#define R_INLINE inline
#define Rprintf printf
#endif

#include <stdint.h>

// This checks whether this is compiled in 64 bit
#if UINTPTR_MAX == 0xffffffffffffffff
#define BM_64
#endif


#ifdef BM_64
// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// This results in storing values in 64-bit (unsigned) integers:
// 8 byte-sized fields are used
// (Also used by gapSiteList)
typedef unsigned long long  bitmap_type; 							// Size at least 64 bits (C99)
const unsigned long long idx[]      ={0,8,16,24,32,40,48,56};
const unsigned long  pat[]={
                           0xFF             , 0xFF00              , 0xFF0000          , 0xFF000000        ,
                           0xFF00000000     , 0xFF0000000000      , 0xFF000000000000  , 0xFF00000000000000
                           };
const unsigned long lpat[]={
		                    0xFF            , 0xFFFF              , 0xFFFFFF          , 0xFFFFFFFF        ,
		                    0xFFFFFFFFFF    , 0xFFFFFFFFFFFF      , 0xFFFFFFFFFFFFFF  , 0xFFFFFFFFFFFFFFFF
                           };
const unsigned long rpat[]={
		                    0xFFFFFFFFFFFFFFFF, 0xFFFFFFFFFFFFFF00, 0xFFFFFFFFFFFF0000, 0xFFFFFFFFFF000000,
							0xFFFFFFFF00000000, 0xFFFFFF0000000000, 0xFFFF000000000000, 0xFF00000000000000
                           };

const unsigned bitmap_size=8;
const int max_bitmap_index=7;
// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#else
// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// This results in storing values in (at least) 32-bit (unsigned) integers
// 4 byte-sized fields are used
// (Also used by gapSiteList)
typedef unsigned long       bitmap_type;							// Size at least 32 bits (C99)
const unsigned long idx[]      ={0,8,16,24};

const unsigned long  pat[]={
                           0xFF             , 0xFF00              , 0xFF0000          , 0xFF000000
                           };
const unsigned long lpat[]={
		                    0xFF            , 0xFFFF              , 0xFFFFFF          , 0xFFFFFFFF
                           };
const unsigned long rpat[]={
		                    0xFFFFFFFF      , 0xFFFFFF00          , 0xFFFF0000        , 0xFF000000
                           };

const unsigned bitmap_size=4;
const int max_bitmap_index=3;
// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#endif

typedef unsigned int       index_type;
typedef unsigned char      value_type;

const unsigned int nQsm=4; // Only nQsm largest values are included in Qs vectors

value_type getByte(const bitmap_type val,const index_type i) {return (value_type) (val>>idx[i])&0xFF;}

index_type bitmask_nPos(const bitmap_type val)
{
	index_type res=0;
	unsigned i;
	for(i=0;i<bitmap_size;++i)
	{
		// getByte
		if((val>>idx[i])&0xFF)
			++res;
	}
	return res;
}

index_type bitmask_sumPos(const bitmap_type val)
{
	index_type res=0;
	unsigned i;
	for(i=0;i<bitmap_size;++i)
		res+=(index_type)(val>>idx[i])&0xFF;
	return res;
}

// Calculates mean of four rightmost values
// (= largest values when r_type is used)
// This value is not monotone increasing on adding aligns
unsigned getQmean(const bitmap_type val)
{
	unsigned res,gb,n;
	res=(unsigned)(val>>idx[0])&0xFF;
	n=1;

	gb=(unsigned)(val>>idx[1])&0xFF;
	if(gb>0)
	{
		res+=gb;
		++n;
	}
	gb=(unsigned)(val>>idx[2])&0xFF;
	if(gb>0)
	{
		res+=gb;
		++n;
	}
	gb=(unsigned)(val>>idx[3])&0xFF;
	if(gb>0)
	{
		res+=gb;
		++n;
	}
	return res/n;
}

// Sums up nr of matching positions
// smcl = sum of minimal cigar length
unsigned getSmcl(const bitmap_type val)
{
	unsigned res=0,i;
	for(i=0;i<nQsm;++i)
		res+=(unsigned)(val>>idx[i])&0xFF;
	return res;
}

void r_insByte(bitmap_type *map,const value_type val, const index_type block)
{
	// left insert Byte: Creates l -> r ascending values
	// block: 0-based
	bitmap_type shift_value=(*map)&rpat[block];
	// reset shift_block
	(*map)&=~rpat[block];
	// reinsert byte-shifted value
	(*map)|=(shift_value<<8);
	// Conversion induces copy of val into larger sized variable
	(*map)|=(((bitmap_type)val)<<idx[block]);
}

void r_addVal(bitmap_type *map,const value_type val)
{
	// right add Value: l -> r ascending values
	// when val is equal to some found value, the function returns without insert
	// block values: 0-based
	index_type block=0;
	value_type byte_val;
	while(block<bitmap_size)
	{
		byte_val=(value_type) ((*map)>>idx[block])&0xFF; //getByte(*map,block);
		if(val==byte_val)
			return;
		if(val>byte_val)
		{
			r_insByte(map,val,block);
			return;
		}
		++block;
	}
}

void r_addVal_f(bitmap_type *map,const value_type val) // force version
{
	// right add Value: l -> r ascending values
	// the function also inserts, when val is equal to some found value
	// block values: 0-based
	index_type block=0;
	value_type byte_val;
	while(block<bitmap_size)
	{
		byte_val=(value_type) ((*map)>>idx[block])&0xFF; //getByte(*map,block);
		if(val>=byte_val)
		{
			r_insByte(map,val,block);
			return;
		}
		++block;
	}
}

void l_insByte(bitmap_type *map,const value_type val, const index_type block)
{
	// right insert Byte: Creates l -> r descending values
	// block values: 0-based
	bitmap_type shift_value,shift_block=0;
	shift_block=lpat[block];
	shift_value=*map&shift_block;
	// reset shift_block
	*map&=~shift_block;
	// reinsert byte-shifted value
	*map|=(shift_value>>8);
	// copy val into larger sized variable
	*map|=(((bitmap_type) val)<<idx[block]);
}

void l_addVal(bitmap_type *map,const value_type val)
{
	// right add Value: l -> r descending values
	// block values: 0-based
	int block=max_bitmap_index;
	value_type byte_val;
	while(block>=0)
	{
		byte_val=(value_type) ((*map)>>idx[block])&0xFF; //getByte(*map,block);
		if(val==byte_val)
			return;
		if(val>byte_val)
		{
			l_insByte(map,val,(index_type)block);
			return;
		}
		--block;
	}
}

void r_zip(bitmap_type lhs, bitmap_type rhs,bitmap_type *res)
{
	// merges two bitmap_types in order to get
	// an ordered array of greatest values

	// When lhs and rhs contain equal values,
	// only one of them is inserted into res.

	// l-type: l -> r ascending
	index_type lblock,rblock,res_block=1;
	value_type lByte,rByte;

	lByte=(value_type)(lhs>>idx[0])&0xFF; //getByte(lhs,0);
	rByte=(value_type)(rhs>>idx[0])&0xFF; //getByte(rhs,0);
	if(lByte>rByte)
	{
		*res=lByte;
		lblock=1;
		rblock=0;
	}
	else if(lByte<rByte)
	{
		*res=rByte;
		rblock=1;
		lblock=0;
	}
	else // lbyte==rbyte
	{
		*res=lByte;
		lblock=1;
		rblock=1;
	}

	while(res_block < bitmap_size)
	{
		//printf("[r_zip] lblock: %u\trblock: %u\tres_block: %u\n",lblock,rblock,res_block);
		//print_bitmask_dec(res);
		lByte=(value_type) (lhs>>idx[lblock])&0xFF; //getByte(lhs,lblock);
		rByte=(value_type) (rhs>>idx[rblock])&0xFF; //getByte(rhs,rblock);
		if(lByte>rByte)
		{
			//setByte(res,lByte,res_block);
			*res&=~pat[res_block];
			*res|=(((bitmap_type)lByte)<<idx[res_block]);
			++lblock;
		}
		else if(rByte>lByte)
		{
			//setByte(res,rByte,res_block);
			*res&=~pat[res_block];
			*res|=(((bitmap_type)rByte)<<idx[res_block]);
			++rblock;
		}
		else // rByte==lbyte: only one value is inserted
		{
			//setByte(res,lByte,res_block);
			*res&=~pat[res_block];
			*res|=(((bitmap_type)lByte)<<idx[res_block]);
			++lblock;
			++rblock;
		}
		++res_block;
	}
}

void r_zip_f(bitmap_type lhs, bitmap_type rhs,bitmap_type *res)  // force version
{
	// merges two bitmap_types in order to get
	// an ordered array of greatest values

	// When lhs and rhs contain equal values,
	// both of them are inserted into res.

	// l-type: l -> r ascending
	index_type lblock,rblock,res_block=1;
	value_type lByte,rByte;

	lByte=(value_type)(lhs>>idx[0])&0xFF; //getByte(lhs,0);
	rByte=(value_type)(rhs>>idx[0])&0xFF; //getByte(rhs,0);
	if(lByte>rByte)
	{
		*res=lByte;
		lblock=1;
		rblock=0;
	}
	else if(lByte<rByte)
	{
		*res=rByte;
		rblock=1;
		lblock=0;
	}
	else // lbyte==rbyte
	{
		*res=lByte;

		// insert second value
		*res&=~pat[res_block];
		*res|=(((bitmap_type)rByte)<<idx[res_block]);
		++res_block;

		lblock=1;
		rblock=1;
	}

	while(res_block < bitmap_size)
	{
		//printf("[r_zip] lblock: %u\trblock: %u\tres_block: %u\n",lblock,rblock,res_block);
		//print_bitmask_dec(res);
		lByte=(value_type) (lhs>>idx[lblock])&0xFF; //getByte(lhs,lblock);
		rByte=(value_type) (rhs>>idx[rblock])&0xFF; //getByte(rhs,rblock);
		if(lByte>rByte)
		{
			//setByte(res,lByte,res_block);
			*res&=~pat[res_block];
			*res|=(((bitmap_type)lByte)<<idx[res_block]);
			++lblock;
		}
		else if(rByte>lByte)
		{
			//setByte(res,rByte,res_block);
			*res&=~pat[res_block];
			*res|=(((bitmap_type)rByte)<<idx[res_block]);
			++rblock;
		}
		else // rByte==lbyte: both values are inserted
		{
			//setByte(res,lByte,res_block);
			*res&=~pat[res_block];
			*res|=(((bitmap_type)lByte)<<idx[res_block]);
			++lblock;

			++res_block;
			if(res_block < bitmap_size)
			{
				*res&=~pat[res_block];
				*res|=(((bitmap_type)rByte)<<idx[res_block]);
				++rblock;
			}
		}
		++res_block;
	}
}

bitmap_type r_zip_val(bitmap_type lhs, bitmap_type rhs)
{
	// merges two bitmap_types in order to get
	// an ordered array of greatest values

	// When lhs and rhs contain equal values,
	// only one of them is inserted into res.

	// l-type: l -> r ascending
	index_type lblock,rblock,res_block=1;
	value_type lByte,rByte;
	bitmap_type res;

	lByte=(value_type)(lhs>>idx[0])&0xFF; //getByte(lhs,0);
	rByte=(value_type)(rhs>>idx[0])&0xFF; //getByte(rhs,0);
	if(lByte>rByte)
	{
		res=lByte;
		lblock=1;
		rblock=0;
	}
	else if(lByte<rByte)
	{
		res=rByte;
		rblock=1;
		lblock=0;
	}
	else
	{
		res=lByte;
		lblock=1;
		rblock=1;
	}

	while(res_block < bitmap_size)
	{
		//printf("[r_zip] lblock: %u\trblock: %u\tres_block: %u\n",lblock,rblock,res_block);
		//print_bitmask_dec(res);

		lByte=(value_type)(lhs>>idx[lblock])&0xFF; //getByte(lhs,lblock);
		rByte=(value_type)(rhs>>idx[rblock])&0xFF; //getByte(rhs,rblock);
		if(lByte>rByte)
		{
			//setByte(&res,lByte,res_block);
			res&=~pat[res_block];
			res|=(((bitmap_type)lByte)<<idx[res_block]);
			++lblock;
		}
		else if(rByte>lByte)
		{
			//setByte(&res,rByte,res_block);
			res&=~pat[res_block];
			res|=(((bitmap_type)rByte)<<idx[res_block]);
			++rblock;
		}
		else // rByte==lByte
		{
			//setByte(&res,lByte,res_block);
			res&=~pat[res_block];
			res|=(((bitmap_type)lByte)<<idx[res_block]);
			++lblock;
			++rblock;
		}
		++res_block;
	}
	return res;
}

bitmap_type r_zip_val_f(bitmap_type lhs, bitmap_type rhs)
{
	// merges two bitmap_types in order to get
	// an ordered array of greatest values

	// When lhs and rhs contain equal values,
	// both of them are inserted into res.

	// l-type: l -> r ascending
	index_type lblock,rblock,res_block=1;
	value_type lByte,rByte;
	bitmap_type res;

	lByte=(value_type)(lhs>>idx[0])&0xFF; //getByte(lhs,0);
	rByte=(value_type)(rhs>>idx[0])&0xFF; //getByte(rhs,0);
	if(lByte>rByte)
	{
		res=lByte;

		// insert second value
		res&=~pat[res_block];
		res|=(((bitmap_type)rByte)<<idx[res_block]);
		++res_block;

		lblock=1;
		rblock=0;
	}
	else if(lByte<rByte)
	{
		res=rByte;
		rblock=1;
		lblock=0;
	}
	else
	{
		res=lByte;
		lblock=1;
		rblock=1;
	}

	while(res_block < bitmap_size)
	{
		//printf("[r_zip] lblock: %u\trblock: %u\tres_block: %u\n",lblock,rblock,res_block);
		//print_bitmask_dec(res);

		lByte=(value_type)(lhs>>idx[lblock])&0xFF; //getByte(lhs,lblock);
		rByte=(value_type)(rhs>>idx[rblock])&0xFF; //getByte(rhs,rblock);
		if(lByte>rByte)
		{
			//setByte(&res,lByte,res_block);
			res&=~pat[res_block];
			res|=(((bitmap_type)lByte)<<idx[res_block]);
			++lblock;
		}
		else if(rByte>lByte)
		{
			//setByte(&res,rByte,res_block);
			res&=~pat[res_block];
			res|=(((bitmap_type)rByte)<<idx[res_block]);
			++rblock;
		}
		else // rByte==lByte
		{
			//setByte(&res,lByte,res_block);
			res&=~pat[res_block];
			res|=(((bitmap_type)lByte)<<idx[res_block]);
			++lblock;

			++res_block;
			if(res_block < bitmap_size)
			{
				res&=~pat[res_block];
				res|=(((bitmap_type)rByte)<<idx[res_block]);
				++rblock;
			}
		}
		++res_block;
	}
	return res;
}


static R_INLINE void l_zip(bitmap_type lhs, bitmap_type rhs, bitmap_type *res)
{
	// merges two bitmap_types in order to get
	// an ordered array of greatest values
	// r-type: r -> l ascending
	int lblock,rblock,res_block=max_bitmap_index-1;
	value_type lByte,rByte;

	lByte=(value_type)(lhs>>idx[max_bitmap_index])&0xFF; //getByte(lhs,max_bitmap_index);
	rByte=(value_type)(rhs>>idx[max_bitmap_index])&0xFF; //getByte(rhs,max_bitmap_index);
	if(lByte>rByte)
	{
		l_insByte(res,lByte,bitmap_size-1);
		lblock=max_bitmap_index-1;
		rblock=(int)max_bitmap_index;
	}
	else if(lByte<rByte)
	{
		l_insByte(res,rByte,bitmap_size-1);
		rblock=max_bitmap_index-1;
		lblock=(int)max_bitmap_index;
	}
	else
	{
		l_insByte(res,lByte,bitmap_size-1);
		lblock=max_bitmap_index-1;
		rblock=lblock;
	}
	while(res_block >=0)
	{
		//printf("[l_zip] lblock: %u\trblock: %u\tres_block: %u\n",lblock,rblock,res_block);
		//print_bitmask_dec(res);

		lByte=(value_type)(lhs>>idx[lblock])&0xFF; //getByte(lhs,lblock);
		rByte=(value_type)(rhs>>idx[rblock])&0xFF; //getByte(rhs,rblock);
		if(lByte>rByte)
		{
			//setByte(res,lByte,res_block);
			*res&=~pat[res_block];
			*res|=(((bitmap_type)lByte)<<idx[res_block]);
			--lblock;
		}
		else if(rByte>lByte)
		{
			//setByte(res,rByte,res_block);
			*res&=~pat[res_block];
			*res|=(((bitmap_type)rByte)<<idx[res_block]);
			--rblock;
		}
		else
		{
			//setByte(res,lByte,res_block);
			*res&=~pat[res_block];
			*res|=(((bitmap_type)lByte)<<idx[res_block]);
			--lblock;
			--rblock;
		}
		--res_block;
	}
}

#endif /* BITMASK_H_ */
