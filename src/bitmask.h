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
 */

#ifndef BITMASK_H_
#define BITMASK_H_

// This checks whether this is compiled in 64 bit
#include <stdint.h>
#if UINTPTR_MAX == 0xffffffffffffffff
#define BM_64
#endif


#ifdef BM_64
// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// This results in storing values in 64-bit (unsigned) integers:
// 8 byte-sized fields are used
// (Also used by gapSiteList)
typedef unsigned long long  bitmap_type; 							// Size at least 64 bits (C99)
const unsigned idx[]      ={0,8,16,24,32,40,48,56};
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
const unsigned max_bitmap_index=7;
// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#else
// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// This results in storing values in (at least) 32-bit (unsigned) integers
// 4 byte-sized fields are used
// (Also used by gapSiteList)
typedef unsigned long       bitmap_type;							// Size at least 32 bits (C99)
const unsigned idx[]      ={0,8,16,24};

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
const unsigned max_bitmap_index=3;
// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#endif

typedef unsigned int       index_type;
typedef unsigned char      value_type;


inline value_type getByte(const bitmap_type val,const index_type i) {return (value_type) (val>>idx[i])&0xFF;}

inline index_type bitmask_nPos(const bitmap_type val)
{
	index_type res=0;
	unsigned i;
	for(i=0;i<bitmap_size;++i)
	{
		if(getByte(val,i))
			++res;
	}
	return res;
}

inline index_type bitmask_sumPos(const bitmap_type val)
{
	index_type res=0;
	unsigned i;
	for(i=0;i<bitmap_size;++i)
		res+=getByte(val,i);
	return res;
}

// Calculates mean of four rightmost values
// (= largest values when r_type is used)
inline unsigned getQmean(const bitmap_type val)
{
	unsigned res,gb,n;
	res =(unsigned)getByte(val,0);
	n=1;

	gb=(unsigned)getByte(val,1);
	if(gb>0)
	{
		res+=gb;
		++n;
	}
	gb=(unsigned)getByte(val,2);
	if(gb>0)
	{
		res+=gb;
		++n;
	}
	gb=(unsigned)getByte(val,3);
	if(gb>0)
	{
		res+=gb;
		++n;
	}
	return res/n;
}


inline void setByte(bitmap_type *val,const value_type ins,const index_type i)
{
	// reset byte
	*val&= ~pat[i];
	// set new value
	// does not work without type conversion
	// = copy into larger sized variable
	*val|= (((bitmap_type)ins)<<idx[i]);
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
	// block values: 0-based
	index_type block=0;
	value_type byte_val=getByte(*map,block);
	while(block<bitmap_size)
	{
		if(val==byte_val)
			return;
		if(val>byte_val)
		{
			r_insByte(map,val,block);
			return;
		}
		++block;
		byte_val=getByte(*map,block);
	}
	//setByte(map,val,bitmap_size-1);
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
		byte_val=getByte(*map,block);
		if(val==byte_val)
			return;
		if(val>byte_val)
		{
			l_insByte(map,val,block);
			return;
		}
		--block;
	}
}

void r_zip(bitmap_type lhs, bitmap_type rhs,bitmap_type *res)
{
	// merges two bitmap_types in order to get
	// an ordered array of greatest values
	// l-type: l -> r ascending
	index_type lblock,rblock,res_block=1;
	value_type lByte,rByte;

	lByte=getByte(lhs,0);
	rByte=getByte(rhs,0);
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
	else
	{
		*res=lByte;
		lblock=1;
		rblock=1;
	}

	while(res_block < bitmap_size)
	{
		//printf("[r_zip] lblock: %u\trblock: %u\tres_block: %u\n",lblock,rblock,res_block);
		//print_bitmask_dec(res);
		lByte=getByte(lhs,lblock);
		rByte=getByte(rhs,rblock);
		if(lByte>rByte)
		{
			setByte(res,lByte,res_block);
			++lblock;
		}
		else if(rByte>lByte)
		{
			setByte(res,rByte,res_block);
			++rblock;
		}
		else
		{
			setByte(res,lByte,res_block);
			++lblock;
			++rblock;
		}
		++res_block;
	}
}


bitmap_type r_zip_val(bitmap_type lhs, bitmap_type rhs)
{
	// merges two bitmap_types in order to get
	// an ordered array of greatest values
	// l-type: l -> r ascending
	index_type lblock,rblock,res_block=1;
	value_type lByte,rByte;
	bitmap_type res;

	lByte=getByte(lhs,0);
	rByte=getByte(rhs,0);
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
		lByte=getByte(lhs,lblock);
		rByte=getByte(rhs,rblock);
		if(lByte>rByte)
		{
			setByte(&res,lByte,res_block);
			++lblock;
		}
		else if(rByte>lByte)
		{
			setByte(&res,rByte,res_block);
			++rblock;
		}
		else
		{
			setByte(&res,lByte,res_block);
			++lblock;
			++rblock;
		}
		++res_block;
	}
	return res;
}


inline void l_zip(bitmap_type lhs, bitmap_type rhs, bitmap_type *res)
{
	// merges two bitmap_types in order to get
	// an ordered array of greatest values
	// r-type: r -> l ascending
	int lblock,rblock,res_block=max_bitmap_index-1;
	value_type lByte,rByte;

	lByte=getByte(lhs,max_bitmap_index);
	rByte=getByte(rhs,max_bitmap_index);
	if(lByte>rByte)
	{
		l_insByte(res,lByte,bitmap_size-1);
		lblock=max_bitmap_index-1;
		rblock=max_bitmap_index;
	}
	else if(lByte<rByte)
	{
		l_insByte(res,rByte,bitmap_size-1);
		rblock=max_bitmap_index-1;
		lblock=max_bitmap_index;
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

		lByte=getByte(lhs,lblock);
		rByte=getByte(rhs,rblock);
		if(lByte>rByte)
		{
			setByte(res,lByte,res_block);
			--lblock;
		}
		else if(rByte>lByte)
		{
			setByte(res,rByte,res_block);
			--rblock;
		}
		else
		{
			setByte(res,lByte,res_block);
			--lblock;
			--rblock;
		}
		--res_block;
	}
}

#endif /* BITMASK_H_ */
