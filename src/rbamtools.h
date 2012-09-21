/*
 *	File:		rbamtools.c
 *
 * 	Created on:	17.06.2011
 *  Author: 	Wolfgang Kaisers
 *	Content:	C Header File for R package rbamtools
 */

#ifndef rbamtools_h
#define rbamtools_h

#include <string.h>

#include <R.h>
#include <Rinternals.h>
#include <Rdefines.h>
#include <R_ext/PrtUtil.h>
#include "samtools/bam.h"
#include "samtools/sam.h"
#include "align_list.h"
#include "gap_list.h"

const char * const CIGAR_TYPES="MIDNSHP=X";

inline int cigar2str(char *c,const bam1_t *align);
bam_header_t* clone_bam_header(bam_header_t *h);
SEXP is_nil_externalptr(SEXP ptr);

///////////////////////////////////////////////////////////////////////////////////////////////////
// bamHeader
///////////////////////////////////////////////////////////////////////////////////////////////////

static void finalize_bam_header(SEXP ptr);
SEXP init_bam_header(SEXP pHeaderText);
SEXP bam_header_get_header_text(SEXP pHeader);

///////////////////////////////////////////////////////////////////////////////////////////////////
// BamWriter
///////////////////////////////////////////////////////////////////////////////////////////////////
static void finalize_bam_writer(SEXP ptr);
SEXP bam_writer_open(SEXP pHeader,SEXP pFilename);
SEXP bam_reader_open_writer(SEXP pReader,SEXP pFilename);
SEXP bam_writer_save_align(SEXP pWriter, SEXP pAlign);
SEXP bam_writer_close(SEXP pWriter);

///////////////////////////////////////////////////////////////////////////////////////////////////
// BamReader
///////////////////////////////////////////////////////////////////////////////////////////////////
static void finalize_bam_reader(SEXP ptr);
static void finalize_bam_index(SEXP ptr);
SEXP bam_reader_open(SEXP filename);
SEXP bam_reader_close(SEXP pReader);
SEXP bam_reader_get_header_text(SEXP pReader);
SEXP bam_reader_get_ref_count(SEXP pReader);
SEXP bam_reader_get_ref_data(SEXP pReader);
SEXP bam_reader_create_index(SEXP bam_file,SEXP idx_file);
SEXP bam_reader_load_index(SEXP idx_file);
SEXP bam_reader_unload_index(SEXP pIdx);
SEXP bam_reader_get_next_align(SEXP pReader);
SEXP bam_reader_sort_file(SEXP pFilename,SEXP pPrefix,SEXP pMaxMem,SEXP pByName);
SEXP bam_reader_get_header(SEXP pReader);

///////////////////////////////////////////////////////////////////////////////////////////////////
// GapList
///////////////////////////////////////////////////////////////////////////////////////////////////

static void finalize_gap_list(SEXP ptr);
SEXP create_gap_list();
static int gap_fetch_func(const bam1_t *b, void *data);
SEXP gap_list_fetch(SEXP pReader,SEXP pIndex,SEXP pCoords);
SEXP gap_list_get_size(SEXP pGapList);
SEXP get_gap_list_df(SEXP pGapList);

///////////////////////////////////////////////////////////////////////////////////////////////////
// BamRange
///////////////////////////////////////////////////////////////////////////////////////////////////

static void finalize_bam_range(SEXP ptr);
static int range_fetch_func(const bam1_t *b, void *data);
static int range_fetch_complex_func(const bam1_t *b,void *data);
SEXP bam_range_init();
SEXP bam_range_fetch(SEXP pReader,SEXP pIndex,SEXP pCoords,SEXP pComplex);
SEXP bam_range_get_size(SEXP pRange);
SEXP bam_range_get_next_align(SEXP pRange);
SEXP bam_range_get_prev_align(SEXP pRange);
SEXP bam_range_step_next_align(SEXP pRange);
SEXP bam_range_step_prev_align(SEXP pRange);
SEXP bam_range_get_align_df(SEXP pRange);
SEXP bam_range_write(SEXP pWriter,SEXP pRange);
SEXP bam_range_wind_back(SEXP pRange);
SEXP bam_range_push_back(SEXP pRange,SEXP pAlign);
SEXP bam_range_pop_back(SEXP pRange);
SEXP bam_range_push_front(SEXP pRange,SEXP pAlign);
SEXP bam_range_pop_front(SEXP pRange);
SEXP bam_range_write_current_align(SEXP pRange,SEXP pAlign);
SEXP bam_range_insert_past_curr_align(SEXP pRange,SEXP pAlign);
SEXP bam_range_insert_pre_curr_align(SEXP pRange,SEXP pAlign);
SEXP bam_range_mv_curr_align(SEXP pSrc, SEXP pTarget);

///////////////////////////////////////////////////////////////////////////////////////////////////
// BamAlignment
///////////////////////////////////////////////////////////////////////////////////////////////////
static void finalize_bam_align(SEXP pAlign);
SEXP bam_align_get_name(SEXP pAlign);
SEXP bam_align_get_refid(SEXP pAlign);
SEXP bam_align_get_position(SEXP pAlign);
SEXP bam_align_get_nCigar(SEXP pAlign);
SEXP bam_align_get_cigar_df(SEXP pAlign);
SEXP bam_align_get_mate_refid(SEXP pAlign);
SEXP bam_align_get_mate_position(SEXP pAlign);
SEXP bam_align_get_insert_size(SEXP pAlign);
SEXP bam_align_get_map_quality(SEXP pAlign);
SEXP bam_align_get_segment_sequence(SEXP pAlign);
SEXP bam_align_get_qualities(SEXP pAlign);

///////////////////////////////////////////////////////////
// alignment flags

// Reading accessors
SEXP bam_align_is_paired(SEXP pAlign);//
SEXP bam_align_mapped_in_proper_pair(SEXP pAlign);//
SEXP bam_align_is_unmapped(SEXP pAlign);//
SEXP bam_align_mate_is_unmapped(SEXP pAlign);//
SEXP bam_align_strand_reverse(SEXP pAlign);//
SEXP bam_align_mate_strand_reverse(SEXP pAlign);//
SEXP bam_align_is_first_in_pair(SEXP pAlign);//
SEXP bam_align_is_second_in_pair(SEXP pAlign);//
SEXP bam_align_is_secondary_align(SEXP pAlign);
SEXP bam_align_fail_qc(SEXP pAlign);//
SEXP bam_align_is_pcr_or_optical_dup(SEXP pAlign);//
SEXP bam_align_get_flag(SEXP pAlign);

// Writing accessors
SEXP bam_align_set_is_paired(SEXP pAlign, SEXP val);
SEXP bam_align_set_mapped_in_proper_pair(SEXP pAlign, SEXP val);
SEXP bam_align_set_is_unmapped(SEXP pAlign, SEXP val);
SEXP bam_align_set_mate_is_unmapped(SEXP pAlign, SEXP val);
SEXP bam_align_set_strand_reverse(SEXP pAlign, SEXP val);
SEXP bam_align_set_mate_strand_reverse(SEXP pAlign, SEXP val);
SEXP bam_align_set_is_first_in_pair(SEXP pAlign, SEXP val);
SEXP bam_align_set_is_second_in_pair(SEXP pAlign, SEXP val);
SEXP bam_align_set_is_secondary_align(SEXP pAlign, SEXP val);
SEXP bam_align_set_fail_qc(SEXP pAlign, SEXP val);
SEXP bam_align_set_is_pcr_or_optical_dup(SEXP pAlign, SEXP val);
SEXP bam_align_set_flag(SEXP pAlign, SEXP val);

#endif
