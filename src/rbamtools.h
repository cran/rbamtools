/*
 *	File:		rbamtools.c
 *
 * 	Created on:	17.06.2011
 * 	Author: 	Wolfgang Kaisers
 *	Content:	C Header File for R package rbamtools
 *
 *	Change log:
 *      29.Okt.12 (nemesis) [gap_list_get_df] Changed cigar_type output to 'factor'.
 *      30.Okt.12 (phoibe)  [gap_list_fetch] Removed gap_list_fetch message.
 *      31.Okt.12 (phoibe)  [bam_reader_save_aligns] function & SAM_TYPE_READ added.
 *		01.Nov.12 (nemesis) [get_const_next_align] added to correct memory leak.
 *		07.Jan.13 (phoibe)  gap_site_list added (contains bitmask functions)
 *		18.Jan.13 (phoibe)  gap_site_list_list added.
 *		01.Feb.13 (phoibe)  [gap_site_ll_add_curr_pp] Testing for empty input list added.
 *		05.Feb.13 (phoibe)  First successful test with reading bamGap objects with transfer to gapProbes
 *		18.Mar.13 (nemesis) Corrected inline declarations; Changed qmm to qsm (mean to sum)
 * 		11.Jun.13 (phoibe)  Added bam_reader_write_fastq, bam_reader_write_fastq_lgl,
 * 									bam_range_write_fastq, bam_range_write_fastq_lgl functions. Valgrind tested.
 *      02.Jul.13 (phoibe)  Added bam_count. Valgrind tested
 *      02.Sep.13 (phoibe)	Added R_init_rbamtools
 *      25.Nov.13 (gaia)
 */

#ifndef rbamtools_h
#define rbamtools_h

#include <string.h>
#include <ctype.h>

#include <R.h>
#include <Rinternals.h>
#include <Rdefines.h>
#include <R_ext/PrtUtil.h>
#include <R_ext/Rdynload.h> // DllInfo
#include "samtools/bam.h"
#include "samtools/sam.h"
#include "align_list.h"
#include "gap_list.h"

// Extended gap_list
#include "gapSiteList.h"
#include "bitmask.h"
#include "gapSiteListList.h"

const char * const CIGAR_TYPES="MIDNSHP=X";
#define SAM_TYPE_READ 2

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
SEXP bam_writer_save_align(SEXP pWriter,SEXP pAlign,SEXP pRefid);
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
SEXP bam_reader_get_target_name(SEXP pReader,SEXP pSeqid);
SEXP bam_reader_get_ref_data(SEXP pReader);
SEXP bam_reader_create_index(SEXP bam_file,SEXP idx_file);
SEXP bam_reader_load_index(SEXP idx_file);
SEXP bam_reader_unload_index(SEXP pIdx);
SEXP bam_reader_get_next_align(SEXP pReader);
SEXP bam_reader_save_aligns(SEXP pReader,SEXP pWriter);
SEXP bam_reader_sort_file(SEXP pFilename,SEXP pPrefix,SEXP pMaxMem,SEXP pByName);
SEXP bam_reader_get_header(SEXP pReader);
SEXP bam_reader_tell(SEXP pReader);
SEXP bam_reader_seek(SEXP pReader, SEXP pPos);
SEXP bam_reader_write_fastq(SEXP pReader,SEXP pFilename,SEXP pAppend);
SEXP bam_reader_write_fastq_index(SEXP pReader,SEXP pFilename,SEXP pWhichWrite,SEXP pAppend);

static int bam_count_fetch_func(const bam1_t *align, void *data);
SEXP bam_count(SEXP pReader,SEXP pIndex,SEXP pCoords);

///////////////////////////////////////////////////////////////////////////////////////////////////
// gap_list
///////////////////////////////////////////////////////////////////////////////////////////////////

static void finalize_gap_list(SEXP ptr);
SEXP create_gap_list();
static int gap_fetch_func(const bam1_t *b, void *data);
SEXP gap_list_fetch(SEXP pReader,SEXP pIndex,SEXP pCoords);
SEXP gap_list_get_df(SEXP pGapList);
SEXP gap_list_get_size(SEXP pGapList);
SEXP gap_list_get_nAligns(SEXP pGapList);
SEXP gap_list_get_nAlignGaps(SEXP pGapList);

///////////////////////////////////////////////////////////////////////////////////////////////////
// gap_site_list
///////////////////////////////////////////////////////////////////////////////////////////////////

static void finalize_gap_site_list(SEXP ptr);
SEXP create_gap_site_list();
static int gap_site_list_fetch_func(const bam1_t *b, void *data);
SEXP gap_site_list_fetch(SEXP pReader,SEXP pIndex,SEXP pCoords);
SEXP gap_site_list_get_df(SEXP pGapList);
SEXP gap_site_list_get_ref_id(SEXP pGapList);
SEXP gap_site_list_get_size(SEXP pGapList);
SEXP gap_site_list_get_nAligns(SEXP pGapList);
SEXP gap_site_list_get_nAlignGaps(SEXP pGapList);
SEXP gap_site_list_merge(SEXP pLhs, SEXP pRhs, SEXP pRef);
SEXP gap_site_list_copy(SEXP pGapList);
SEXP bitmask_r_zip(SEXP lhs, SEXP rhs);

///////////////////////////////////////////////////////////////////////////////////////////////////
// gap_site_list_list
///////////////////////////////////////////////////////////////////////////////////////////////////

static void finialize_gap_site_ll(SEXP ptr);
SEXP gap_site_ll_init();
SEXP gap_site_ll_fetch(SEXP pReader, SEXP pIndex, SEXP pRefid, SEXP pStart, SEXP pEnd);
SEXP gap_site_ll_get_df(SEXP pGapList,SEXP pRefNames);
SEXP gap_site_ll_get_size(SEXP pGapList);
SEXP gap_site_ll_get_nAligns(SEXP pGapList);
SEXP gap_site_ll_get_nAlignGaps(SEXP pGapList);
SEXP gap_site_ll_add_curr_pp(SEXP pSrc,SEXP pTrg,SEXP pRefid);
SEXP gap_site_ll_add_merge_pp(SEXP plSrc,SEXP prSrc,SEXP pTrg,SEXP pRefid);
SEXP gap_site_ll_reset_refid(SEXP pGapList);
SEXP gap_site_ll_get_summary_df(SEXP pGapList);
SEXP gap_site_ll_set_curr_first(SEXP pGapList);

///////////////////////////////////////////////////////////////////////////////////////////////////
// bam_range
///////////////////////////////////////////////////////////////////////////////////////////////////

static void finalize_bam_range(SEXP ptr);
static int range_fetch_func(const bam1_t *b, void *data);
static int range_fetch_complex_func(const bam1_t *b,void *data);
SEXP bam_range_init();
SEXP bam_range_fetch(SEXP pReader,SEXP pIndex,SEXP pCoords,SEXP pComplex);
SEXP bam_range_get_size(SEXP pRange);
SEXP bam_range_get_coords(SEXP pRange);
SEXP bam_range_get_params(SEXP pRange);
SEXP bam_range_get_refname(SEXP pRange);
SEXP bam_range_get_align_range(SEXP pRange);
SEXP bam_range_get_next_align(SEXP pRange);
SEXP bam_range_get_prev_align(SEXP pRange);
SEXP bam_range_step_next_align(SEXP pRange);
SEXP bam_range_step_prev_align(SEXP pRange);
SEXP bam_range_get_align_df(SEXP pRange);
SEXP bam_range_write(SEXP pWriter,SEXP pRange,SEXP pRefid);
SEXP bam_range_wind_back(SEXP pRange);
SEXP bam_range_push_back(SEXP pRange,SEXP pAlign);
SEXP bam_range_pop_back(SEXP pRange);
SEXP bam_range_push_front(SEXP pRange,SEXP pAlign);
SEXP bam_range_pop_front(SEXP pRange);
SEXP bam_range_write_current_align(SEXP pRange,SEXP pAlign);
SEXP bam_range_insert_past_curr_align(SEXP pRange,SEXP pAlign);
SEXP bam_range_insert_pre_curr_align(SEXP pRange,SEXP pAlign);
SEXP bam_range_mv_curr_align(SEXP pSrc, SEXP pTarget);
SEXP bam_range_write_fastq(SEXP pRange,SEXP pFilename,SEXP pAppend);
SEXP bam_range_write_fastq_index(SEXP pRange,SEXP pFilename,SEXP pWhichWrite,SEXP pAppend);
SEXP bam_range_get_seqlen(SEXP pRange);
SEXP bam_range_get_qual_df(SEXP pRange);
SEXP bam_range_get_align_depth(SEXP pRange,SEXP pGap);
SEXP bam_range_count_nucs(SEXP pRange);
SEXP bam_range_idx_copy(SEXP pRange, SEXP pIndex);



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
SEXP bam_align_get_qual_values(SEXP pAlign);
SEXP bam_align_count_nucs(SEXP pAlign);

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
SEXP bam_align_set_refid(SEXP pAlign,SEXP pRefid);
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

// Create new bamAlign structure from scratch
SEXP bam_align_create(SEXP pStrVals, SEXP pIntVals);

///////////////////////////////////////////////////////////////////////////////////////////////////
// Miscellaneous functions
///////////////////////////////////////////////////////////////////////////////////////////////////
SEXP copy_fastq_records(SEXP pInfile,SEXP pOutfile,SEXP pWhichCopy,SEXP pAppend);
SEXP count_fastq(SEXP pInfile,SEXP pMaxCol);
SEXP get_col_quantiles(SEXP pQuant, SEXP pDf);
SEXP count_text_lines(SEXP pInfile);

///////////////////////////////////////////////////////////////////////////////////////////////////
// Declarations for R_registerRoutines
///////////////////////////////////////////////////////////////////////////////////////////////////

void R_init_rbamtools(DllInfo *info);

#endif
