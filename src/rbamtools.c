/*
 *	File:		rbamtools.c
 *
 *  Created on:	17.06.2011
 *  Author: 	W. Kaisers
 *	Content:	Function definitions in C for R package rbamtools
 */

#ifndef rbamtools_c
#define rbamtools_c
#include "rbamtools.h"

inline void clear_buf(char *c,unsigned n)
{
	int i;
	for(i=0;i<n;++i)
		c[i]=(char)0;
}

inline void set_flag(bam1_t *align,_Bool val,unsigned pattern)
{
	if(val)
		align->core.flag=(align->core.flag) | pattern;
	else
		align->core.flag=(align->core.flag) & !pattern;
}

inline int cigar2str(char *c,const bam1_t *align)
{
	if(align==NULL)
		return 0;

	uint32_t len=align->core.n_cigar;
	uint32_t *cigar=bam1_cigar(align);
	char buf[128];

	sprintf(buf,"%lu",(unsigned long) (cigar[0] >> BAM_CIGAR_SHIFT));
	strcpy(c,buf);
	if((cigar[0]&BAM_CIGAR_MASK)>=strlen(CIGAR_TYPES))	// Error
		return 0;
	strncat(c,&(CIGAR_TYPES[cigar[0] & BAM_CIGAR_MASK]),1);


	uint32_t i;
	for(i=1;i<len;++i)
	{
		sprintf(buf,"%lu",(unsigned long) (cigar[i] >> BAM_CIGAR_SHIFT));
		strncat(c,buf,strlen(buf));

		if((cigar[i]&BAM_CIGAR_MASK)>=strlen(CIGAR_TYPES))	// Error
			return 0;

		strncat(c,&(CIGAR_TYPES[cigar[i] & BAM_CIGAR_MASK]),1);
	}
	return strlen(c);
}

bam_header_t* clone_bam_header(bam_header_t *h)
{
	bam_header_t *ans=(bam_header_t*)calloc(1, sizeof(bam_header_t));
	ans->n_targets=h->n_targets;
	ans->l_text=h->l_text;
	ans->n_text=h->n_text;

	ans->text=(char*) calloc(1,(h->l_text)+1);
	strncpy(ans->text,h->text,h->l_text);
	sam_header_parse(ans);
	bam_init_header_hash(ans);
	return ans;
}


SEXP is_nil_externalptr(SEXP ptr)
{
	if(TYPEOF(ptr)!=EXTPTRSXP)
	{
		error("[is_nil_externalptr] No external pointer");
		return R_NilValue;
	}
	SEXP ans;
	PROTECT(ans=allocVector(LGLSXP,1));
	LOGICAL(ans)[0]=(R_ExternalPtrAddr(ptr)==NULL);
	UNPROTECT(1);
	return ans;
}

///////////////////////////////////////////////////////////////////////////////////////////////
// bam_reader
///////////////////////////////////////////////////////////////////////////////////////////////
static void finalize_bam_reader(SEXP ptr)
{
	if(TYPEOF(ptr)!=EXTPTRSXP)
	{
		error("[finalize_bam_reader] No external pointer!");
		return;
	}
	samfile_t *reader=(samfile_t*)(R_ExternalPtrAddr(ptr));
	samclose(reader);	// checks for 0
	R_SetExternalPtrAddr(ptr,NULL);
}

SEXP bam_reader_open(SEXP filename)
{
	if(TYPEOF(filename)!=STRSXP)
	{
		error("[bam_reader_open] Filename must be a string.\n");
		return R_NilValue;
	}
	const char* _filename=CHAR(STRING_ELT(filename,0));
	samfile_t *reader=samopen(_filename,"rb",0);
	if(!reader)
		error("[bam_reader_open] Opening bam_file \"%s\" failed!",_filename);
	else
		Rprintf("[bamReader] Opened file \"%s\"\n",_filename);

	SEXP pReader;
	PROTECT(pReader=R_MakeExternalPtr( (void*)(reader),R_NilValue,R_NilValue));
	R_RegisterCFinalizer(pReader,finalize_bam_reader);
	UNPROTECT(1);
	return pReader;
}

SEXP bam_reader_close(SEXP pReader)
{
	if(TYPEOF(pReader)!=EXTPTRSXP)
	{
		error("[bam_reader_close] No external pointer!");
		return R_NilValue;
	}
	samfile_t *reader=(samfile_t*) (R_ExternalPtrAddr(pReader));
	samclose(reader);
	R_SetExternalPtrAddr(pReader,NULL);
	Rprintf("[BamReader] closed.\n");
	return R_NilValue;
}


SEXP bam_reader_get_header_text(SEXP pReader)
{
	if(TYPEOF(pReader)!=EXTPTRSXP)
	{
		error("[bam_reader_get_header_text] No external pointer!");
		return R_NilValue;
	}
	samfile_t *reader=(samfile_t*) (R_ExternalPtrAddr(pReader));
	SEXP ans;
	PROTECT(ans=Rf_allocVector(STRSXP,1));
	SET_STRING_ELT(ans,0,mkChar(reader->header->text));
	UNPROTECT(1);
	return ans;
}

SEXP bam_reader_get_ref_count(SEXP pReader)
{
	if(TYPEOF(pReader)!=EXTPTRSXP)
	{
		error("[bam_reader_get_ref_count] No external pointer!");
		return R_NilValue;
	}
	samfile_t *reader=(samfile_t*) (R_ExternalPtrAddr(pReader));
	SEXP ans;
	PROTECT(ans=Rf_allocVector(INTSXP,1));
	INTEGER(ans)[0]=reader->header->n_targets;
	UNPROTECT(1);
	return ans;
}

SEXP bam_reader_get_ref_data(SEXP pReader)
{
	if(TYPEOF(pReader)!=EXTPTRSXP)
	{
		error("[bam_reader_get_ref_data] No external pointer!");
		return R_NilValue;
	}
	samfile_t *reader=(samfile_t*)(R_ExternalPtrAddr(pReader));
	bam_header_t* header=reader->header;

	// create data.frame
	int nProtected=0;
	int nCols=3;
	SEXP dflist;
	PROTECT(dflist=allocVector(VECSXP,nCols));
	++nProtected;
	int nRows=header->n_targets;

	// Column 0: ID (RefID)
	SEXP RefID_vector;
	PROTECT(RefID_vector=allocVector(INTSXP,nRows));
	++nProtected;

	// Column 1: SN (RefName)
	SEXP RefName_vector;
	PROTECT(RefName_vector=allocVector(STRSXP,nRows));
	++nProtected;

	// Column 2: LN (RefLength)
	SEXP RefLength_vector;
	PROTECT(RefLength_vector=allocVector(INTSXP,nRows));
	++nProtected;

	int i;
	for(i=0;i<nRows;++i)
	{
		INTEGER(RefID_vector)[i]=i;
		SET_STRING_ELT(RefName_vector,i,mkChar(header->target_name[i]));
		INTEGER(RefLength_vector)[i]=header->target_len[i];
	}
	SET_VECTOR_ELT(dflist,0,RefID_vector);
	SET_VECTOR_ELT(dflist,1,RefName_vector);
	SET_VECTOR_ELT(dflist,2,RefLength_vector);

	// Column Names
	SEXP col_names;
	PROTECT(col_names=allocVector(STRSXP,nCols));
	++nProtected;
	SET_STRING_ELT(col_names,0,mkChar("ID"));
	SET_STRING_ELT(col_names,1,mkChar("SN"));
	SET_STRING_ELT(col_names,2,mkChar("LN"));
	setAttrib(dflist,R_NamesSymbol,col_names);

	SEXP row_names;
    PROTECT(row_names=allocVector(STRSXP,nRows));
    ++nProtected;
    char c[20];
    for(i=1;i<=nRows;++i)
    {
    	sprintf(c,"%i",i);
    	SET_STRING_ELT(row_names,i-1,mkChar(c));
    }
    setAttrib(dflist,R_RowNamesSymbol,row_names);
	setAttrib(dflist,R_ClassSymbol,mkString("data.frame"));
	UNPROTECT(nProtected);
	return dflist;
}

SEXP bam_reader_create_index(SEXP pBamFile,SEXP pIdxFile)
{
	if(TYPEOF(pBamFile)!=STRSXP)
	{
		error("[bam_reader_create_index] BamFile must be a string!\n");
		return R_NilValue;
	}
	if(TYPEOF(pIdxFile)!=STRSXP)
	{
		error("[bam_reader_create_index] IndexFile must be a string!\n");
		return R_NilValue;
	}
	const char *bamFile=CHAR(STRING_ELT(pBamFile,0));
	const char *idxFile=CHAR(STRING_ELT(pIdxFile,0));

	SEXP ans;
	PROTECT(ans=Rf_allocVector(INTSXP,1));
	INTEGER(ans)[0]=bam_index_build2(bamFile,idxFile);
	UNPROTECT(1);
	return ans;
}

static void finalize_bam_index(SEXP ptr)
{
	if(TYPEOF(ptr)!=EXTPTRSXP)
	{
		error("[finalize_bam_index] No external pointer!");
		return;
	}
	bam_index_t *index=(bam_index_t *)(R_ExternalPtrAddr(ptr));
	bam_index_destroy(index);	// checks for zero
	R_SetExternalPtrAddr(ptr,NULL);
}

SEXP bam_reader_load_index(SEXP pIdxFile)
{
	if(TYPEOF(pIdxFile)!=STRSXP)
	{
		error("[bam_reader_load_index] pIdxFile must be a string!\n");
		return R_NilValue;
	}
	const char *idxFile=CHAR(STRING_ELT(pIdxFile,0));
	FILE *f=fopen(idxFile,"rb");

	bam_index_t *index = bam_index_load_core(f);
	SEXP idx;
	PROTECT(idx=R_MakeExternalPtr( (void*)(index),R_NilValue,R_NilValue));
	R_RegisterCFinalizer(idx,finalize_bam_index);
	UNPROTECT(1);
	return idx;
}

SEXP bam_reader_unload_index(SEXP pIdx)
{
	if(TYPEOF(pIdx)!=EXTPTRSXP)
	{
		error("[bam_reader_unload_index] No external pointer!\n");
		return R_NilValue;
	}
	bam_index_t *idx=(bam_index_t *)(R_ExternalPtrAddr(pIdx));
	bam_index_destroy(idx);
	R_SetExternalPtrAddr(pIdx,NULL);
	return R_NilValue;
}

SEXP bam_reader_get_next_align(SEXP pReader)
{
	if(TYPEOF(pReader)!=EXTPTRSXP)
	{
		error("[bam_reader_get_next_align] No external pointer!\n");
		return R_NilValue;
	}

	samfile_t *reader=(samfile_t*)(R_ExternalPtrAddr(pReader));
	bam1_t *align=bam_init1();
	int res=samread(reader,align);
	if(res==-1)
	{
		Rprintf("[getNextAlign] samread found EOF.\n");
		return R_NilValue;
	}
	if(res==-2)
	{
		error("[getNextAlign] samread found truncated BAM-file.\n");
		return R_NilValue;
	}

	SEXP ptr;
	PROTECT(ptr=R_MakeExternalPtr((void*)(align),R_NilValue,R_NilValue));
	R_RegisterCFinalizer(ptr,finalize_bam_align);
	UNPROTECT(1);
	return ptr;
}


SEXP bam_reader_sort_file(SEXP pFilename,SEXP pPrefix,SEXP pMaxMem,SEXP pByName)
{
	if(TYPEOF(pFilename)!=STRSXP)
	{
		error("[bam_writer_sort_file] Filename must be a string\n");
		return R_NilValue;
	}
	if(TYPEOF(pPrefix)!=STRSXP)
	{
		error("[bam_writer_sort_file] Prefix must be a string\n");
		return R_NilValue;
	}
	if(TYPEOF(pMaxMem)!=REALSXP)
	{
		error("[bam_writer_sort_file] MaxMem must be integer value!\n");
		return R_NilValue;
	}
	if(TYPEOF(pByName)!=LGLSXP)
	{
		error("[bam_writer_sort_file] ByName must be bool value!\n");
		return R_NilValue;
	}
	const char *filename=CHAR(STRING_ELT(pFilename,0));
	const char *prefix=CHAR(STRING_ELT(pPrefix,0));
	size_t max_mem=*REAL(pMaxMem);
	_Bool sort_by_name =*(LOGICAL(AS_LOGICAL(pByName)));
	if(sort_by_name)
	{
		bam_sort_core_ext(1, filename, prefix, max_mem, 0);
	}
	else
	{
		bam_sort_core_ext(0, filename, prefix, max_mem, 0);
	}
	return R_NilValue;
}


SEXP bam_reader_get_header(SEXP pReader)
{
	if(TYPEOF(pReader)!=EXTPTRSXP)
	{
		error("[bam_reader_get_header] No external pointer!\n");
		return R_NilValue;
	}

	samfile_t *reader=(samfile_t*)(R_ExternalPtrAddr(pReader));
	SEXP ptr;
	PROTECT(ptr=R_MakeExternalPtr((void*)clone_bam_header(reader->header),R_NilValue,R_NilValue));
	R_RegisterCFinalizer(ptr,finalize_bam_header);
	UNPROTECT(1);
	return ptr;
}

///////////////////////////////////////////////////////////////////////////////////////////////////
// GapList
///////////////////////////////////////////////////////////////////////////////////////////////////

static void finalize_gap_list(SEXP ptr)
{
	if(TYPEOF(ptr)!=EXTPTRSXP)
	{
		error("[finalize_gap_list] No external pointer!\n");
		return;
	}
	gap_list *l=(gap_list*)(R_ExternalPtrAddr(ptr));
	destroy_gap_list(l);
	l=NULL;
	R_SetExternalPtrAddr(ptr,NULL);
}

SEXP create_gap_list()
{
	gap_list *l=init_gap_list();
	SEXP list;
	PROTECT(list=R_MakeExternalPtr((void*)(l),R_NilValue,R_NilValue));
	R_RegisterCFinalizer(list,finalize_gap_list);
	UNPROTECT(1);
	return list;
}

static int gap_fetch_func(const bam1_t *b, void *data)
{
	gap_list *l=(gap_list*)data;
	list_gaps(l,b);
	return 0;
}

SEXP gap_list_fetch(SEXP pReader,SEXP pIndex,SEXP pCoords)
{
	if(TYPEOF(pReader)!=EXTPTRSXP)
	{
		error("[gap_list_fetch] pReader is no external pointer!\n");
		return R_NilValue;
	}
	if(TYPEOF(pIndex)!=EXTPTRSXP)
	{
		error("[gap_list_fetch] pIndex is no external pointer!\n");
		return R_NilValue;
	}
	if(TYPEOF(pCoords)!=REALSXP)
	{
		error("[gap_list_fetch] pCoords is no REAL!\n");
		return R_NilValue;
	}
	if(LENGTH(pCoords)!=3)
	{
		error("[gap_list_fetch] pCoords must contain three values (refid,begin,end)!\n");
		return R_NilValue;
	}

	samfile_t *reader=(samfile_t*)(R_ExternalPtrAddr(pReader));
	bam_index_t *index=(bam_index_t*)(R_ExternalPtrAddr(pIndex));
	if(reader==NULL)
	{
		error("[gap_list_fetch] Reader must not be NULL pointer!\n");
		return R_NilValue;
	}
	if(index==NULL)
	{
		error("[gap_list_fetch] Index must not be NULL pointer!\n");
		return R_NilValue;
	}

	double *pi=REAL(pCoords);
	int refid=(int) pi[0];
	int begin=(int) pi[1];
	int end=(int) pi[2];

	if(refid<0 || refid >=(reader->header->n_targets))
	{
		error("[gap_list_fetch] refid out of range!\n");
		return R_NilValue;
	}
	if(begin<0 || begin>=end || end>(reader->header->target_len[refid]))
	{
		error("[gap_list_fetch] Begin or end out of range!\n");
		return R_NilValue;
	}

	gap_list *l=init_gap_list();
    bam_fetch(reader->x.bam, index, refid, begin, end, (void*)l, gap_fetch_func);
    SEXP list;
	PROTECT(list=R_MakeExternalPtr( (void*)(l),R_NilValue,R_NilValue));
	R_RegisterCFinalizer(list,finalize_gap_list);
	UNPROTECT(1);
	Rprintf("[gap_list_fetch] Fetched list of size %i.\n",l->size);
	return list;
}

SEXP gap_list_get_size(SEXP pGapList)
{
	if(TYPEOF(pGapList)!=EXTPTRSXP)
		error("[gap_list_get_size] no external pointer!");

	gap_list *l=(gap_list*)(R_ExternalPtrAddr(pGapList));
	SEXP ans;
	PROTECT(ans=allocVector(INTSXP,1));
	INTEGER(ans)[0]=(l->size);
	UNPROTECT(1);
	return ans;
}

SEXP get_gap_list_df(SEXP pGapList)
{
	if(TYPEOF(pGapList)!=EXTPTRSXP)
		error("[get_gap_list_df] no external pointer!");

	gap_list *l=(gap_list*)(R_ExternalPtrAddr(pGapList));

	// create data.frame
	int nProtected=0;
	int nCols=8;
	SEXP dflist;
	PROTECT(dflist=allocVector(VECSXP,nCols));
	++nProtected;
	int nRows=(l->size);
	int i;
	unsigned *pos_data;

	// Column 0: refid
	SEXP ref_vector;
	PROTECT(ref_vector=allocVector(INTSXP,nRows));
	++nProtected;

	// Column 1: position
	SEXP pos_vector;
	PROTECT(pos_vector=allocVector(INTSXP,nRows));
	++nProtected;

	// Column 2: left_cigar_len
	SEXP left_cigar_len_vector;
	PROTECT(left_cigar_len_vector=allocVector(INTSXP,nRows));
	++nProtected;

	// Column 3: left_cigar_type
	SEXP left_cigar_type_vector;
	PROTECT(left_cigar_type_vector=allocVector(INTSXP,nRows));
	++nProtected;

	// Column 4: left_stop
	SEXP left_stop_vector;
	PROTECT(left_stop_vector=allocVector(INTSXP,nRows));
	++nProtected;

	// Column 5: right_start
	SEXP right_start_vector;
	PROTECT(right_start_vector=allocVector(INTSXP,nRows));
	++nProtected;

	// Column 6: right_cigar_len
	SEXP right_cigar_len_vector;
	PROTECT(right_cigar_len_vector=allocVector(INTSXP,nRows));
	++nProtected;

	// Column 7: right_cigar_type
	SEXP right_cigar_type_vector;
	PROTECT(right_cigar_type_vector=allocVector(INTSXP,nRows));
	++nProtected;

	// start position
	l->curr_el=l->first_el;

	for(i=0;i<nRows;++i)
	{
		pos_data=l->curr_el->pos_data;
		INTEGER(ref_vector)        		[i]=pos_data[0];
		INTEGER(pos_vector)        		[i]=pos_data[1];
		INTEGER(left_cigar_len_vector)	[i]=pos_data[2];
		INTEGER(left_cigar_type_vector)	[i]=pos_data[3];
		INTEGER(left_stop_vector)  		[i]=pos_data[4];
		INTEGER(right_start_vector)		[i]=pos_data[5]+1;
		INTEGER(right_cigar_len_vector)	[i]=pos_data[6];
		INTEGER(right_cigar_type_vector)[i]=pos_data[7];
		l->curr_el=l->curr_el->next_el;
	}

	SET_VECTOR_ELT(dflist,0,ref_vector);
	SET_VECTOR_ELT(dflist,1,pos_vector);
	SET_VECTOR_ELT(dflist,2,left_cigar_len_vector);
	SET_VECTOR_ELT(dflist,3,left_cigar_type_vector);
	SET_VECTOR_ELT(dflist,4,left_stop_vector);
	SET_VECTOR_ELT(dflist,5,right_start_vector);
	SET_VECTOR_ELT(dflist,6,right_cigar_len_vector);
	SET_VECTOR_ELT(dflist,7,right_cigar_type_vector);

	// Column Names
	SEXP col_names;
	PROTECT(col_names=allocVector(STRSXP,nCols));
	++nProtected;

	SET_STRING_ELT(col_names,0,mkChar("refid"));
	SET_STRING_ELT(col_names,1,mkChar("position"));
	SET_STRING_ELT(col_names,2,mkChar("left_cigar_len"));
	SET_STRING_ELT(col_names,3,mkChar("left_cigar_type"));
	SET_STRING_ELT(col_names,4,mkChar("left_stop"));
	SET_STRING_ELT(col_names,5,mkChar("right_start"));
	SET_STRING_ELT(col_names,6,mkChar("right_cigar_len"));
	SET_STRING_ELT(col_names,7,mkChar("right_cigar_type"));
	setAttrib(dflist,R_NamesSymbol,col_names);

	// Row Names
	SEXP row_names;
    PROTECT(row_names=allocVector(STRSXP,nRows));
    ++nProtected;

	int buf_size=128;
	char *buf=(char*) calloc(buf_size,sizeof(char));
    for(i=0;i<nRows;++i)
    {
    	sprintf(buf,"%i",i);
    	SET_STRING_ELT(row_names,i,mkChar(buf));
    }
    free(buf);
    setAttrib(dflist,R_RowNamesSymbol,row_names);
	setAttrib(dflist,R_ClassSymbol,mkString("data.frame"));
	UNPROTECT(nProtected);
	return dflist;
}

///////////////////////////////////////////////////////////////////////////////////////////////////
// BamRange
///////////////////////////////////////////////////////////////////////////////////////////////////

static void finalize_bam_range(SEXP ptr)
{
	if(TYPEOF(ptr)!=EXTPTRSXP)
	{
		error("[finalize_bam_range] No external pointer!\n");
		return;
	}
	align_list *l=(align_list *)(R_ExternalPtrAddr(ptr));
	destroy_align_list(l);
	l=NULL;
	R_SetExternalPtrAddr(ptr,NULL);
}


SEXP bam_range_init()
{
	align_list *l=init_align_list();
	SEXP list;
	PROTECT(list=R_MakeExternalPtr( (void*)(l),R_NilValue,R_NilValue));
	R_RegisterCFinalizer(list,finalize_bam_range);
	UNPROTECT(1);
	return list;
}

static int range_fetch_func(const bam1_t *b, void *data)
{
	align_list *l=(align_list*)data;
	align_list_push_back(l,b);
	return 0;
}

static int range_fetch_complex_func(const bam1_t *b,void *data)
{
	align_list *l=(align_list*)data;
	if(b->core.n_cigar>1)
		align_list_push_back(l,b);
	return 0;
}

SEXP bam_range_fetch(SEXP pReader,SEXP pIndex,SEXP pCoords,SEXP pComplex)
{
	if(TYPEOF(pReader)!=EXTPTRSXP)
	{
		error("[bam_range_fetch] pReader is no external pointer!\n");
		return R_NilValue;
	}
	if(TYPEOF(pIndex)!=EXTPTRSXP)
	{
		error("[bam_range_fetch] pIndex is no external pointer!\n");
		return R_NilValue;
	}
	if(TYPEOF(pCoords)!=REALSXP)
	{
		error("[bam_range_fetch] pCoords is no REAL!\n");
		return R_NilValue;
	}
	if(LENGTH(pCoords)!=3)
	{
		error("[bam_range_fetch] pCoords must contain three values (refid,begin,end)!\n");
		return R_NilValue;
	}
	if(TYPEOF(pComplex)!=LGLSXP)
	{
		error("[bam_range_fetch] pComplex must be logical!\n");
		return R_NilValue;
	}

	samfile_t *reader=(samfile_t*)(R_ExternalPtrAddr(pReader));
	bam_index_t *index=(bam_index_t*)(R_ExternalPtrAddr(pIndex));
	if(reader==NULL)
	{
		error("[bam_range_fetch] Reader must not be NULL pointer!\n");
		return R_NilValue;
	}
	if(index==NULL)
	{
		error("[bam_range_fetch] Index must not be NULL pointer!\n");
		return R_NilValue;
	}

	double *pi=REAL(pCoords);
	int refid=(int) pi[0];
	int begin=(int) pi[1];
	int end=(int) pi[2];

	if(refid<0 || refid >=(reader->header->n_targets))
	{
		error("[bam_range_fetch] refid out of range!\n");
		return R_NilValue;
	}
	if(begin<0 || begin>=end || end>(reader->header->target_len[refid]))
	{
		error("[bam_range_fetch] Begin or end out of range!\n");
		return R_NilValue;
	}

	align_list *l=init_align_list();

	// Retrieve only complex aligns (nCigar>1) when pComplex is set:
	if(LOGICAL(pComplex)[0]==TRUE)
		bam_fetch(reader->x.bam, index, refid, begin, end, (void*)l, range_fetch_complex_func);
	else
		bam_fetch(reader->x.bam, index, refid, begin, end, (void*)l, range_fetch_func);

	// Wind back
    l->curr_el=NULL;

    SEXP list;
	PROTECT(list=R_MakeExternalPtr( (void*)(l),R_NilValue,R_NilValue));
	R_RegisterCFinalizer(list,finalize_bam_range);
	UNPROTECT(1);
	return list;
}

SEXP bam_range_get_next_align(SEXP pRange)
{
	if(TYPEOF(pRange)!=EXTPTRSXP)
		error("[bam_range_get_next_align] No external pointer!");
	align_list *l=(align_list*)(R_ExternalPtrAddr(pRange));
	bam1_t *align=get_next_align(l);
	if(align==NULL)
		return R_NilValue;

	SEXP ptr;
	PROTECT(ptr=R_MakeExternalPtr((void*)(align),R_NilValue,R_NilValue));
	R_RegisterCFinalizer(ptr,finalize_bam_align);
	UNPROTECT(1);
	return ptr;
}

SEXP bam_range_get_prev_align(SEXP pRange)
{
	if(TYPEOF(pRange)!=EXTPTRSXP)
		error("[bam_range_get_prev_align] No external pointer!\n");
	align_list *l=(align_list*)(R_ExternalPtrAddr(pRange));
	bam1_t *align=get_prev_align(l);
	if(align==NULL)
		return R_NilValue;

	SEXP ptr;
	PROTECT(ptr=R_MakeExternalPtr((void*)(align),R_NilValue,R_NilValue));
	R_RegisterCFinalizer(ptr,finalize_bam_align);
	UNPROTECT(1);
	return ptr;
}

SEXP bam_range_step_next_align(SEXP pRange)
{
	if(TYPEOF(pRange)!=EXTPTRSXP)
		error("[bam_range_step_next_align] No external pointer!\n");
	pp_curr_align((align_list*)(R_ExternalPtrAddr(pRange)));
	return R_NilValue;
}

SEXP bam_range_step_prev_align(SEXP pRange)
{
	if(TYPEOF(pRange)!=EXTPTRSXP)
		error("[bam_range_step_prev_align] No external pointer!\n");
	mm_curr_align((align_list*)(R_ExternalPtrAddr(pRange)));
	return R_NilValue;
}


SEXP bam_range_write_current_align(SEXP pRange,SEXP pAlign)
{
	if(TYPEOF(pRange)!=EXTPTRSXP)
	{
		error("[bam_range_push_back] pRange is no external pointer!\n");
		return R_NilValue;
	}
	if(TYPEOF(pAlign)!=EXTPTRSXP)
	{
		error("[bam_range_push_back] pAlign is no external pointer!\n");
		return R_NilValue;
	}
	write_current_align((align_list*)(R_ExternalPtrAddr(pRange)),(bam1_t*) (R_ExternalPtrAddr(pAlign)));
	return R_NilValue;
}

SEXP bam_range_insert_past_curr_align(SEXP pRange,SEXP pAlign)
{
	if(TYPEOF(pRange)!=EXTPTRSXP)
	{
		error("[bam_range_push_back] pRange is no external pointer!\n");
		return R_NilValue;
	}
	if(TYPEOF(pAlign)!=EXTPTRSXP)
	{
		error("[bam_range_push_back] pAlign is no external pointer!\n");
		return R_NilValue;
	}
	insert_past_curr_align((align_list*)(R_ExternalPtrAddr(pRange)),(bam1_t*) (R_ExternalPtrAddr(pAlign)));
	return R_NilValue;
}

SEXP bam_range_insert_pre_curr_align(SEXP pRange,SEXP pAlign)
{
	if(TYPEOF(pRange)!=EXTPTRSXP)
	{
		error("[bam_range_push_back] pRange is no external pointer!\n");
		return R_NilValue;
	}
	if(TYPEOF(pAlign)!=EXTPTRSXP)
	{
		error("[bam_range_push_back] pAlign is no external pointer!\n");
		return R_NilValue;
	}
	insert_pre_curr_align((align_list*)(R_ExternalPtrAddr(pRange)),(bam1_t*) (R_ExternalPtrAddr(pAlign)));
	return R_NilValue;
}

SEXP bam_range_mv_curr_align(SEXP pSrc, SEXP pTarget)
{
	// Moves current align in src to end of target list
	// and Moves current align in src to next align
	if(TYPEOF(pSrc)!=EXTPTRSXP)
	{
		error("[bam_range_mv_curr_align] pSrc is no external pointer!\n");
		return R_NilValue;
	}
	if(TYPEOF(pTarget)!=EXTPTRSXP)
	{
		error("[bam_range_mv_curr_align] pTarget is no external pointer!\n");
		return R_NilValue;
	}

	align_list_mv_curr_elem((align_list*)(R_ExternalPtrAddr(pSrc)),(align_list*)(R_ExternalPtrAddr(pTarget)));
	return R_NilValue;
}

SEXP bam_range_get_align_df(SEXP pRange)
{
	if(TYPEOF(pRange)!=EXTPTRSXP)
		error("[bam_range_get_align_df] no external pointer!");

	align_list *l=(align_list*)(R_ExternalPtrAddr(pRange));
	// Read Adress of current align->for reconstitution at the end
	align_element *e=l->curr_el;
	wind_back(l);
	bam1_t *align;

	// create data.frame
	int nProtected=0;
	int nCols=7;
	SEXP dflist;
	PROTECT(dflist=allocVector(VECSXP,nCols));
	++nProtected;
	int nRows=(l->size);
	int i,j;

	// Column 0: refid
	SEXP ref_vector;
	PROTECT(ref_vector=allocVector(INTSXP,nRows));
	++nProtected;

	// Column 1: position
	SEXP pos_vector;
	PROTECT(pos_vector=allocVector(INTSXP,nRows));
	++nProtected;

	// Column 2: nCigar
	SEXP nCigar_vector;
	PROTECT(nCigar_vector=allocVector(INTSXP,nRows));
	++nProtected;

	// Column 3: cigar
	SEXP cig_vector;
	PROTECT(cig_vector=allocVector(STRSXP,nRows));
	++nProtected;

	// Column 4: flag
	SEXP flag_vector;
	PROTECT(flag_vector=allocVector(INTSXP,nRows));
	++nProtected;

	// Column 5: seq
	SEXP seq_vector;
	PROTECT(seq_vector=allocVector(STRSXP,nRows));
	++nProtected;

	// Column 6: qual
	SEXP qual_vector;
	PROTECT(qual_vector=allocVector(STRSXP,nRows));
	++nProtected;

	// seq+cigar
	unsigned char *raw_seq;
	int32_t seq_len;
	int buf_size=1024;
	char *buf=(char*) calloc(buf_size,sizeof(char));
	uint8_t *quals;

	for(i=0;i<nRows;++i)
	{
		align=get_next_align(l);
		INTEGER(ref_vector)[i]=(align->core.tid);
		INTEGER(pos_vector)[i]=(align->core.pos);
		INTEGER(nCigar_vector)[i]=(align->core.n_cigar);

		/////////////////////////////////////////
		// Cigar String
		if(cigar2str(buf,align)==0)
		{
			error("[bam_align_get_align_df] Cigar error!\n");
			return R_NilValue;
		}
		SET_STRING_ELT(cig_vector,i,mkChar(buf));
		//Rprintf("%i\t%s\n",i,buf);
		clear_buf(buf,buf_size);
		/////////////////////////////////////////

		INTEGER(flag_vector)[i]=(align->core.flag);
		/////////////////////////////////////////
		// seq
		seq_len=align->core.l_qseq;
		if(seq_len>buf_size)
		{
			buf_size=2*(seq_len+1);
			free(buf);
			buf= (char*) calloc(buf_size,sizeof(char));
		}
		raw_seq=bam1_seq(align);
		for(j=0;j<seq_len;++j)
			buf[j]=bam_nt16_rev_table[bam1_seqi(raw_seq,j)];
		buf[j]=0;

		SET_STRING_ELT(seq_vector,i,mkChar(buf));
		////////////////////////////////////////
		// quals
			quals=bam1_qual(align);
			for(j=0;j<seq_len;++j)
				buf[j]=(char) (quals[j]+33);
			buf[j]=0;
			SET_STRING_ELT(qual_vector,i,mkChar(buf));
	}

	// Reset curr_el pointer
	l->curr_el=e;

	SET_VECTOR_ELT(dflist,0,ref_vector);
	SET_VECTOR_ELT(dflist,1,pos_vector);
	SET_VECTOR_ELT(dflist,2,nCigar_vector);
	SET_VECTOR_ELT(dflist,3,cig_vector);
	SET_VECTOR_ELT(dflist,4,flag_vector);
	SET_VECTOR_ELT(dflist,5,seq_vector);
	SET_VECTOR_ELT(dflist,6,qual_vector);

	// Column Names
	SEXP col_names;
	PROTECT(col_names=allocVector(STRSXP,nCols));
	++nProtected;

	SET_STRING_ELT(col_names,0,mkChar("refid"));
	SET_STRING_ELT(col_names,1,mkChar("position"));
	SET_STRING_ELT(col_names,2,mkChar("nCigar"));
	SET_STRING_ELT(col_names,3,mkChar("cigar"));
	SET_STRING_ELT(col_names,4,mkChar("flag"));
	SET_STRING_ELT(col_names,5,mkChar("seq"));
	SET_STRING_ELT(col_names,6,mkChar("qual"));
	setAttrib(dflist,R_NamesSymbol,col_names);

	SEXP row_names;
    PROTECT(row_names=allocVector(STRSXP,nRows));
    ++nProtected;

    for(i=0;i<nRows;++i)
    {
    	sprintf(buf,"%i",i);
    	SET_STRING_ELT(row_names,i,mkChar(buf));
    }
    free(buf);
    setAttrib(dflist,R_RowNamesSymbol,row_names);
	setAttrib(dflist,R_ClassSymbol,mkString("data.frame"));
	UNPROTECT(nProtected);
	return dflist;
}

SEXP bam_range_write(SEXP pWriter,SEXP pRange)
{
	if(TYPEOF(pRange)!=EXTPTRSXP)
	{
		error("[bam_range_write] pRange is no external pointer!\n");
		return R_NilValue;
	}
	if(TYPEOF(pWriter)!=EXTPTRSXP)
	{
		error("[bam_range_write] pWriter is no external pointer!\n");
		return R_NilValue;
	}

	unsigned long bytes_written=0;
	unsigned long range_size, i;

	samfile_t *writer=(samfile_t*) R_ExternalPtrAddr(pWriter);
	align_list *l=(align_list*)(R_ExternalPtrAddr(pRange));
	bam1_t *align;

	range_size=l->size;
	// For restoration
	align_element *e=l->curr_el;
	wind_back(l);

	for(i=0;i<range_size;++i)
	{
		align=get_next_align(l);
		bytes_written+=samwrite(writer,align);
	}
	// restore curr_el
	l->curr_el=e;

	SEXP ans;
	PROTECT(ans=allocVector(INTSXP,1));
	INTEGER(ans)[0]=bytes_written;
	UNPROTECT(1);
	return ans;
}

SEXP bam_range_wind_back(SEXP pRange)
{
	if(TYPEOF(pRange)!=EXTPTRSXP)
	{
		error("[bam_range_wind_back] pRange is no external pointer!\n");
		return R_NilValue;
	}
	align_list *l=(align_list*)(R_ExternalPtrAddr(pRange));
	wind_back(l);
	return R_NilValue;
}

SEXP bam_range_get_size(SEXP pRange)
{
	if(TYPEOF(pRange)!=EXTPTRSXP)
	{
		error("[bam_range_get_size] pRange is no external pointer!\n");
		return R_NilValue;
	}
	align_list *l=(align_list*)(R_ExternalPtrAddr(pRange));
	SEXP ans;
	PROTECT(ans=allocVector(INTSXP,1));
	INTEGER(ans)[0]=(l->size);
	UNPROTECT(1);
	return ans;
}

SEXP bam_range_push_back(SEXP pRange,SEXP pAlign)
{
	if(TYPEOF(pRange)!=EXTPTRSXP)
	{
		error("[bam_range_push_back] pRange is no external pointer!\n");
		return R_NilValue;
	}
	if(TYPEOF(pAlign)!=EXTPTRSXP)
	{
		error("[bam_range_push_back] pAlign is no external pointer!\n");
		return R_NilValue;
	}
	align_list_push_back((align_list*)(R_ExternalPtrAddr(pRange)),(bam1_t*) (R_ExternalPtrAddr(pAlign)));
	return R_NilValue;
}

SEXP bam_range_push_front(SEXP pRange,SEXP pAlign)
{
	if(TYPEOF(pRange)!=EXTPTRSXP)
	{
		error("[bam_range_push_front] pRange is no external pointer!\n");
		return R_NilValue;
	}
	if(TYPEOF(pAlign)!=EXTPTRSXP)
	{
		error("[bam_range_push_front] pAlign is no external pointer!\n");
		return R_NilValue;
	}
	align_list_push_front((align_list*)(R_ExternalPtrAddr(pRange)),(bam1_t*) (R_ExternalPtrAddr(pAlign)));
	return R_NilValue;
}

SEXP bam_range_pop_back(SEXP pRange)
{
	if(TYPEOF(pRange)!=EXTPTRSXP)
	{
		error("[bam_range_pop_back] pRange is no external pointer!\n");
		return R_NilValue;
	}
	align_list_pop_back((align_list*)(R_ExternalPtrAddr(pRange)));
	return R_NilValue;
}
SEXP bam_range_pop_front(SEXP pRange)
{
	if(TYPEOF(pRange)!=EXTPTRSXP)
	{
		error("[bam_range_pop_front] pRange is no external pointer!\n");
		return R_NilValue;
	}
	align_list_pop_front((align_list*)(R_ExternalPtrAddr(pRange)));
	return R_NilValue;
}


///////////////////////////////////////////////////////////////////////////////////////////////
// bam_header
///////////////////////////////////////////////////////////////////////////////////////////////

static void finalize_bam_header(SEXP ptr)
{
	if(TYPEOF(ptr)!=EXTPTRSXP)
	{
		error("[finalize_bam_header] No external pointer!");
		return;
	}
	bam_header_t* header=(bam_header_t*)(R_ExternalPtrAddr(ptr));
	if(header)
		bam_header_destroy(header);
	R_SetExternalPtrAddr(ptr,NULL);
}

SEXP init_bam_header(SEXP pHeaderText)
{
	if(TYPEOF(pHeaderText)==NILSXP)
		return R_NilValue;

	if(TYPEOF(pHeaderText)!=STRSXP)
	{
		error("[init_bam_header] Header Text must be a string.\n");
		return R_NilValue;
	}
	bam_header_t *h=(bam_header_t*)calloc(1, sizeof(bam_header_t));

	// Copy header text
	const char* header_text=CHAR(STRING_ELT(pHeaderText,0));
	h->l_text=strlen(header_text);
	char *txt=(char*) calloc(1,(h->l_text)+1);
	strncpy(txt,CHAR(STRING_ELT(pHeaderText,0)),h->l_text);
	h->text=txt;

	sam_header_parse(h);
	bam_init_header_hash(h);

	SEXP ptr;
	PROTECT(ptr=R_MakeExternalPtr((void*)h,R_NilValue,R_NilValue));
	R_RegisterCFinalizer(ptr,finalize_bam_header);
	UNPROTECT(1);
	return ptr;
}

SEXP bam_header_get_header_text(SEXP pHeader)
{
	if(TYPEOF(pHeader)!=EXTPTRSXP)
	{
		error("[bam_header_get_header_text] No external pointer!");
		return R_NilValue;
	}

	bam_header_t* header=(bam_header_t*)(R_ExternalPtrAddr(pHeader));
	SEXP ans;
	PROTECT(ans=Rf_allocVector(STRSXP,1));
	SET_STRING_ELT(ans,0,mkChar(header->text));
	UNPROTECT(1);
	return ans;
}
///////////////////////////////////////////////////////////////////////////////////////////////
// bam_writer
///////////////////////////////////////////////////////////////////////////////////////////////

static void finalize_bam_writer(SEXP ptr)
{
	if(TYPEOF(ptr)!=EXTPTRSXP)
	{
		error("[finalize_bam_writer] No external pointer!");
		return;
	}
	samfile_t *writer=(samfile_t*)(R_ExternalPtrAddr(ptr));
	if(writer)
	{
		samclose(writer);
		R_SetExternalPtrAddr(ptr,NULL);
		Rprintf("[bamWriter] finalized.\n");
	}
}

SEXP bam_writer_open(SEXP pHeader,SEXP pFilename)
{
	if(TYPEOF(pHeader)!=EXTPTRSXP)
	{
		error("[bam_writer_open] pHeader no external pointer!\n");
		return R_NilValue;
	}
	if(TYPEOF(pFilename)!=STRSXP)
	{
		error("[bam_writer_open] pFilename no string!\n");
		return R_NilValue;
	}

	bam_header_t *header=(bam_header_t*) (R_ExternalPtrAddr(pHeader));
	samfile_t    *writer=samopen(CHAR(STRING_ELT(pFilename,0)),"wb",header);
	if(!writer)
		error("[bam_writer_open] samopen returned NULL pointer!\n");

	SEXP ptr;
	PROTECT(ptr=R_MakeExternalPtr( (void*) (writer),R_NilValue,R_NilValue));
	R_RegisterCFinalizer(ptr,finalize_bam_writer);
	UNPROTECT(1);
	return ptr;
}

SEXP bam_reader_open_writer(SEXP pReader,SEXP pFilename)
{
	if(TYPEOF(pReader)!=EXTPTRSXP)
	{
		error("[bam_writer_open] pReader no external pointer!\n");
		return R_NilValue;
	}
	if(TYPEOF(pFilename)!=STRSXP)
	{
		error("[bam_writer_open] pFilename no string!\n");
		return R_NilValue;
	}
	samfile_t *reader=(samfile_t*) (R_ExternalPtrAddr(pReader));
	samfile_t *writer=samopen(CHAR(STRING_ELT(pFilename,0)),"wb",reader->header);

	SEXP ptr;
	PROTECT(ptr=R_MakeExternalPtr( (void*) (writer),R_NilValue,R_NilValue));
	R_RegisterCFinalizer(ptr,finalize_bam_writer);
	UNPROTECT(1);
	return ptr;
}

SEXP bam_writer_save_align(SEXP pWriter, SEXP pAlign)
{
	if(TYPEOF(pWriter)!=EXTPTRSXP)
	{
		error("[bam_writer_save_align] No external pointer!\n");
		return R_NilValue;
	}
	if(TYPEOF(pAlign)!=EXTPTRSXP)
	{
		error("[bam_writer_save_align] No external pointer!\n");
		return R_NilValue;
	}
	samfile_t *writer=(samfile_t*)R_ExternalPtrAddr(pWriter);
	bam1_t *align=(bam1_t*)R_ExternalPtrAddr(pAlign);

	SEXP ans;
	PROTECT(ans=allocVector(INTSXP,1));
	INTEGER(ans)[0]=samwrite(writer,align);
	UNPROTECT(1);
	return ans;
}

SEXP bam_writer_close(SEXP pWriter)
{
	if(TYPEOF(pWriter)!=EXTPTRSXP)
	{
		error("[bam_writer_close] No exteranl pointer!\n");
		return R_NilValue;
	}
	samfile_t *writer= (samfile_t*) (R_ExternalPtrAddr(pWriter));
	samclose(writer);
	R_SetExternalPtrAddr(pWriter,NULL);
	Rprintf("[bamWriter] closed.\n");
	return R_NilValue;
}

///////////////////////////////////////////////////////////////////////////////////////////////
// bam_align
///////////////////////////////////////////////////////////////////////////////////////////////

static void finalize_bam_align(SEXP pAlign)
{
	if(TYPEOF(pAlign)!=EXTPTRSXP)
		error("[finalize_bam_align] No external pointer!");
	bam1_t *align= (bam1_t*)(R_ExternalPtrAddr(pAlign));
	bam_destroy1(align);	// checks for >0!
}

SEXP bam_align_get_name(SEXP pAlign)
{
	if(TYPEOF(pAlign)!=EXTPTRSXP)
		error("[bam_align_get_name] No external pointer!");
	bam1_t *align= (bam1_t*)(R_ExternalPtrAddr(pAlign));

	SEXP ans;
	PROTECT(ans=Rf_allocVector(STRSXP,1));
	SET_STRING_ELT(ans,0,mkChar(bam1_qname(align)));
	UNPROTECT(1);
	return ans;
}

SEXP bam_align_get_refid(SEXP pAlign)
{
	if(TYPEOF(pAlign)!=EXTPTRSXP)
		error("[bam_align_getRefID] No external pointer!");
	bam1_t *align= (bam1_t*) (R_ExternalPtrAddr(pAlign));

	SEXP ans;
	PROTECT(ans=Rf_allocVector(INTSXP,1));
	INTEGER(ans)[0]=align->core.tid;
	UNPROTECT(1);
	return ans;
}

SEXP bam_align_get_position(SEXP pAlign)
{
	if(TYPEOF(pAlign)!=EXTPTRSXP)
		error("[bam_align_get_position] No external pointer!");
	bam1_t *align= (bam1_t*) (R_ExternalPtrAddr(pAlign));

	SEXP ans;
	PROTECT(ans=Rf_allocVector(INTSXP,1));
	INTEGER(ans)[0]=align->core.pos;
	UNPROTECT(1);
	return ans;
}

SEXP bam_align_get_nCigar(SEXP pAlign)
{
	if(TYPEOF(pAlign)!=EXTPTRSXP)
		error("[bam_align_get_nCigar] No external pointer!");
	bam1_t *align= (bam1_t*) (R_ExternalPtrAddr(pAlign));
	SEXP ans;
	PROTECT(ans=Rf_allocVector(INTSXP,1));
	INTEGER(ans)[0]=align->core.n_cigar;
	UNPROTECT(1);
	return ans;
}

SEXP bam_align_get_cigar_df(SEXP pAlign)
{
	if(TYPEOF(pAlign)!=EXTPTRSXP)
		error("[bam_align_get_cigar_df] No external pointer!");
	bam1_t *align= (bam1_t*) (R_ExternalPtrAddr(pAlign));
		// create data.frame
	int nProtected=0;
	int nCols=2;
	SEXP dflist;
	PROTECT(dflist=allocVector(VECSXP,nCols));
	++nProtected;
	int nRows=align->core.n_cigar;
	int i;

	// Column 0: Length
	SEXP Length_vector;
	PROTECT(Length_vector=allocVector(INTSXP,nRows));
	++nProtected;

	// Column 1: Type
	SEXP Type_vector;
	PROTECT(Type_vector=allocVector(STRSXP,nRows));
	++nProtected;

	uint32_t *cigar=bam1_cigar(align);
	for(i=0;i<nRows;++i)
	{
		if((cigar[i]&BAM_CIGAR_MASK)>=strlen(CIGAR_TYPES))
		{
			error("[bam_align_getCigar_df] Cigar_type not in defined range!");
			return R_NilValue;
		}
		INTEGER(Length_vector)[i]=cigar[i] >> BAM_CIGAR_SHIFT;
		SET_STRING_ELT(Type_vector,i,mkCharLen(CIGAR_TYPES+(cigar[i]&BAM_CIGAR_MASK),1));
	}

	SET_VECTOR_ELT(dflist,0,Length_vector);
	SET_VECTOR_ELT(dflist,1,Type_vector);

	// Column Names
	SEXP col_names;
	PROTECT(col_names=allocVector(STRSXP,nCols));
	++nProtected;

	SET_STRING_ELT(col_names,0,mkChar("Length"));
	SET_STRING_ELT(col_names,1,mkChar("Type"));
	setAttrib(dflist,R_NamesSymbol,col_names);

	SEXP row_names;
    PROTECT(row_names=allocVector(STRSXP,nRows));
    ++nProtected;

    char c[20];
    for(i=0;i<nRows;++i)
    {
    	sprintf(c,"%i",i);
    	SET_STRING_ELT(row_names,i,mkChar(c));
    }
    setAttrib(dflist,R_RowNamesSymbol,row_names);

	setAttrib(dflist,R_ClassSymbol,mkString("data.frame"));
	UNPROTECT(nProtected);
	return dflist;
}

SEXP bam_align_get_mate_refid(SEXP pAlign)
{
	if(TYPEOF(pAlign)!=EXTPTRSXP)
		error("[bam_align_get_mate_refid] No external pointer!");
	bam1_t *align= (bam1_t*) (R_ExternalPtrAddr(pAlign));

	SEXP ans;
	PROTECT(ans=Rf_allocVector(INTSXP,1));
	INTEGER(ans)[0]=align->core.mtid;
	UNPROTECT(1);
	return ans;
}

SEXP bam_align_get_mate_position(SEXP pAlign)
{
	if(TYPEOF(pAlign)!=EXTPTRSXP)
		error("[bam_align_get_mate_position] No external pointer!");
	bam1_t *align= (bam1_t*) (R_ExternalPtrAddr(pAlign));

	SEXP ans;
	PROTECT(ans=Rf_allocVector(INTSXP,1));
	INTEGER(ans)[0]=align->core.mpos;
	UNPROTECT(1);
	return ans;
}

SEXP bam_align_get_insert_size(SEXP pAlign)
{
	if(TYPEOF(pAlign)!=EXTPTRSXP)
		error("[bam_align_get_insert_size] No external pointer!");
	bam1_t *align= (bam1_t*) (R_ExternalPtrAddr(pAlign));

	SEXP ans;
	PROTECT(ans=Rf_allocVector(INTSXP,1));
	INTEGER(ans)[0]=align->core.isize;
	UNPROTECT(1);
	return ans;
}

SEXP bam_align_get_map_quality(SEXP pAlign)
{
	if(TYPEOF(pAlign)!=EXTPTRSXP)
		error("[bam_align_get_map_quality] No external pointer!");
	bam1_t *align= (bam1_t*) (R_ExternalPtrAddr(pAlign));

	SEXP ans;
	PROTECT(ans=Rf_allocVector(INTSXP,1));
	INTEGER(ans)[0]=align->core.qual;
	UNPROTECT(1);
	return ans;
}

SEXP bam_align_get_segment_sequence(SEXP pAlign)
{
	if(TYPEOF(pAlign)!=EXTPTRSXP)
	{
		error("[bam_align_get_segment_sequence] No external pointer!");
		return R_NilValue;
	}

	// Extract char* sequence with samtools
	bam1_t *align= (bam1_t*) (R_ExternalPtrAddr(pAlign));
	int32_t seq_len=align->core.l_qseq;
	char *seq= (char*) calloc(seq_len+1,sizeof(char));
	unsigned char *raw_seq=bam1_seq(align);
	int32_t i;
	for(i=0;i<seq_len;++i)
		seq[i]=bam_nt16_rev_table[bam1_seqi(raw_seq,i)];
	seq[i]=0;

	SEXP ans;
	PROTECT(ans=allocVector(STRSXP,1));
	SET_STRING_ELT(ans,0,mkChar(seq));
	free(seq);
	UNPROTECT(1);
	return ans;
}

SEXP bam_align_get_qualities(SEXP pAlign)
{
	if(TYPEOF(pAlign)!=EXTPTRSXP)
	{
		error("[bam_align_get_qualities] No external pointer!");
		return R_NilValue;
	}
	bam1_t *align= (bam1_t*) (R_ExternalPtrAddr(pAlign));

	SEXP ans;
	PROTECT(ans=allocVector(STRSXP,1));
	SET_STRING_ELT(ans,0,mkChar((char*)bam1_qual(align)));
	UNPROTECT(1);
	return ans;
}

///////////////////////////////////////////////////////////////////////////////////////////////
// Alignment flags

SEXP bam_align_is_paired(SEXP pAlign)
{
	if(TYPEOF(pAlign)!=EXTPTRSXP)
	{
		error("[bam_align_is_paired] No external pointer!");
		return R_NilValue;
	}
	bam1_t *align= (bam1_t*) (R_ExternalPtrAddr(pAlign));
	SEXP ans;
	PROTECT(ans=allocVector(LGLSXP,1));
	LOGICAL(ans)[0]=align->core.flag & BAM_FPAIRED;
	UNPROTECT(1);
	return ans;
}

SEXP bam_align_mapped_in_proper_pair(SEXP pAlign)
{
	if(TYPEOF(pAlign)!=EXTPTRSXP)
	{
		error("[bam_align_mapped_in_proper_pair] No external pointer!");
		return R_NilValue;
	}
	bam1_t *align=(bam1_t*) (R_ExternalPtrAddr(pAlign));
	SEXP ans;
	PROTECT(ans=allocVector(LGLSXP,1));
	LOGICAL(ans)[0]=align->core.flag & BAM_FPROPER_PAIR;
	UNPROTECT(1);
	return ans;
}

SEXP bam_align_is_unmapped(SEXP pAlign)
{
	if(TYPEOF(pAlign)!=EXTPTRSXP)
	{
		error("[bam_align_is_unmapped] No external pointer!");
		return R_NilValue;
	}
	bam1_t *align=(bam1_t*) (R_ExternalPtrAddr(pAlign));

	SEXP ans;
	PROTECT(ans=allocVector(LGLSXP,1));
	LOGICAL(ans)[0]=align->core.flag & BAM_FUNMAP;
	UNPROTECT(1);
	return ans;
}

SEXP bam_align_mate_is_unmapped(SEXP pAlign)
{
	if(TYPEOF(pAlign)!=EXTPTRSXP)
	{
		error("[bam_align_mate_is_unmapped] No external pointer!");
		return R_NilValue;
	}
	bam1_t *align=(bam1_t*) (R_ExternalPtrAddr(pAlign));

	SEXP ans;
	PROTECT(ans=allocVector(LGLSXP,1));
	LOGICAL(ans)[0]=align->core.flag & BAM_FUNMAP;
	UNPROTECT(1);
	return ans;
}

SEXP bam_align_strand_reverse(SEXP pAlign)
{
	if(TYPEOF(pAlign)!=EXTPTRSXP)
	{
		error("[bam_align_strand_reverse] No external pointer!");
		return R_NilValue;
	}
	bam1_t *align=(bam1_t*) (R_ExternalPtrAddr(pAlign));

	SEXP ans;
	PROTECT(ans=allocVector(LGLSXP,1));
	LOGICAL(ans)[0]=bam1_strand(align);
	UNPROTECT(1);
	return ans;
}

SEXP bam_align_mate_strand_reverse(SEXP pAlign)
{
	if(TYPEOF(pAlign)!=EXTPTRSXP)
	{
		error("[bam_align_strand_reverse] No external pointer!");
		return R_NilValue;
	}
	bam1_t *align=(bam1_t*)(R_ExternalPtrAddr(pAlign));

	SEXP ans;
	PROTECT(ans=allocVector(LGLSXP,1));
	LOGICAL(ans)[0]=bam1_mstrand(align);
	UNPROTECT(1);
	return ans;
}

SEXP bam_align_is_first_in_pair(SEXP pAlign)
{
	if(TYPEOF(pAlign)!=EXTPTRSXP)
	{
		error("[bam_align_is_first_in_pair] No external pointer!");
		return R_NilValue;
	}
	bam1_t *align=(bam1_t*)(R_ExternalPtrAddr(pAlign));

	SEXP ans;
	PROTECT(ans=allocVector(LGLSXP,1));
	LOGICAL(ans)[0]=align->core.flag & BAM_FREAD1;
	UNPROTECT(1);
	return ans;
}

SEXP bam_align_is_second_in_pair(SEXP pAlign)
{
	if(TYPEOF(pAlign)!=EXTPTRSXP)
	{
		error("[bam_align_is_second_in_pair] No external pointer!");
		return R_NilValue;
	}
	bam1_t *align=(bam1_t*)(R_ExternalPtrAddr(pAlign));

	SEXP ans;
	PROTECT(ans=allocVector(LGLSXP,1));
	LOGICAL(ans)[0]=align->core.flag & BAM_FREAD1;
	UNPROTECT(1);
	return ans;
}

SEXP bam_align_is_secondary_align(SEXP pAlign)
{
	if(TYPEOF(pAlign)!=EXTPTRSXP)
	{
		error("[bam_align_is_secondary_align] No external pointer!");
		return R_NilValue;
	}
	bam1_t *align=(bam1_t*)(R_ExternalPtrAddr(pAlign));

	SEXP ans;
	PROTECT(ans=allocVector(LGLSXP,1));
	LOGICAL(ans)[0]=align->core.flag & BAM_FSECONDARY;
	UNPROTECT(1);
	return ans;
}

SEXP bam_align_fail_qc(SEXP pAlign)
{
	if(TYPEOF(pAlign)!=EXTPTRSXP)
	{
		error("[bam_align_fail_qc] No external pointer!");
		return R_NilValue;
	}
	bam1_t *align=(bam1_t*)(R_ExternalPtrAddr(pAlign));

	SEXP ans;
	PROTECT(ans=allocVector(LGLSXP,1));
	LOGICAL(ans)[0]=align->core.flag & BAM_FQCFAIL;
	UNPROTECT(1);
	return ans;
}

SEXP bam_align_is_pcr_or_optical_dup(SEXP pAlign)
{
	if(TYPEOF(pAlign)!=EXTPTRSXP)
	{
		error("[bam_align_is_pcr_or_optical_dup] No external pointer!");
		return R_NilValue;
	}
	bam1_t *align=(bam1_t*)(R_ExternalPtrAddr(pAlign));

	SEXP ans;
	PROTECT(ans=allocVector(LGLSXP,1));
	LOGICAL(ans)[0]=align->core.flag & BAM_FDUP;
	UNPROTECT(1);
	return ans;
}

SEXP bam_align_get_flag(SEXP pAlign)
{
	if(TYPEOF(pAlign)!=EXTPTRSXP)
	{
		error("[bam_align_get_flag] No external pointer!");
		return R_NilValue;
	}
	bam1_t *align=(bam1_t*)(R_ExternalPtrAddr(pAlign));

	SEXP ans;
	PROTECT(ans=allocVector(INTSXP,1));
	INTEGER(ans)[0]=align->core.flag;
	UNPROTECT(1);
	return ans;
}

///////////////////////////////////////////////////////////////////////////////////////////////
// Writing accessors
SEXP bam_align_set_is_paired(SEXP pAlign, SEXP val)
{
	if(TYPEOF(pAlign)!=EXTPTRSXP)
	{
		error("[bam_align_set_is_paired] No external pointer!");
		return R_NilValue;
	}
	if(TYPEOF(val)!=LGLSXP)
	{
		error("[bam_align_set_is_paired] No bool value!");
		return R_NilValue;
	}
	bam1_t *align=(bam1_t*)(R_ExternalPtrAddr(pAlign));
	set_flag(align,*(LOGICAL(AS_LOGICAL(val))),BAM_FPAIRED);
	return R_NilValue;
}

SEXP bam_align_set_mapped_in_proper_pair(SEXP pAlign, SEXP val)
{
	if(TYPEOF(pAlign)!=EXTPTRSXP)
	{
		error("[bam_align_set_mapped_in_proper_pair] No external pointer!");
		return R_NilValue;
	}
	if(TYPEOF(val)!=LGLSXP)
	{
		error("[bam_align_set_mapped_in_proper_pair] No bool value!");
		return R_NilValue;
	}
	bam1_t *align=(bam1_t*)(R_ExternalPtrAddr(pAlign));
	set_flag(align,*(LOGICAL(AS_LOGICAL(val))),BAM_FPROPER_PAIR);
	return R_NilValue;
}

SEXP bam_align_set_is_unmapped(SEXP pAlign, SEXP val)
{
	if(TYPEOF(pAlign)!=EXTPTRSXP)
	{
		error("[bam_align_set_is_unmapped] No external pointer!");
		return R_NilValue;
	}
	if(TYPEOF(val)!=LGLSXP)
	{
		error("[bam_align_set_is_unmapped] No bool value!");
		return R_NilValue;
	}
	bam1_t *align=(bam1_t*)(R_ExternalPtrAddr(pAlign));
	set_flag(align,*(LOGICAL(AS_LOGICAL(val))),BAM_FUNMAP);
	return R_NilValue;
}

SEXP bam_align_set_mate_is_unmapped(SEXP pAlign, SEXP val)
{
	if(TYPEOF(pAlign)!=EXTPTRSXP)
	{
		error("[bam_align_set_mate_is_unmapped] No external pointer!");
		return R_NilValue;
	}
	if(TYPEOF(val)!=LGLSXP)
	{
		error("[bam_align_set_mate_is_unmapped] No bool value!");
		return R_NilValue;
	}
	bam1_t *align=(bam1_t*)(R_ExternalPtrAddr(pAlign));
	set_flag(align,*(LOGICAL(AS_LOGICAL(val))),BAM_FMUNMAP);
	return R_NilValue;
}

SEXP bam_align_set_strand_reverse(SEXP pAlign, SEXP val)
{
	if(TYPEOF(pAlign)!=EXTPTRSXP)
	{
		error("[bam_align_set_strand_reverse] No external pointer!");
		return R_NilValue;
	}
	if(TYPEOF(val)!=LGLSXP)
	{
		error("[bam_align_set_strand_reverse] No bool value!");
		return R_NilValue;
	}
	bam1_t *align=(bam1_t*)(R_ExternalPtrAddr(pAlign));
	set_flag(align,*(LOGICAL(AS_LOGICAL(val))),16);
	return R_NilValue;
}

SEXP bam_align_set_mate_strand_reverse(SEXP pAlign, SEXP val)
{
	if(TYPEOF(pAlign)!=EXTPTRSXP)
	{
		error("[bam_align_set_mate_strand_reverse] No external pointer!");
		return R_NilValue;
	}
	if(TYPEOF(val)!=LGLSXP)
	{
		error("[bam_align_set_mate_strand_reverse] No bool value!");
		return R_NilValue;
	}
	bam1_t *align=(bam1_t*)(R_ExternalPtrAddr(pAlign));
	set_flag(align,*(LOGICAL(AS_LOGICAL(val))),32);
	return R_NilValue;
}

SEXP bam_align_set_is_first_in_pair(SEXP pAlign, SEXP val)
{
	if(TYPEOF(pAlign)!=EXTPTRSXP)
	{
		error("[bam_align_set_is_first_in_pair] No external pointer!");
		return R_NilValue;
	}
	if(TYPEOF(val)!=LGLSXP)
	{
		error("[bam_align_set_is_first_in_pair] No bool value!");
		return R_NilValue;
	}
	bam1_t *align=(bam1_t*)(R_ExternalPtrAddr(pAlign));
	set_flag(align,*(LOGICAL(AS_LOGICAL(val))),BAM_FREAD1);
	return R_NilValue;
}

SEXP bam_align_set_is_second_in_pair(SEXP pAlign, SEXP val)
{
	if(TYPEOF(pAlign)!=EXTPTRSXP)
	{
		error("[bam_align_set_is_second_in_pair] No external pointer!");
		return R_NilValue;
	}
	if(TYPEOF(val)!=LGLSXP)
	{
		error("[bam_align_set_is_second_in_pair] No bool value!");
		return R_NilValue;
	}
	bam1_t *align=(bam1_t*)(R_ExternalPtrAddr(pAlign));
	set_flag(align,*(LOGICAL(AS_LOGICAL(val))),BAM_FREAD2);
	return R_NilValue;
}

SEXP bam_align_set_is_secondary_align(SEXP pAlign, SEXP val)
{
	if(TYPEOF(pAlign)!=EXTPTRSXP)
	{
		error("[bam_align_set_is_secondary_align] No external pointer!");
		return R_NilValue;
	}
	if(TYPEOF(val)!=LGLSXP)
	{
		error("[bam_align_set_is_secondary_align] No bool value!");
		return R_NilValue;
	}
	bam1_t *align=(bam1_t*)(R_ExternalPtrAddr(pAlign));
	set_flag(align,*(LOGICAL(AS_LOGICAL(val))),BAM_FSECONDARY);
	return R_NilValue;
}

SEXP bam_align_set_fail_qc(SEXP pAlign, SEXP val)
{
	if(TYPEOF(pAlign)!=EXTPTRSXP)
	{
		error("[bam_align_set_fail_qc] No external pointer!");
		return R_NilValue;
	}
	if(TYPEOF(val)!=LGLSXP)
	{
		error("[bam_align_set_fail_qc] No bool value!");
		return R_NilValue;
	}
	bam1_t *align=(bam1_t*)(R_ExternalPtrAddr(pAlign));
	set_flag(align,*(LOGICAL(AS_LOGICAL(val))),BAM_FQCFAIL);
	return R_NilValue;
}

SEXP bam_align_set_is_pcr_or_optical_dup(SEXP pAlign, SEXP val)
{
	if(TYPEOF(pAlign)!=EXTPTRSXP)
	{
		error("[bam_align_set_is_pcr_or_optical_dup] No external pointer!");
		return R_NilValue;
	}
	if(TYPEOF(val)!=LGLSXP)
	{
		error("[bam_align_set_is_pcr_or_optical_dup] No bool value!");
		return R_NilValue;
	}
	bam1_t *align=(bam1_t*)(R_ExternalPtrAddr(pAlign));
	set_flag(align,*(LOGICAL(AS_LOGICAL(val))),BAM_FDUP);
	return R_NilValue;
}

SEXP bam_align_set_flag(SEXP pAlign, SEXP val)
{
	if(TYPEOF(pAlign)!=EXTPTRSXP)
	{
		error("[bam_align_set_flag] No external pointer!");
		return R_NilValue;
	}
	if(TYPEOF(val)!=INTSXP)
	{
		error("[bam_align_set_flag] No integer value!");
		return R_NilValue;
	}
	bam1_t *align=(bam1_t*)(R_ExternalPtrAddr(pAlign));
	align->core.flag=*INTEGER(val);
	return R_NilValue;
}
#endif
