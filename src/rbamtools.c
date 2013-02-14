/*
 *  File      : rbamtools.c
 *
 *  Created on: 17.06.2011
 *  Author    : W. Kaisers
 *  Content   : Function definitions in C for R package rbamtools
 *
 *  Change log: 
 *              29.Okt12  (nemesis) Changed gap_list_get_df cigar_type output to 'factor'
 *              30.Okt.12 (phoibe)  Removed gap_list_fetch message
 *              31.Okt.12 (phoibe)  [bam_reader_save_aligns] function & SAM_TYPE_READ added
 */

#ifndef rbamtools_c
#define rbamtools_c
#include "rbamtools.h"

void print_bitmask_dec(bitmap_type *value)
{
	int i,last_block;
	last_block=bitmap_size-1;
	Rprintf("%3u",getByte(*value,last_block));
	for(i=last_block-1;i>=0;--i)
		Rprintf(" | %3u",getByte(*value,i));
	Rprintf("\n");
}

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
	//else
	//	Rprintf("[bamReader] Opened file \"%s\"\n",_filename);

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
	Rprintf("[bamReader] closed.\n");
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

SEXP bam_reader_save_aligns(SEXP pReader,SEXP pWriter)
{
	if(TYPEOF(pReader)!=EXTPTRSXP)
	{
		error("[bam_reader_save_aligns] pReader: No external pointer!\n");
		return R_NilValue;
	}
	if(TYPEOF(pWriter)!=EXTPTRSXP)
	{
		error("[bam_reader_save_aligns] pWriter: No external pointer!\n");
		return R_NilValue;
	}

	// reader
	samfile_t *reader=(samfile_t*)(R_ExternalPtrAddr(pReader));
	if (reader == 0 || !(reader->type & SAM_TYPE_READ))
		error("[bam_reader_save_aligns] Reader not open for reading!");

	// writer
	samfile_t *writer=(samfile_t*)R_ExternalPtrAddr(pWriter);
	if (writer == 0 || (writer->type & SAM_TYPE_READ))
		error("[bam_reader_save_aligns] Writer not open for writing!");


	bam1_t *align=bam_init1();
	unsigned long nAligns=0;
	int res=samread(reader,align);
	while(res>0)
	{
		samwrite(writer,align);
		++nAligns;
		if(nAligns % 1000000 ==0)
			Rprintf("[bam_reader_save_aligns] %u aligns written.\n",nAligns);
		res=samread(reader,align);
	}

	if(res==-2)
	{
		error("[bam_reader_save_aligns] samread found truncated BAM-file.\n");
		return R_NilValue;
	}

	SEXP ans;
	PROTECT(ans=allocVector(INTSXP,1));
	INTEGER(ans)[0]=nAligns;
	UNPROTECT(1);
	return ans;
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

SEXP bam_reader_tell(SEXP pReader)
{
	if(TYPEOF(pReader)!=EXTPTRSXP)
	{
		error("[bam_reader_tell] No external pointer!\n");
		return R_NilValue;
	}
	return Rf_ScalarReal(bam_tell(((samfile_t*)(R_ExternalPtrAddr(pReader)))->x.bam));
}

SEXP bam_reader_seek(SEXP pReader, SEXP pPos)
{
	if(TYPEOF(pReader)!=EXTPTRSXP)
	{
		error("[bam_reader_seek] No external pointer!\n");
		return R_NilValue;
	}
	if(TYPEOF(pPos)!=REALSXP)
	{
		error("[bam_reader_seek] Position must be numeric!\n");
		return R_NilValue;
	}
	if(LENGTH(pPos)>1)
	{
		error("[bam_reader_seek] Length of position must be 1!\n");
		return R_NilValue;
	}

	samfile_t *reader=(samfile_t*)(R_ExternalPtrAddr(pReader));
	double *pos=REAL(pPos);
    int64_t res=bam_seek(reader->x.bam,(int64_t)*pos,SEEK_SET);
    if(res<0)
    	Rprintf("[bam_reader_seek] bam_seek fails!\n");
	return R_NilValue;
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
		error("[gap_list_fetch] pReader is No external pointer!\n");
		return R_NilValue;
	}
	if(TYPEOF(pIndex)!=EXTPTRSXP)
	{
		error("[gap_list_fetch] pIndex is No external pointer!\n");
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
	//Rprintf("[gap_list_fetch] Fetched list of size %i.\n",l->size);
	return list;
}

SEXP gap_list_get_df(SEXP pGapList)
{
	if(TYPEOF(pGapList)!=EXTPTRSXP)
		error("[get_gap_list_df] No external pointer!");

	gap_list *l=(gap_list*)(R_ExternalPtrAddr(pGapList));

	// create data.frame
	int nProtected=0;
	int nCols=9;
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

	// Column 5: gap_len
	SEXP gap_len_vector;
	PROTECT(gap_len_vector=allocVector(INTSXP,nRows));
	++nProtected;

	// Column 6: right_start
	SEXP right_start_vector;
	PROTECT(right_start_vector=allocVector(INTSXP,nRows));
	++nProtected;

	// Column 7: right_cigar_len
	SEXP right_cigar_len_vector;
	PROTECT(right_cigar_len_vector=allocVector(INTSXP,nRows));
	++nProtected;

	// Column 8: right_cigar_type
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
		INTEGER(gap_len_vector)			[i]=pos_data[5];
		INTEGER(right_start_vector)		[i]=pos_data[6];
		INTEGER(right_cigar_len_vector)	[i]=pos_data[7];
		INTEGER(right_cigar_type_vector)[i]=pos_data[8];
		l->curr_el=l->curr_el->next_el;
	}

	// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
	// Convert left_cigar_type to vector

	SEXP levs;
	int nLevels=9;
	PROTECT(levs=allocVector(STRSXP,nLevels));
	++nProtected;

	SET_STRING_ELT(levs,0,mkChar("M"));
	SET_STRING_ELT(levs,1,mkChar("I"));
	SET_STRING_ELT(levs,2,mkChar("D"));
	SET_STRING_ELT(levs,3,mkChar("N"));
	SET_STRING_ELT(levs,4,mkChar("S"));
	SET_STRING_ELT(levs,5,mkChar("H"));
	SET_STRING_ELT(levs,6,mkChar("P"));
	SET_STRING_ELT(levs,7,mkChar("="));
	SET_STRING_ELT(levs,8,mkChar("X"));
	setAttrib(left_cigar_type_vector,R_LevelsSymbol,levs);

	SEXP csymb;
	PROTECT(csymb=mkString("factor"));
	++nProtected;
	setAttrib(left_cigar_type_vector,R_ClassSymbol,csymb);

	// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
	// Convert right_cigar_type_vector to vector

	nLevels=9;
	PROTECT(levs=allocVector(STRSXP,nLevels));
	++nProtected;

	SET_STRING_ELT(levs,0,mkChar("M"));
	SET_STRING_ELT(levs,1,mkChar("I"));
	SET_STRING_ELT(levs,2,mkChar("D"));
	SET_STRING_ELT(levs,3,mkChar("N"));
	SET_STRING_ELT(levs,4,mkChar("S"));
	SET_STRING_ELT(levs,5,mkChar("H"));
	SET_STRING_ELT(levs,6,mkChar("P"));
	SET_STRING_ELT(levs,7,mkChar("="));
	SET_STRING_ELT(levs,8,mkChar("X"));
	setAttrib(right_cigar_type_vector,R_LevelsSymbol,levs);

	PROTECT(csymb=mkString("factor"));
	++nProtected;
	setAttrib(right_cigar_type_vector,R_ClassSymbol,csymb);
	// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

	SET_VECTOR_ELT(dflist,0,ref_vector);
	SET_VECTOR_ELT(dflist,1,pos_vector);
	SET_VECTOR_ELT(dflist,2,left_cigar_len_vector);
	SET_VECTOR_ELT(dflist,3,left_cigar_type_vector);
	SET_VECTOR_ELT(dflist,4,left_stop_vector);
	SET_VECTOR_ELT(dflist,5,gap_len_vector);
	SET_VECTOR_ELT(dflist,6,right_start_vector);
	SET_VECTOR_ELT(dflist,7,right_cigar_len_vector);
	SET_VECTOR_ELT(dflist,8,right_cigar_type_vector);

	///////////////////////////////////////////////////////////////////////////
	// Column Names
	SEXP col_names;
	PROTECT(col_names=allocVector(STRSXP,nCols));
	++nProtected;

	SET_STRING_ELT(col_names,0,mkChar("refid"));
	SET_STRING_ELT(col_names,1,mkChar("position"));
	SET_STRING_ELT(col_names,2,mkChar("left_cigar_len"));
	SET_STRING_ELT(col_names,3,mkChar("left_cigar_type"));
	SET_STRING_ELT(col_names,4,mkChar("left_stop"));
	SET_STRING_ELT(col_names,5,mkChar("gaplen"));
	SET_STRING_ELT(col_names,6,mkChar("right_start"));
	SET_STRING_ELT(col_names,7,mkChar("right_cigar_len"));
	SET_STRING_ELT(col_names,8,mkChar("right_cigar_type"));
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

SEXP gap_list_get_size(SEXP pGapList)
{
	if(TYPEOF(pGapList)!=EXTPTRSXP)
		error("[gap_list_get_size] No external pointer!");

	gap_list *l=(gap_list*)(R_ExternalPtrAddr(pGapList));
	SEXP ans;
	PROTECT(ans=allocVector(INTSXP,1));
	INTEGER(ans)[0]=(l->size);
	UNPROTECT(1);
	return ans;
}

SEXP gap_list_get_nAligns(SEXP pGapList)
{
	if(TYPEOF(pGapList)!=EXTPTRSXP)
		error("[gap_list_get_size] No external pointer!");

	gap_list *l=(gap_list*)(R_ExternalPtrAddr(pGapList));
	SEXP ans;
	PROTECT(ans=allocVector(INTSXP,1));
	INTEGER(ans)[0]=(l->nAligns);
	UNPROTECT(1);
	return ans;
}
SEXP gap_list_get_nGapAligns(SEXP pGapList)
{
	if(TYPEOF(pGapList)!=EXTPTRSXP)
		error("[gap_list_get_size] No external pointer!");

	gap_list *l=(gap_list*)(R_ExternalPtrAddr(pGapList));
	SEXP ans;
	PROTECT(ans=allocVector(INTSXP,1));
	INTEGER(ans)[0]=(l->nGapAligns);
	UNPROTECT(1);
	return ans;
}


///////////////////////////////////////////////////////////////////////////////////////////////////
// gap_site_list
///////////////////////////////////////////////////////////////////////////////////////////////////

static void finalize_gap_site_list(SEXP ptr)
{
	if(TYPEOF(ptr)!=EXTPTRSXP)
	{
		error("[finalize_gap_list] No external pointer!\n");
		return;
	}
	struct site_list *l=(site_list*)(R_ExternalPtrAddr(ptr));
	site_list_destroy(l);
	l=NULL;
	R_SetExternalPtrAddr(ptr,NULL);
}

SEXP create_gap_site_list()
{
	site_list *l=site_list_init();
	SEXP list;
	PROTECT(list=R_MakeExternalPtr((void*)(l),R_NilValue,R_NilValue));
	R_RegisterCFinalizer(list,finalize_gap_site_list);
	UNPROTECT(1);
	return list;
}

static int gap_site_list_fetch_func(const bam1_t *b, void *data)
{
	site_list *l=(site_list*)data;
	list_gap_sites(l,b);
	return 0;
}

SEXP gap_site_list_fetch(SEXP pReader,SEXP pIndex,SEXP pCoords)
{
	if(TYPEOF(pReader)!=EXTPTRSXP)
	{
		error("[gap_site_list_fetch] pReader is No external pointer!\n");
		return R_NilValue;
	}
	if(TYPEOF(pIndex)!=EXTPTRSXP)
	{
		error("[gap_site_list_fetch] pIndex is No external pointer!\n");
		return R_NilValue;
	}
	if(TYPEOF(pCoords)!=REALSXP)
	{
		error("[gap_site_list_fetch] pCoords is no REAL!\n");
		return R_NilValue;
	}
	if(LENGTH(pCoords)!=3)
	{
		error("[gap_site_list_fetch] pCoords must contain three values (refid,begin,end)!\n");
		return R_NilValue;
	}

	samfile_t *reader=(samfile_t*)(R_ExternalPtrAddr(pReader));
	bam_index_t *index=(bam_index_t*)(R_ExternalPtrAddr(pIndex));
	if(reader==NULL)
	{
		error("[gap_site_list_fetch] Reader must not be NULL pointer!\n");
		return R_NilValue;
	}
	if(index==NULL)
	{
		error("[gap_site_list_fetch] Index must not be NULL pointer!\n");
		return R_NilValue;
	}

	double *pi=REAL(pCoords);
	int refid=(int) pi[0];
	int begin=(int) pi[1];
	int end=(int) pi[2];

	if(refid<0 || refid >=(reader->header->n_targets))
	{
		error("[gap_site_list_fetch] refid out of range!\n");
		return R_NilValue;
	}
	if(begin<0 || begin>=end || end>(reader->header->target_len[refid]))
	{
		error("[gap_site_list_fetch] Begin or end out of range!\n");
		return R_NilValue;
	}

	site_list *l=site_list_init();
	l->refid=refid;
	bam_fetch(reader->x.bam,index,refid,begin,end,(void*)l,gap_site_list_fetch_func);

	SEXP list;
	PROTECT(list=R_MakeExternalPtr( (void*)(l),R_NilValue,R_NilValue));
	R_RegisterCFinalizer(list,finalize_gap_site_list);
	UNPROTECT(1);
	Rprintf("[gap_site_list_fetch] Fetched list of size %i.\n",l->size);
	return list;
}

SEXP gap_site_list_get_df(SEXP pGapList)
{
	if(TYPEOF(pGapList)!=EXTPTRSXP)
		error("[get_gap_list_df] No external pointer!");

	site_list *l=(site_list*)(R_ExternalPtrAddr(pGapList));

	// create data.frame
	int nProtected=0;
	int nCols=13;
	SEXP dflist;
	PROTECT(dflist=allocVector(VECSXP,nCols));
	++nProtected;
	int nRows=(l->size);
	int i;

	// Column 0: id
	SEXP id_vector;
	PROTECT(id_vector=allocVector(INTSXP,nRows));
	++nProtected;
	// Column 1: refid
	SEXP refid_vector;
	PROTECT(refid_vector=allocVector(INTSXP,nRows));
	++nProtected;
	// Column 2: lstart
	SEXP lstart_vector;
	PROTECT(lstart_vector=allocVector(INTSXP,nRows));
	++nProtected;
	// Column 3: lend
	SEXP lend_vector;
	PROTECT(lend_vector=allocVector(INTSXP,nRows));
	++nProtected;
	// Column 4: rstart
	SEXP rstart_vector;
	PROTECT(rstart_vector=allocVector(INTSXP,nRows));
	++nProtected;
	// Column 5: rend
	SEXP rend_vector;
	PROTECT(rend_vector=allocVector(INTSXP,nRows));
	++nProtected;
	// Column 6: gap_len
	SEXP gap_len_vector;
	PROTECT(gap_len_vector=allocVector(INTSXP,nRows));
	++nProtected;
	// Column 7: nAligns
	SEXP nAligns_vector;
	PROTECT(nAligns_vector=allocVector(INTSXP,nRows));
	++nProtected;
	// Column 8: nProbes
	SEXP nProbes_vector;
	PROTECT(nProbes_vector=allocVector(INTSXP,nRows));
	++nProtected;
	// Column 9: nlstart
	SEXP nlstart_vector;
	PROTECT(nlstart_vector=allocVector(INTSXP,nRows));
	++nProtected;
	// Column 10: sum_left_cigar
	SEXP qmm;
	PROTECT(qmm=allocVector(INTSXP,nRows));
	++nProtected;
	// Column 11: lcs
	SEXP nmcs_vector;
	PROTECT(nmcs_vector=allocVector(INTSXP,nRows));
	++nProtected;
	// Column 12: mcs
	SEXP mcs_vector;
	PROTECT(mcs_vector=allocVector(INTSXP,nRows));
	++nProtected;


	// start position
	site_list_element *el,*curr;
	curr=l->curr;
	site_list_set_curr_first(l);

	for(i=0;i<nRows;++i)
	{
		el=site_list_get_curr_pp(l);
		INTEGER(id_vector)                  [i]=i+1;
		INTEGER(refid_vector)        		[i]=l->refid;
		// lend - max(left_cigar_len)+1
		INTEGER(lstart_vector)        		[i]=el->lend-getByte(el->lcs,0)+1;
		INTEGER(lend_vector)	            [i]=el->lend;
		INTEGER(rstart_vector)	            [i]=el->rstart;
		INTEGER(rend_vector)  		        [i]=el->rstart+el->r_cigar_size-1;
		INTEGER(gap_len_vector)			    [i]=el->gap_len;
		INTEGER(nAligns_vector)		        [i]=el->nAligns;
		INTEGER(nProbes_vector)             [i]=el->nProbes;
		INTEGER(nlstart_vector)	            [i]=bitmask_nPos(el->lcs);
		INTEGER(qmm)      [i]=bitmask_sumPos(el->lcs);
		INTEGER(nmcs_vector)					[i]=el->lcs;
		INTEGER(mcs_vector)					[i]=el->mcs;
	}
	// Reset curr position
	l->curr=curr;

	SET_VECTOR_ELT(dflist,0,id_vector);
	SET_VECTOR_ELT(dflist,1,refid_vector);
	SET_VECTOR_ELT(dflist,2,lstart_vector);
	SET_VECTOR_ELT(dflist,3,lend_vector);
	SET_VECTOR_ELT(dflist,4,rstart_vector);
	SET_VECTOR_ELT(dflist,5,rend_vector);
	SET_VECTOR_ELT(dflist,6,gap_len_vector);
	SET_VECTOR_ELT(dflist,7,nAligns_vector);
	SET_VECTOR_ELT(dflist,8,nProbes_vector);
	SET_VECTOR_ELT(dflist,9,nlstart_vector);
	SET_VECTOR_ELT(dflist,10,qmm);
	SET_VECTOR_ELT(dflist,11,nmcs_vector);
	SET_VECTOR_ELT(dflist,12,mcs_vector);

	///////////////////////////////////////////////////////////////////////////
	// Column Names
	SEXP col_names;
	PROTECT(col_names=allocVector(STRSXP,nCols));
	++nProtected;


	SET_STRING_ELT(col_names,0,mkChar("id"));
	SET_STRING_ELT(col_names,1,mkChar("refid"));
	SET_STRING_ELT(col_names,2,mkChar("lstart"));
	SET_STRING_ELT(col_names,3,mkChar("lend"));
	SET_STRING_ELT(col_names,4,mkChar("rstart"));
	SET_STRING_ELT(col_names,5,mkChar("rend"));
	SET_STRING_ELT(col_names,6,mkChar("gaplen"));
	SET_STRING_ELT(col_names,7,mkChar("nAligns"));
	SET_STRING_ELT(col_names,8,mkChar("nProbes"));
	SET_STRING_ELT(col_names,9,mkChar("nlstart"));
	SET_STRING_ELT(col_names,10,mkChar("lm_sum"));
	SET_STRING_ELT(col_names,11,mkChar("lcs"));
	SET_STRING_ELT(col_names,12,mkChar("mcs"));
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

SEXP gap_site_list_get_ref_id(SEXP pGapList)
{
	if(TYPEOF(pGapList)!=EXTPTRSXP)
		error("[gap_site_list_get_ref_id] No external pointer!");

	site_list *l=(site_list*)(R_ExternalPtrAddr(pGapList));
	SEXP ans;
	PROTECT(ans=allocVector(INTSXP,1));
	INTEGER(ans)[0]=(l->refid);
	UNPROTECT(1);
	return ans;
}

SEXP gap_site_list_get_size(SEXP pGapList)
{
	if(TYPEOF(pGapList)!=EXTPTRSXP)
		error("[gap_site_list_get_size] No external pointer!");

	site_list *l=(site_list*)(R_ExternalPtrAddr(pGapList));
	SEXP ans;
	PROTECT(ans=allocVector(INTSXP,1));
	INTEGER(ans)[0]=(l->size);
	UNPROTECT(1);
	return ans;
}

SEXP gap_site_list_get_nAligns(SEXP pGapList)
{
	if(TYPEOF(pGapList)!=EXTPTRSXP)
		error("[gap_site_list_get_nAligns] No external pointer!");

	site_list *l=(site_list*)(R_ExternalPtrAddr(pGapList));
	SEXP ans;
	PROTECT(ans=allocVector(INTSXP,1));
	INTEGER(ans)[0]=(l->nAligns);
	UNPROTECT(1);
	return ans;
}

SEXP gap_site_list_get_nGapAligns(SEXP pGapList)
{
	if(TYPEOF(pGapList)!=EXTPTRSXP)
		error("[gap_site_list_get_nGapAligns] No external pointer!");

	site_list *l=(site_list*)(R_ExternalPtrAddr(pGapList));
	SEXP ans;
	PROTECT(ans=allocVector(INTSXP,1));
	INTEGER(ans)[0]=(l->nGapAligns);
	UNPROTECT(1);
	return ans;
}

SEXP gap_site_list_merge(SEXP pLhs, SEXP pRhs, SEXP pRef)
{
	if(TYPEOF(pLhs)!=EXTPTRSXP)
		error("[gap_site_list_merge] pLhs: No external pointer");
	if(TYPEOF(pRhs)!=EXTPTRSXP)
		error("[gap_site_list_merge] pRhs: No external pointer");
	if(TYPEOF(pRef)!=INTSXP)
		error("[gap_site_list_merge] pRef must be Integer!");

	site_list *lhs=(site_list*)(R_ExternalPtrAddr(pLhs));
	site_list *rhs=(site_list*)(R_ExternalPtrAddr(pRhs));
	site_list *mrg=site_list_merge(lhs,rhs,INTEGER(pRef)[0]);

	SEXP list;
	PROTECT(list=R_MakeExternalPtr( (void*)(mrg),R_NilValue,R_NilValue));
	R_RegisterCFinalizer(list,finalize_gap_site_list);
	UNPROTECT(1);
	Rprintf("[gap_site_list_merge] Merge list of size %i.\n",mrg->size);
	return list;
}

SEXP gap_site_list_copy(SEXP pGapList)
{
	if(TYPEOF(pGapList)!=EXTPTRSXP)
		error("[gap_site_list_copy] No external pointer!");
	site_list *src=(site_list*)(R_ExternalPtrAddr(pGapList));
	site_list *tar=site_list_copy(src);

	SEXP list;
	PROTECT(list=R_MakeExternalPtr((void*)(tar),R_NilValue,R_NilValue));
	R_RegisterCFinalizer(list,finalize_gap_site_list);
	UNPROTECT(1);
	return list;

}

SEXP bitmask_r_zip(SEXP lhs, SEXP rhs)
{
	if(TYPEOF(lhs)!=INTSXP)
		error("[bitmask_r_zip] lhs must be integer!");
	if(TYPEOF(rhs)!=INTSXP)
		error("[bitmask_r_zip] rhs must be integer!");
	if(LENGTH(lhs)!=LENGTH(rhs))
		error("[bitmask_r_zip] lhs and rhs must have same length!");

	const int len=LENGTH(lhs);
	SEXP res;
	PROTECT(res=allocVector(INTSXP,len));
	int i;
	for(i=0;i<len;++i)
		INTEGER(res)[i]=r_zip_val(INTEGER(lhs)[i],INTEGER(rhs)[i]);

	UNPROTECT(1);
	return(res);
}

///////////////////////////////////////////////////////////////////////////////////////////////////
// gap_site_ll
///////////////////////////////////////////////////////////////////////////////////////////////////

static void finialize_gap_site_ll(SEXP ptr)
{
	if(TYPEOF(ptr)!=EXTPTRSXP)
		error("[finalize_gap_site_ll] No external pointer!\n");
	site_ll *l=(site_ll*) (R_ExternalPtrAddr(ptr));
	site_ll_destroy(l);
	l=NULL;
	R_SetExternalPtrAddr(ptr,NULL);
}

SEXP gap_site_ll_init()
{
	site_ll *l=site_ll_init();
	SEXP list;
	PROTECT(list=R_MakeExternalPtr( (void*)(l),R_NilValue,R_NilValue));
	R_RegisterCFinalizer(list,finialize_gap_site_ll);
	UNPROTECT(1);
	return list;
}

SEXP gap_site_ll_fetch(SEXP pReader, SEXP pIndex, SEXP pRefid, SEXP pStart, SEXP pEnd)
{
	if(TYPEOF(pReader)!=EXTPTRSXP)
		error("[gap_site_ll_fetch] pReader is no external pointer!\n");
	if(TYPEOF(pIndex) !=EXTPTRSXP)
		error("[gap_site_ll_fetch] pIndex is no external pointer!\n");
	if(TYPEOF(pRefid) !=INTSXP)
		error("[gap_site_ll_fetch] pRefid is no INT!\n");
	if(TYPEOF(pStart) !=INTSXP)
		error("[gap_site_ll_fetch] pStart is no INT!\n");
	if(TYPEOF(pEnd)   !=INTSXP)
		error("[gap_site_ll_fetch] pEnd is no INT!\n");
	if(LENGTH(pRefid)!=LENGTH(pStart) || LENGTH(pStart)!=LENGTH(pEnd))
		error("[gap_site_ll_fetch] Unequal Length of pRefid, pStart and pEnd!\n");

	samfile_t *reader=(samfile_t*)(R_ExternalPtrAddr(pReader));
	bam_index_t *index=(bam_index_t*)(R_ExternalPtrAddr(pIndex));
	if(reader==NULL)
	{
		error("[gap_site_ll_fetch] Reader must not be NULL pointer!\n");
		return R_NilValue;
	}
	if(index==NULL)
	{
		error("[gap_site_ll_fetch] Index must not be NULL pointer!\n");
		return R_NilValue;
	}

	int i,refid,begin,end;
	site_list *sl;

	int len=LENGTH(pStart);
	site_ll *l=site_ll_init();

	for(i=0;i<len;++i,++sl)
	{
		sl=site_list_init();
		refid=INTEGER(pRefid)[i];
		begin=INTEGER(pStart)[i];
		end  =INTEGER(pEnd)  [i];
		sl->refid=refid;

		if(refid<0 || refid>=(reader->header->n_targets))
			error("[gap_site_ll_fetch] refid out of range!");
		if(begin<0 || begin>=end || end>(reader->header->target_len[refid]))
			error("[gap_site_ll_fetch] Begin or end out of range!");
		bam_fetch(reader->x.bam,index,refid,begin,end,(void*)sl,gap_site_list_fetch_func);
		if(sl->size>0)
			site_ll_add_site_list(l,sl);
	}

	SEXP list;
	PROTECT(list=R_MakeExternalPtr( (void*)(l),R_NilValue,R_NilValue));
	R_RegisterCFinalizer(list,finialize_gap_site_ll);
	UNPROTECT(1);
	return list;
}

SEXP gap_site_ll_set_curr_first(SEXP pGapList)
{
	if(TYPEOF(pGapList)!=EXTPTRSXP)
		error("[gap_site_ll_set_curr_first] pGapList must be external pointer!");
	site_ll *sll=(site_ll*)(R_ExternalPtrAddr(pGapList));
	site_ll_set_curr_first(sll);
	return R_NilValue;
}

SEXP gap_site_ll_add_curr_pp(SEXP pSrc,SEXP pTrg,SEXP pRefid)
{
	if(TYPEOF(pSrc)!=EXTPTRSXP)
		error("[gap_site_ll_add_curr_pp] pSrc must be external pointer!");
	if(TYPEOF(pTrg)!=EXTPTRSXP)
		error("[gap_site_ll_add_curr_pp] pTrg must be external pointer!");
	if(TYPEOF(pRefid)!=INTSXP)
		error("[gap_site_ll_add_curr_pp] pRefid must be Integer!");

	site_ll *src=(site_ll*)(R_ExternalPtrAddr(pSrc));
	site_ll *trg=(site_ll*)(R_ExternalPtrAddr(pTrg));
	site_list *l=site_ll_get_curr_site_list_pp(src);

	site_list *ins=site_list_copy(l);
	ins->refid=INTEGER(pRefid)[0];
	site_ll_add_site_list(trg,ins);

	return R_NilValue;
}

SEXP gap_site_ll_add_merge_pp(SEXP plSrc,SEXP prSrc,SEXP pTrg,SEXP pRefid)
{
	if(TYPEOF(plSrc)!=EXTPTRSXP)
		error("[gap_site_ll_add_merge_pp] plSrc must be external pointer!");
	if(TYPEOF(prSrc)!=EXTPTRSXP)
		error("[gap_site_ll_add_merge_pp] prSrc must be external pointer!");
	if(TYPEOF(pTrg)!=EXTPTRSXP)
		error("[gap_site_ll_add_merge_pp] pTrg must be external pointer!");
	if(TYPEOF(pRefid)!=INTSXP)
		error("[gap_site_ll_add_merge_pp] pRefid must be Integer!");

	unsigned refid=INTEGER(pRefid)[0];

	site_ll *l_src=(site_ll*)(R_ExternalPtrAddr(plSrc));
	site_ll *r_src=(site_ll*)(R_ExternalPtrAddr(prSrc));
	site_ll *trg  =(site_ll*)(R_ExternalPtrAddr(pTrg));

	//Rprintf("[gap_site_ll_add_merge_pp] plSrc size: %u\tprSrc size: %u\tpTrg size: %u\trefid: %u\n",l_src->size,r_src->size,trg->size,refid);

	// Get pointer to current site_list
	site_list *ls=site_ll_get_curr_site_list_pp(l_src);
	if(ls==0)
		error("[gap_site_ll_add_merge_pp] l_src curr returned 0! Use 'set_curr_first' or check size>0!");
	//Rprintf("[gap_site_ll_add_merge_pp] l_src size: %u\n",ls->size);
	site_list *rs=site_ll_get_curr_site_list_pp(r_src);
	if(rs==0)
		error("[gap_site_ll_add_merge_pp] r_src curr returned 0! Use 'set_curr_first' or check size>0!");
	//Rprintf("[gap_site_ll_add_merge_pp] r_src size: %u\n",rs->size);

	// Do merging and add to target list
	site_list *l=site_list_merge(ls,rs,refid);
	l->nAligns=ls->nAligns+rs->nAligns;
	l->nGapAligns=ls->nGapAligns+rs->nGapAligns;

	site_ll_add_site_list(trg,l);
	//Rprintf("[gap_site_ll_add_merge_pp] Adding site_list of size: %u\trefid: %u\n",l->size,l->refid);
	return R_NilValue;
}


SEXP gap_site_ll_reset_refid(SEXP pGapList)
{
	if(TYPEOF(pGapList)!=EXTPTRSXP)
		error("[gap_site_ll_reset_refid] No external pointer!");
	site_ll *sll=(site_ll*)(R_ExternalPtrAddr(pGapList));

	unsigned size=sll->size;
	SEXP res;
	PROTECT(res=allocVector(INTSXP,size));

	site_list *l;
	site_ll_set_curr_first(sll);
	unsigned i;
	for(i=0;i<size;++i)
	{
		l=site_ll_get_curr_site_list_pp(sll);
		l->refid=i;
		INTEGER(res)[i]=i;
	}
	UNPROTECT(1);
	return(res);
}


SEXP gap_site_ll_get_summary_df(SEXP pGapList)
{
	if(TYPEOF(pGapList)!=EXTPTRSXP)
		error("[gap_site_ll_get_summary_df] No external pointer!");
	site_ll *sll=(site_ll*)(R_ExternalPtrAddr(pGapList));

	// create data.frame
	unsigned nProtected=0;
	unsigned nCols=4;
	unsigned nRows=sll->size;

	SEXP dflist;
	PROTECT(dflist=allocVector(VECSXP,nCols));
	++nProtected;

	// Column 0: refid
	SEXP refid_vector;
	PROTECT(refid_vector=allocVector(INTSXP,nRows));
	++nProtected;

	// Column 1: size
	SEXP size_vector;
	PROTECT(size_vector=allocVector(INTSXP,nRows));
	++nProtected;

	// Column 2: nAligns
	SEXP nalign_vector;
	PROTECT(nalign_vector=allocVector(INTSXP,nRows));
	++nProtected;

	// Column 3: nGapAligns
	SEXP ngap_vector;
	PROTECT(ngap_vector=allocVector(INTSXP,nRows));
	++nProtected;

	site_list *l;
	site_ll_set_curr_first(sll);
	unsigned i;
	for(i=0;i<nRows;++i)
	{
		l=site_ll_get_curr_site_list_pp(sll);
		INTEGER(refid_vector) [i]=l->refid;
		INTEGER(size_vector)  [i]=l->size;
		INTEGER(nalign_vector)[i]=l->nAligns;
		INTEGER(ngap_vector)  [i]=l->nGapAligns;
	}

	SET_VECTOR_ELT(dflist, 0,refid_vector);
	SET_VECTOR_ELT(dflist, 1,size_vector);
	SET_VECTOR_ELT(dflist, 2,nalign_vector);
	SET_VECTOR_ELT(dflist, 3,ngap_vector);

	///////////////////////////////////////////////////////////////////////////
	// Column Names
	SEXP col_names;
	PROTECT(col_names=allocVector(STRSXP,nCols));
	++nProtected;

	SET_STRING_ELT(col_names, 0,mkChar("ID"));
	SET_STRING_ELT(col_names, 1,mkChar("size"));
	SET_STRING_ELT(col_names, 2,mkChar("nAligns"));
	SET_STRING_ELT(col_names, 3,mkChar("nGapAligns"));
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

SEXP gap_site_ll_get_df(SEXP pGapList,SEXP pRefNames)
{
	if(TYPEOF(pGapList)!=EXTPTRSXP)
		error("[gap_site_ll_get_df] No external pointer!");
	if(TYPEOF(pRefNames)!=STRSXP)
		error("[gap_site_ll_get_df] pRefNames must be character!");

	site_ll *sll=(site_ll*)(R_ExternalPtrAddr(pGapList));

	// create data.frame
	unsigned nProtected=0;
	unsigned nCols=13;
	unsigned nRefs=sll->size;
	unsigned nRows=sum_ll_sizes(sll);

	SEXP dflist;
	PROTECT(dflist=allocVector(VECSXP,nCols));
	++nProtected;

	// Column 0: id
	SEXP id_vector;
	PROTECT(id_vector=allocVector(INTSXP,nRows));
	++nProtected;
	// Column 1: refid
	SEXP refid_vector;
	PROTECT(refid_vector=allocVector(INTSXP,nRows));
	++nProtected;
	// Column 2: lstart
	SEXP lstart_vector;
	PROTECT(lstart_vector=allocVector(INTSXP,nRows));
	++nProtected;
	// Column 3: lend
	SEXP lend_vector;
	PROTECT(lend_vector=allocVector(INTSXP,nRows));
	++nProtected;
	// Column 4: rstart
	SEXP rstart_vector;
	PROTECT(rstart_vector=allocVector(INTSXP,nRows));
	++nProtected;
	// Column 5: rend
	SEXP rend_vector;
	PROTECT(rend_vector=allocVector(INTSXP,nRows));
	++nProtected;
	// Column 6: gap_len
	SEXP gap_len_vector;
	PROTECT(gap_len_vector=allocVector(INTSXP,nRows));
	++nProtected;
	// Column 7: nAligns
	SEXP nAligns_vector;
	PROTECT(nAligns_vector=allocVector(INTSXP,nRows));
	++nProtected;
	// Column 8: nAligns
	SEXP nProbes_vector;
	PROTECT(nProbes_vector=allocVector(INTSXP,nRows));
	++nProtected;
	// Column 9: nlstart: number of different lstart positions on gap_site
	SEXP nlstart_vector;
	PROTECT(nlstart_vector=allocVector(INTSXP,nRows));
	++nProtected;

	// Column 10: qmm=Quadrupel mcs mean; mcs contains minimum cigar size values
	// (of flanking left and right M=match segment)
	// Mean of Quadrupel of rightmost values (=4 largest values) in mcs
	SEXP qmm_vector;
	PROTECT(qmm_vector=allocVector(INTSXP,nRows));
	++nProtected;

	// Column 11: nMcs
	SEXP nmcs_vector;
	PROTECT(nmcs_vector=allocVector(INTSXP,nRows));
	++nProtected;

	// Column 12: gqs
	SEXP gqs_vector;
	PROTECT(gqs_vector=allocVector(INTSXP,nRows));
	++nProtected;


	//site_ll_element *ell;
	site_list *l;
	site_list_element *el;

	site_ll_set_curr_first(sll);
	unsigned i,j,k=0;
	for(i=0;i<nRefs;++i)
	{
		l=site_ll_get_curr_site_list_pp(sll);
		site_list_set_curr_first(l);
		for(j=0;j<l->size;++j,++k)
		{
			el=site_list_get_curr_pp(l);
			INTEGER(id_vector)                  [k]=k+1; // row-id: should start with 1
			INTEGER(refid_vector)        		[k]=l->refid+1;
			INTEGER(lstart_vector)        		[k]=el->lend-getByte(el->lcs,0)+1;
			INTEGER(lend_vector)	            [k]=el->lend;
			INTEGER(rstart_vector)	            [k]=el->rstart;
			INTEGER(rend_vector)  		        [k]=el->rstart+el->r_cigar_size-1;
			INTEGER(gap_len_vector)			    [k]=el->gap_len;
			INTEGER(nAligns_vector)		        [k]=el->nAligns;
			INTEGER(nProbes_vector)             [k]=el->nProbes;
			INTEGER(nlstart_vector)	            [k]=bitmask_nPos(el->lcs);
			INTEGER(qmm_vector)					[k]=getQmean(el->mcs);
			INTEGER(nmcs_vector)				[k]=bitmask_nPos(el->mcs);
			INTEGER(gqs_vector)                 [k]=INTEGER(nlstart_vector)[k]*INTEGER(qmm_vector)[k];
		}
	}

	// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
	// Convert left_cigar_type to vector

	SEXP levs;
	int nLevels=LENGTH(pRefNames);
	PROTECT(levs=allocVector(STRSXP,nLevels));
	++nProtected;

	for(i=0;i<nLevels;++i)
		SET_STRING_ELT(levs,i,STRING_ELT(pRefNames,i));
	setAttrib(refid_vector,R_LevelsSymbol,levs);

	SEXP csymb;
	PROTECT(csymb=mkString("factor"));
	++nProtected;
	setAttrib(refid_vector,R_ClassSymbol,csymb);

	//UNPROTECT(nProtected);
	//return(refid_vector);

	SET_VECTOR_ELT(dflist,0,id_vector);
	SET_VECTOR_ELT(dflist,1,refid_vector);
	SET_VECTOR_ELT(dflist,2,lstart_vector);
	SET_VECTOR_ELT(dflist,3,lend_vector);
	SET_VECTOR_ELT(dflist,4,rstart_vector);
	SET_VECTOR_ELT(dflist,5,rend_vector);
	SET_VECTOR_ELT(dflist,6,gap_len_vector);
	SET_VECTOR_ELT(dflist,7,nAligns_vector);
	SET_VECTOR_ELT(dflist,8,nProbes_vector);
	SET_VECTOR_ELT(dflist,9,nlstart_vector);
	SET_VECTOR_ELT(dflist,10,qmm_vector);
	SET_VECTOR_ELT(dflist,11,nmcs_vector);
	SET_VECTOR_ELT(dflist,12,gqs_vector);

	///////////////////////////////////////////////////////////////////////////
	// Column Names
	SEXP col_names;
	PROTECT(col_names=allocVector(STRSXP,nCols));
	++nProtected;

	SET_STRING_ELT(col_names, 0,mkChar("id"));
	SET_STRING_ELT(col_names, 1,mkChar("seqid")); // Needed for spliceSites
	SET_STRING_ELT(col_names, 2,mkChar("lstart"));
	SET_STRING_ELT(col_names, 3,mkChar("lend"));
	SET_STRING_ELT(col_names, 4,mkChar("rstart"));
	SET_STRING_ELT(col_names, 5,mkChar("rend"));
	SET_STRING_ELT(col_names, 6,mkChar("gaplen"));
	SET_STRING_ELT(col_names, 7,mkChar("nAligns"));
	SET_STRING_ELT(col_names, 8,mkChar("nProbes"));
	SET_STRING_ELT(col_names, 9,mkChar("nlstart"));
	SET_STRING_ELT(col_names,10,mkChar("qmm"));
	SET_STRING_ELT(col_names,11,mkChar("nMcs"));
	SET_STRING_ELT(col_names,12,mkChar("gqs"));
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

SEXP gap_site_ll_get_size(SEXP pGapList)
{
	if(TYPEOF(pGapList)!=EXTPTRSXP)
		error("[gap_site_ll_get_size] No external pointer!");
	site_ll *l=(site_ll*)(R_ExternalPtrAddr(pGapList));
	SEXP ans;
	PROTECT(ans=allocVector(REALSXP,1));
	REAL(ans)[0]=(double) sum_ll_sizes(l);
	UNPROTECT(1);
	return ans;
}

SEXP gap_site_ll_get_nAligns(SEXP pGapList)
{
	if(TYPEOF(pGapList)!=EXTPTRSXP)
		error("[gap_site_ll_get_nAligns] No external pointer!");
	site_ll *l=(site_ll*)(R_ExternalPtrAddr(pGapList));
	SEXP ans;
	PROTECT(ans=allocVector(REALSXP,1));
	REAL(ans)[0]=(double)get_nAligns(l);
	UNPROTECT(1);
	return ans;
}

SEXP gap_site_ll_get_nGapAligns(SEXP pGapList)
{
	if(TYPEOF(pGapList)!=EXTPTRSXP)
		error("gap_site_ll_get_nGapAligns] No external pointer!");
	site_ll *l=(site_ll*) (R_ExternalPtrAddr(pGapList));
	SEXP ans;
	PROTECT(ans=allocVector(REALSXP,1));
	REAL(ans)[0]=(double)get_nGapAligns(l);
	UNPROTECT(1);
	return ans;
}


/*



SEXP gap_site_list_get_size(SEXP pGapList)
{
	if(TYPEOF(pGapList)!=EXTPTRSXP)
		error("[gap_site_list_get_size] No external pointer!");

	site_list *l=(site_list*)(R_ExternalPtrAddr(pGapList));
	SEXP ans;
	PROTECT(ans=allocVector(INTSXP,1));
	INTEGER(ans)[0]=(l->size);
	UNPROTECT(1);
	return ans;
}

SEXP gap_site_list_get_nAligns(SEXP pGapList)
{
	if(TYPEOF(pGapList)!=EXTPTRSXP)
		error("[gap_site_list_get_nAligns] No external pointer!");

	site_list *l=(site_list*)(R_ExternalPtrAddr(pGapList));
	SEXP ans;
	PROTECT(ans=allocVector(INTSXP,1));
	INTEGER(ans)[0]=(l->nAligns);
	UNPROTECT(1);
	return ans;
}

SEXP gap_site_list_get_nGapAligns(SEXP pGapList)
{
	if(TYPEOF(pGapList)!=EXTPTRSXP)
		error("[gap_site_list_get_nGapAligns] No external pointer!");

	site_list *l=(site_list*)(R_ExternalPtrAddr(pGapList));
	SEXP ans;
	PROTECT(ans=allocVector(INTSXP,1));
	INTEGER(ans)[0]=(l->nGapAligns);
	UNPROTECT(1);
	return ans;
}
*/

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
		error("[bam_range_fetch] pReader is No external pointer!\n");
		return R_NilValue;
	}
	if(TYPEOF(pIndex)!=EXTPTRSXP)
	{
		error("[bam_range_fetch] pIndex is No external pointer!\n");
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
		error("[bam_range_push_back] pRange is No external pointer!\n");
		return R_NilValue;
	}
	if(TYPEOF(pAlign)!=EXTPTRSXP)
	{
		error("[bam_range_push_back] pAlign is No external pointer!\n");
		return R_NilValue;
	}
	write_current_align((align_list*)(R_ExternalPtrAddr(pRange)),(bam1_t*) (R_ExternalPtrAddr(pAlign)));
	return R_NilValue;
}

SEXP bam_range_insert_past_curr_align(SEXP pRange,SEXP pAlign)
{
	if(TYPEOF(pRange)!=EXTPTRSXP)
	{
		error("[bam_range_push_back] pRange is No external pointer!\n");
		return R_NilValue;
	}
	if(TYPEOF(pAlign)!=EXTPTRSXP)
	{
		error("[bam_range_push_back] pAlign is No external pointer!\n");
		return R_NilValue;
	}
	insert_past_curr_align((align_list*)(R_ExternalPtrAddr(pRange)),(bam1_t*) (R_ExternalPtrAddr(pAlign)));
	return R_NilValue;
}

SEXP bam_range_insert_pre_curr_align(SEXP pRange,SEXP pAlign)
{
	if(TYPEOF(pRange)!=EXTPTRSXP)
	{
		error("[bam_range_push_back] pRange is No external pointer!\n");
		return R_NilValue;
	}
	if(TYPEOF(pAlign)!=EXTPTRSXP)
	{
		error("[bam_range_push_back] pAlign is No external pointer!\n");
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
		error("[bam_range_mv_curr_align] pSrc is No external pointer!\n");
		return R_NilValue;
	}
	if(TYPEOF(pTarget)!=EXTPTRSXP)
	{
		error("[bam_range_mv_curr_align] pTarget is No external pointer!\n");
		return R_NilValue;
	}

	align_list_mv_curr_elem((align_list*)(R_ExternalPtrAddr(pSrc)),(align_list*)(R_ExternalPtrAddr(pTarget)));
	return R_NilValue;
}

SEXP bam_range_get_align_df(SEXP pRange)
{
	if(TYPEOF(pRange)!=EXTPTRSXP)
		error("[bam_range_get_align_df] No external pointer!");

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
		error("[bam_range_write] pRange is No external pointer!\n");
		return R_NilValue;
	}
	if(TYPEOF(pWriter)!=EXTPTRSXP)
	{
		error("[bam_range_write] pWriter is No external pointer!\n");
		return R_NilValue;
	}

	unsigned long bytes_written=0;
	unsigned long range_size, i;

	samfile_t *writer=(samfile_t*) R_ExternalPtrAddr(pWriter);
	align_list *l=(align_list*)(R_ExternalPtrAddr(pRange));

	range_size=l->size;
	// For restoration
	align_element *e=l->curr_el;
	wind_back(l);

	// Retrieve const pointer (i.e. no copy)
	for(i=0;i<range_size;++i)
		bytes_written+=samwrite(writer,get_const_next_align(l));

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
		error("[bam_range_wind_back] pRange is No external pointer!\n");
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
		error("[bam_range_get_size] pRange is No external pointer!\n");
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
		error("[bam_range_push_back] pRange is No external pointer!\n");
		return R_NilValue;
	}
	if(TYPEOF(pAlign)!=EXTPTRSXP)
	{
		error("[bam_range_push_back] pAlign is No external pointer!\n");
		return R_NilValue;
	}
	align_list_push_back((align_list*)(R_ExternalPtrAddr(pRange)),(bam1_t*) (R_ExternalPtrAddr(pAlign)));
	return R_NilValue;
}

SEXP bam_range_push_front(SEXP pRange,SEXP pAlign)
{
	if(TYPEOF(pRange)!=EXTPTRSXP)
	{
		error("[bam_range_push_front] pRange is No external pointer!\n");
		return R_NilValue;
	}
	if(TYPEOF(pAlign)!=EXTPTRSXP)
	{
		error("[bam_range_push_front] pAlign is No external pointer!\n");
		return R_NilValue;
	}
	align_list_push_front((align_list*)(R_ExternalPtrAddr(pRange)),(bam1_t*) (R_ExternalPtrAddr(pAlign)));
	return R_NilValue;
}

SEXP bam_range_pop_back(SEXP pRange)
{
	if(TYPEOF(pRange)!=EXTPTRSXP)
	{
		error("[bam_range_pop_back] pRange is No external pointer!\n");
		return R_NilValue;
	}
	align_list_pop_back((align_list*)(R_ExternalPtrAddr(pRange)));
	return R_NilValue;
}
SEXP bam_range_pop_front(SEXP pRange)
{
	if(TYPEOF(pRange)!=EXTPTRSXP)
	{
		error("[bam_range_pop_front] pRange is No external pointer!\n");
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
		error("[bam_writer_open] pHeader No external pointer!\n");
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
		error("[bam_writer_open] pReader No external pointer!\n");
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
