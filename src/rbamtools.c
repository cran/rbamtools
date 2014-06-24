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
 *				18.Mar.13 (nemesis) Corrected inline declarations; Changed qmm to qsm (mean to sum)
 */

#ifndef rbamtools_c
#define rbamtools_c
#include "rbamtools.h"

static const int buf_size=2048;  /* buffer size for printing ints into chars */

/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 *
 * Character encoding for letters which might occur in BAM file sequence segments
 *
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
# define a_val 0
# define c_val 1
# define g_val 2
# define t_val 3
# define zvl   4

static const unsigned char bam_nt16_dna_table[16] = {
//		=    A      C      M    G  	   R	S    V    T      W    Y    H    K    D    B    N
		zvl, a_val, c_val, zvl, g_val, zvl, zvl, zvl, t_val, zvl, zvl, zvl,	zvl, zvl, zvl, zvl	// 0
};
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 *
 * Character encoding for letters which might occur in FASTQ files
 *
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */



# define a_fq 0
# define c_fq 1
# define g_fq 2
# define t_fq 3
# define m_fq 4
# define r_fq 5
# define w_fq 6
# define s_fq 7
# define y_fq 8
# define k_fq 9
# define v_fq 10
# define h_fq 11
# define d_fq 12
# define b_fq 13
# define n_fq 14
# define plfq 15
# define mnfq 16
# define eqfq 17
# define zvfq 18
static const int nFastq=19;


static const unsigned char FASTQ_LETTERS[256] = {
		// 64 - block																					// start index
		zvfq,zvfq,zvfq,zvfq,	zvfq,zvfq,zvfq,zvfq,	zvfq,zvfq,zvfq,zvfq,	zvfq,zvfq,zvfq,zvfq,	// 0
		zvfq,zvfq,zvfq,zvfq,	zvfq,zvfq,zvfq,zvfq,	zvfq,zvfq,zvfq,zvfq,	zvfq,zvfq,zvfq,zvfq,	// 16
		zvfq,zvfq,zvfq,zvfq,	zvfq,zvfq,zvfq,zvfq,	zvfq,zvfq,zvfq,plfq,	zvfq,mnfq,zvfq,zvfq,	// 32
		zvfq,zvfq,zvfq,zvfq,	zvfq,zvfq,zvfq,zvfq,	zvfq,zvfq,zvfq,zvfq,	zvfq,eqfq,zvfq,zvfq,	// 48

		zvfq,a_fq,b_fq,c_fq,	d_fq,zvfq,zvfq,g_fq,	h_fq,zvfq,zvfq,k_fq,	zvfq,m_fq,n_fq,zvfq,	// 64
		zvfq,s_fq,r_fq,zvfq,	t_fq,zvfq,v_fq,w_fq,	zvfq,y_fq,zvfq,zvfq,	zvfq,zvfq,zvfq,zvfq,	// 80
		zvfq,a_fq,b_fq,c_fq,	d_fq,zvfq,zvfq,g_fq,	h_fq,zvfq,zvfq,k_fq,	zvfq,m_fq,n_fq,zvfq,	// 96
		zvfq,s_fq,r_fq,zvfq,	t_fq,zvfq,v_fq,w_fq,	zvfq,y_fq,zvfq,zvfq,	zvfq,zvfq,zvfq,zvfq,	// 112

		zvfq,zvfq,zvfq,zvfq,	zvfq,zvfq,zvfq,zvfq,	zvfq,zvfq,zvfq,zvfq,	zvfq,zvfq,zvfq,zvfq,	// 128
		zvfq,zvfq,zvfq,zvfq,	zvfq,zvfq,zvfq,zvfq,	zvfq,zvfq,zvfq,zvfq,	zvfq,zvfq,zvfq,zvfq,	// 144
		zvfq,zvfq,zvfq,zvfq,	zvfq,zvfq,zvfq,zvfq,	zvfq,zvfq,zvfq,zvfq,	zvfq,zvfq,zvfq,zvfq,	// 160
		zvfq,zvfq,zvfq,zvfq,	zvfq,zvfq,zvfq,zvfq,	zvfq,zvfq,zvfq,zvfq,	zvfq,zvfq,zvfq,zvfq,	// 176

		zvfq,zvfq,zvfq,zvfq,	zvfq,zvfq,zvfq,zvfq,	zvfq,zvfq,zvfq,zvfq,	zvfq,zvfq,zvfq,zvfq,	// 192
		zvfq,zvfq,zvfq,zvfq,	zvfq,zvfq,zvfq,zvfq,	zvfq,zvfq,zvfq,zvfq,	zvfq,zvfq,zvfq,zvfq,	// 208
		zvfq,zvfq,zvfq,zvfq,	zvfq,zvfq,zvfq,zvfq,	zvfq,zvfq,zvfq,zvfq,	zvfq,zvfq,zvfq,zvfq,	// 224
		zvfq,zvfq,zvfq,zvfq,	zvfq,zvfq,zvfq,zvfq,	zvfq,zvfq,zvfq,zvfq,	zvfq,zvfq,zvfq,zvfq	    // 240
};
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////



static R_INLINE void clear_buf(char *c,unsigned n)
{
	unsigned i;
	for(i=0;i<n;++i)
		c[i]=(char)0;
}

static R_INLINE void set_flag(bam1_t *align,_Bool val,unsigned pattern)
{
	if(val)
		align->core.flag=(align->core.flag) | pattern;
	else
		align->core.flag=(align->core.flag) & !pattern;
}

static R_INLINE int cigar2str(char *c,const bam1_t *align)
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

static R_INLINE uint8_t *alloc_align_data(bam1_t *b, int size)
{
	if (b->m_data < size) {
		b->m_data = size;
		kroundup32(b->m_data);
		b->data = (uint8_t*)realloc(b->data, b->m_data);
	}
	return b->data;
}

void print_bitmask_dec(bitmap_type *value)
{
	int i,last_block;
	last_block=bitmap_size-1;
	Rprintf("%3u",getByte(*value,last_block));
	for(i=last_block-1;i>=0;--i)
		Rprintf(" | %3u",getByte(*value,i));
	Rprintf("\n");
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
		error("[is_nil_externalptr] No external pointer");

	SEXP ans;
	PROTECT(ans=allocVector(LGLSXP,1));
	LOGICAL(ans)[0]=(R_ExternalPtrAddr(ptr)==NULL);
	UNPROTECT(1);
	return ans;
}

/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 *
 * bam_reader
 *
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
static void finalize_bam_reader(SEXP ptr)
{
	if(TYPEOF(ptr)!=EXTPTRSXP)
		error("[finalize_bam_reader] No external pointer!");
	samfile_t *reader=(samfile_t*)(R_ExternalPtrAddr(ptr));
	samclose(reader);	// checks for 0
	R_SetExternalPtrAddr(ptr,NULL);
}

SEXP bam_reader_open(SEXP filename)
{
	if(TYPEOF(filename)!=STRSXP)
		error("[bam_reader_open] Filename must be a string.\n");

	const char* _filename=CHAR(STRING_ELT(filename,0));
	samfile_t *reader=samopen(_filename,"rb",0);
	if(!reader)
		error("[bam_reader_open] Opening bam_file \"%s\" failed!",_filename);

	SEXP pReader;
	PROTECT(pReader=R_MakeExternalPtr( (void*)(reader),R_NilValue,R_NilValue));
	R_RegisterCFinalizer(pReader,finalize_bam_reader);
	UNPROTECT(1);
	return pReader;
}

SEXP bam_reader_close(SEXP pReader)
{
	if(TYPEOF(pReader)!=EXTPTRSXP)
		error("[bam_reader_close] No external pointer!");

	samfile_t *reader=(samfile_t*) (R_ExternalPtrAddr(pReader));
	samclose(reader);
	R_SetExternalPtrAddr(pReader,NULL);
	return R_NilValue;
}


SEXP bam_reader_get_header_text(SEXP pReader)
{
	if(TYPEOF(pReader)!=EXTPTRSXP)
		error("[bam_reader_get_header_text] No external pointer!");

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
		error("[bam_reader_get_ref_count] No external pointer!");

	samfile_t *reader=(samfile_t*) (R_ExternalPtrAddr(pReader));
	SEXP ans;
	PROTECT(ans=Rf_allocVector(INTSXP,1));
	INTEGER(ans)[0]=reader->header->n_targets;
	UNPROTECT(1);
	return ans;
}

SEXP bam_reader_get_target_name(SEXP pReader, SEXP pSeqid)
{
	if(TYPEOF(pReader)!=EXTPTRSXP)
		error("[bam_reader_get_target_name] No external pointer!");
	if(TYPEOF(pSeqid)!=INTSXP)
		error("[bam_reader_get_target_name] pSeqid must be INT!");

	samfile_t *reader=(samfile_t*)(R_ExternalPtrAddr(pReader));
	bam_header_t* header=reader->header;

	int seqid=INTEGER(pSeqid)[0];
	if((seqid<0) | (seqid >=header->n_targets))
		error("[bam_reader_get_target_name] pSeqid out of range!");

	SEXP res;
	PROTECT(res=allocVector(STRSXP,1));
	SET_STRING_ELT(res,0,mkChar(header->target_name[seqid]));
	UNPROTECT(1);
	return res;
}

SEXP bam_reader_get_ref_data(SEXP pReader)
{
	if(TYPEOF(pReader)!=EXTPTRSXP)
		error("[bam_reader_get_ref_data] No external pointer!");

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
		error("[bam_reader_create_index] BamFile must be a string!\n");
	if(TYPEOF(pIdxFile)!=STRSXP)
		error("[bam_reader_create_index] IndexFile must be a string!\n");

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
		error("[finalize_bam_index] No external pointer!");

	bam_index_t *index=(bam_index_t *)(R_ExternalPtrAddr(ptr));
	bam_index_destroy(index);	// checks for zero
	R_SetExternalPtrAddr(ptr,NULL);
}

SEXP bam_reader_load_index(SEXP pIdxFile)
{
	if(TYPEOF(pIdxFile)!=STRSXP)
		error("[bam_reader_load_index] pIdxFile must be a string!\n");

	const char *idxFile=CHAR(STRING_ELT(pIdxFile,0));
	FILE *f=fopen(idxFile,"rb");
	bam_index_t *index = bam_index_load_core(f);
	fclose(f);

	SEXP idx;
	PROTECT(idx=R_MakeExternalPtr( (void*)(index),R_NilValue,R_NilValue));
	R_RegisterCFinalizer(idx,finalize_bam_index);
	UNPROTECT(1);
	return idx;
}

SEXP bam_reader_unload_index(SEXP pIdx)
{
	if(TYPEOF(pIdx)!=EXTPTRSXP)
		error("[bam_reader_unload_index] No external pointer!\n");

	bam_index_t *idx=(bam_index_t *)(R_ExternalPtrAddr(pIdx));
	bam_index_destroy(idx);
	R_SetExternalPtrAddr(pIdx,NULL);
	return R_NilValue;
}

SEXP bam_reader_get_next_align(SEXP pReader)
{
	if(TYPEOF(pReader)!=EXTPTRSXP)
		error("[bam_reader_get_next_align] No external pointer!\n");

	samfile_t *reader=(samfile_t*)(R_ExternalPtrAddr(pReader));
	bam1_t *align=bam_init1();
	int res=samread(reader,align);
	if(res==-1)
	{
		Rprintf("[getNextAlign] samread found EOF.\n");
		return R_NilValue;
	}
	if(res==-2)
		error("[getNextAlign] samread found truncated BAM-file.\n");


	SEXP ptr;
	PROTECT(ptr=R_MakeExternalPtr((void*)(align),R_NilValue,R_NilValue));
	R_RegisterCFinalizer(ptr,finalize_bam_align);
	UNPROTECT(1);
	return ptr;
}

SEXP bam_reader_save_aligns(SEXP pReader,SEXP pWriter)
{
	if(TYPEOF(pReader)!=EXTPTRSXP)
		error("[bam_reader_save_aligns] pReader: No external pointer!\n");
	if(TYPEOF(pWriter)!=EXTPTRSXP)
		error("[bam_reader_save_aligns] pWriter: No external pointer!\n");

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

	bam_destroy1(align);	// checks for >0!
	if(res==-2)
		error("[bam_reader_save_aligns] samread found truncated BAM-file.\n");


	SEXP ans;
	PROTECT(ans=allocVector(INTSXP,1));
	INTEGER(ans)[0]=nAligns;
	UNPROTECT(1);
	return ans;
}

SEXP bam_reader_sort_file(SEXP pFilename,SEXP pPrefix,SEXP pMaxMem,SEXP pByName)
{
	if(TYPEOF(pFilename)!=STRSXP)
		error("[bam_writer_sort_file] Filename must be a string\n");
	if(TYPEOF(pPrefix)!=STRSXP)
		error("[bam_writer_sort_file] Prefix must be a string\n");
	if(TYPEOF(pMaxMem)!=REALSXP)
		error("[bam_writer_sort_file] MaxMem must be integer value!\n");
	if(TYPEOF(pByName)!=LGLSXP)
		error("[bam_writer_sort_file] ByName must be bool value!\n");

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
		error("[bam_reader_get_header] No external pointer!\n");

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
		error("[bam_reader_tell] No external pointer!\n");
	return Rf_ScalarReal(bam_tell(((samfile_t*)(R_ExternalPtrAddr(pReader)))->x.bam));
}

SEXP bam_reader_seek(SEXP pReader, SEXP pPos)
{
	if(TYPEOF(pReader)!=EXTPTRSXP)
		error("[bam_reader_seek] No external pointer!\n");
	if(TYPEOF(pPos)!=REALSXP)
		error("[bam_reader_seek] Position must be numeric!\n");
	if(LENGTH(pPos)>1)
		error("[bam_reader_seek] Length of position must be 1!\n");

	samfile_t *reader=(samfile_t*)(R_ExternalPtrAddr(pReader));
	double *pos=REAL(pPos);
    int64_t res=bam_seek(reader->x.bam,(int64_t)*pos,SEEK_SET);
    if(res<0)
    	Rprintf("[bam_reader_seek] bam_seek fails!\n");
	return R_NilValue;
}

SEXP bam_reader_write_fastq(SEXP pReader,SEXP pFilename,SEXP pAppend)
{
	if(TYPEOF(pReader)!=EXTPTRSXP)
		error("[bam_reader_write_fastq] No external pointer!\n");
	if(TYPEOF(pFilename)!=STRSXP)
		error("[bam_reader_write_fastq] Filename must be a string!\n");
	if(TYPEOF(pAppend)!=LGLSXP)
		error("[bam_reader_write_fastq] pAppend must be logical!\n");

	_Bool append =*(LOGICAL(AS_LOGICAL(pAppend)));

	gzFile gz;
	if(append)
		gz =gzopen(CHAR(STRING_ELT(pFilename,0)),"ab");
	else
		gz =gzopen(CHAR(STRING_ELT(pFilename,0)),"wb");
	if(!gz)
		error("[bam_reader_write_fastq] Cannot open file '%s'!",CHAR(STRING_ELT(pFilename,0)));

	/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
	unsigned char *raw_seq;
	int32_t seq_len;
	int buf_size=2048, nuc_val;
	char *buf=R_alloc(buf_size,sizeof(char));
	uint8_t *quals;
	unsigned i=0,j;

	/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
	samfile_t *reader=(samfile_t*)(R_ExternalPtrAddr(pReader));
	bam1_t *align=bam_init1();
	int res=samread(reader,align);

	while(res>=0)
	{
		// Read Sequence
		seq_len=align->core.l_qseq;
		if(seq_len>buf_size)
		{
			buf_size=2*seq_len;
			buf=R_alloc(buf_size,sizeof(char));
		}
		raw_seq=bam1_seq(align);
		for(j=0;j<seq_len;++j)
		{
			nuc_val=bam1_seqi(raw_seq,j);
			if(nuc_val<0 || nuc_val>15)
			{
				Rprintf("[bam_reader_write_fastq] bam_seq out of range at align %u\n",i+1);
				nuc_val=15;
			}
			buf[j]=bam_nt16_rev_table[nuc_val];
		}
		buf[j]='\0';

		// Write header line and sequence
		gzprintf(gz,"@%s\n%s\n+\n",bam1_qname(align),buf);

		// Print quality string
		quals=bam1_qual(align);
		for(j=0;j<seq_len;++j)
		{
			// Ensure that printed character is in range
			quals[j]= (quals[j]<0 ? 0 : quals[j]);
			quals[j]= (quals[j]>93 ? 93 : quals[j]);
			buf[j]=(char) (quals[j]+33);
		}
		gzprintf(gz,"%s\n",buf);
		++i;
		res=samread(reader,align);
	}
	gzclose(gz);
	bam_destroy1(align);	// checks for >0!
	if(res==-2)
		error("[bam_reader_write_fastq] samread found truncated BAM-file.\n");

	SEXP ans;
	PROTECT(ans=allocVector(INTSXP,1));
	INTEGER(ans)[0]=i;
	UNPROTECT(1);
	return ans;
}

SEXP bam_reader_write_fastq_index(SEXP pReader,SEXP pFilename,SEXP pWhichCopy,SEXP pAppend)
{
	if(TYPEOF(pReader)!=EXTPTRSXP)
		error("[bam_reader_write_fastq_index] No external pointer!\n");
	if(TYPEOF(pFilename)!=STRSXP)
		error("[bam_reader_write_fastq_index] Filename must be a string!\n");
	if(TYPEOF(pWhichCopy)!=INTSXP)
		error("[bam_reader_write_fastq_index] pWichWrite must be integer!\n");
	if(TYPEOF(pAppend)!=LGLSXP)
		error("[bam_reader_write_fastq_index] pAppend must be logical!\n");

	/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
	 *
	 * The function reads remaining aligns from reader
	 * until EOF is reached or length(pWichWrite) aligns are checked.
	 * For each align, pWhichCopy[i] is checked.
	 * When TRUE, then align is written to file, otherwise skipped.
	 *
	 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

	unsigned nCheck=INTEGER(pWhichCopy)[LENGTH(pWhichCopy)-1];
	_Bool append =*(LOGICAL(AS_LOGICAL(pAppend)));

	gzFile gz;
	if(append)
		gz =gzopen(CHAR(STRING_ELT(pFilename,0)),"a+");
	else
		gz =gzopen(CHAR(STRING_ELT(pFilename,0)),"w");
	if(!gz)
		error("[bam_reader_write_fastq_index] Cannot open file '%s'!",CHAR(STRING_ELT(pFilename,0)));


	/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
	 * char buffer for copying sequence and qualities.
	 */
	unsigned char *raw_seq;
	int32_t seq_len;
	int buf_size=2048,nuc_val;
	char *buf=R_alloc(buf_size,sizeof(char));
	uint8_t *quals;
	unsigned nWritten=0,j,nChecked=0;

	/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
	 */
	samfile_t *reader=(samfile_t*)(R_ExternalPtrAddr(pReader));
	bam1_t *align=bam_init1();
	int res=samread(reader,align);

	while(res>=0 && nChecked < nCheck)
	{
		++nChecked;
		if(nChecked==(INTEGER(pWhichCopy)[nWritten]))
		{

			// read sequence
			seq_len=align->core.l_qseq;
			if(seq_len>buf_size)
			{
				buf_size=2*seq_len;
				buf= R_alloc(buf_size,sizeof(char));
			}
			raw_seq=bam1_seq(align);
			for(j=0;j<seq_len;++j)
			{
				nuc_val=bam1_seqi(raw_seq,j);
				if(nuc_val<0 || nuc_val>15)
				{
					Rprintf("[bam_reader_write_fastq_index] bam_seq out of range at align %u\n",nChecked);
					nuc_val=15;
				}
				buf[j]=bam_nt16_rev_table[nuc_val];
			}
			buf[j]='\0';

			// write header line and sequence
			gzprintf(gz,"@%s\n%s\n+\n",bam1_qname(align),buf);

			// Print quality string
			quals=bam1_qual(align);
			for(j=0;j<seq_len;++j)
			{
				// Ensure that printed character is in range
				quals[j]= (quals[j]<0 ? 0 : quals[j]);
				quals[j]= (quals[j]>93 ? 93 : quals[j]);
				buf[j]=(char) (quals[j]+33);
			}
			gzprintf(gz,"%s\n",buf);
			++nWritten;
		}
		res=samread(reader,align);
	}

	gzclose(gz);
	bam_destroy1(align);	// checks for >0!
	align=0;

	if(res==-2)
		error("[bam_reader_write_fastq_index] samread found truncated BAM-file.\n");

	SEXP ans;
	PROTECT(ans=allocVector(INTSXP,2));
	INTEGER(ans)[0]=nWritten;
	INTEGER(ans)[1]=nChecked;
	UNPROTECT(1);
	return ans;
}


/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 *
 * GapList
 *
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

static void finalize_gap_list(SEXP ptr)
{
	if(TYPEOF(ptr)!=EXTPTRSXP)
		error("[finalize_gap_list] No external pointer!\n");
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
		error("[gap_list_fetch] pReader is No external pointer!\n");
	if(TYPEOF(pIndex)!=EXTPTRSXP)
		error("[gap_list_fetch] pIndex is No external pointer!\n");
	if(TYPEOF(pCoords)!=REALSXP)
		error("[gap_list_fetch] pCoords is no REAL!\n");
	if(LENGTH(pCoords)!=3)
		error("[gap_list_fetch] pCoords must contain three values (refid,begin,end)!\n");

	samfile_t *reader=(samfile_t*)(R_ExternalPtrAddr(pReader));
	bam_index_t *index=(bam_index_t*)(R_ExternalPtrAddr(pIndex));
	if(reader==NULL)
		error("[gap_list_fetch] Reader must not be NULL pointer!\n");

	if(index==NULL)
		error("[gap_list_fetch] Index must not be NULL pointer!\n");

	double *pi=REAL(pCoords);
	int refid=(int) pi[0];
	int begin=(int) pi[1];
	int end=(int) pi[2];

	if(refid<0 || refid >=(reader->header->n_targets))
		error("[gap_list_fetch] refid out of range!\n");
	if(begin<0 || begin>=end || end>(reader->header->target_len[refid]))
		error("[gap_list_fetch] Begin or end out of range!\n");


	gap_list *l=init_gap_list();
    bam_fetch(reader->x.bam, index, refid, begin, end, (void*)l, gap_fetch_func);
    SEXP list;
	PROTECT(list=R_MakeExternalPtr( (void*)(l),R_NilValue,R_NilValue));
	R_RegisterCFinalizer(list,finalize_gap_list);
	UNPROTECT(1);
	return list;
}

SEXP gap_list_get_df(SEXP pGapList)
{
	if(TYPEOF(pGapList)!=EXTPTRSXP)
		error("[get_gap_list_df] No external pointer!");

	gap_list *l=(gap_list*)(R_ExternalPtrAddr(pGapList));

	/* create data.frame 									*/
	int nProtected=0;
	int nCols=9;
	SEXP dflist;
	PROTECT(dflist=allocVector(VECSXP,nCols));
	++nProtected;
	int nRows=(l->size);
	int i;
	unsigned *pos_data;

	/* Column 0: refid 										*/
	SEXP ref_vector;
	PROTECT(ref_vector=allocVector(INTSXP,nRows));
	++nProtected;

	/* Column 1: position 									*/
	SEXP pos_vector;
	PROTECT(pos_vector=allocVector(INTSXP,nRows));
	++nProtected;

	/* Column 2: left_cigar_len								*/
	SEXP left_cigar_len_vector;
	PROTECT(left_cigar_len_vector=allocVector(INTSXP,nRows));
	++nProtected;

	/* Column 3: left_cigar_type							*/
	SEXP left_cigar_type_vector;
	PROTECT(left_cigar_type_vector=allocVector(INTSXP,nRows));
	++nProtected;

	/* Column 4: left_stop									*/
	SEXP left_stop_vector;
	PROTECT(left_stop_vector=allocVector(INTSXP,nRows));
	++nProtected;

	/* Column 5: gap_len									*/
	SEXP gap_len_vector;
	PROTECT(gap_len_vector=allocVector(INTSXP,nRows));
	++nProtected;

	/* Column 6: right_start								*/
	SEXP right_start_vector;
	PROTECT(right_start_vector=allocVector(INTSXP,nRows));
	++nProtected;

	/* Column 7: right_cigar_len							*/
	SEXP right_cigar_len_vector;
	PROTECT(right_cigar_len_vector=allocVector(INTSXP,nRows));
	++nProtected;

	/* Column 8: right_cigar_type							*/
	SEXP right_cigar_type_vector;
	PROTECT(right_cigar_type_vector=allocVector(INTSXP,nRows));
	++nProtected;

	/* start position										*/
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

	/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
	 * Convert left_cigar_type to vector
	 */

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

	/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
	 * Convert right_cigar_type_vector to vector
	 */

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
	/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

	SET_VECTOR_ELT(dflist,0,ref_vector);
	SET_VECTOR_ELT(dflist,1,pos_vector);
	SET_VECTOR_ELT(dflist,2,left_cigar_len_vector);
	SET_VECTOR_ELT(dflist,3,left_cigar_type_vector);
	SET_VECTOR_ELT(dflist,4,left_stop_vector);
	SET_VECTOR_ELT(dflist,5,gap_len_vector);
	SET_VECTOR_ELT(dflist,6,right_start_vector);
	SET_VECTOR_ELT(dflist,7,right_cigar_len_vector);
	SET_VECTOR_ELT(dflist,8,right_cigar_type_vector);

	/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
	 * Column Names
	 */
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

	/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
	 * Row Names
	 */
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
SEXP gap_list_get_nAlignGaps(SEXP pGapList)
{
	if(TYPEOF(pGapList)!=EXTPTRSXP)
		error("[gap_list_get_size] No external pointer!");

	gap_list *l=(gap_list*)(R_ExternalPtrAddr(pGapList));
	SEXP ans;
	PROTECT(ans=allocVector(INTSXP,1));
	INTEGER(ans)[0]=(l->nAlignGaps);
	UNPROTECT(1);
	return ans;
}


/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 * gap_site_list
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

static void finalize_gap_site_list(SEXP ptr)
{
	if(TYPEOF(ptr)!=EXTPTRSXP)
		error("[finalize_gap_list] No external pointer!\n");
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
		error("[gap_site_list_fetch] pReader is No external pointer!\n");
	if(TYPEOF(pIndex)!=EXTPTRSXP)
		error("[gap_site_list_fetch] pIndex is No external pointer!\n");
	if(TYPEOF(pCoords)!=REALSXP)
		error("[gap_site_list_fetch] pCoords is no REAL!\n");
	if(LENGTH(pCoords)!=3)
		error("[gap_site_list_fetch] pCoords must contain three values (refid,begin,end)!\n");

	samfile_t *reader=(samfile_t*)(R_ExternalPtrAddr(pReader));
	bam_index_t *index=(bam_index_t*)(R_ExternalPtrAddr(pIndex));
	if(reader==NULL)
		error("[gap_site_list_fetch] Reader must not be NULL pointer!\n");
	if(index==NULL)
		error("[gap_site_list_fetch] Index must not be NULL pointer!\n");

	double *pi=REAL(pCoords);
	int refid=(int) pi[0];
	int begin=(int) pi[1];
	int end=(int) pi[2];

	if(refid<0 || refid >=(reader->header->n_targets))
		error("[gap_site_list_fetch] refid out of range!\n");
	if(begin<0 || begin>=end || end>(reader->header->target_len[refid]))
		error("[gap_site_list_fetch] Begin or end out of range!\n");

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

	/* Column 0: id												*/
	SEXP id_vector;
	PROTECT(id_vector=allocVector(INTSXP,nRows));
	++nProtected;
	// Column 1: refid											*/
	SEXP refid_vector;
	PROTECT(refid_vector=allocVector(INTSXP,nRows));
	++nProtected;
	/* Column 2: lstart											*/
	SEXP lstart_vector;
	PROTECT(lstart_vector=allocVector(INTSXP,nRows));
	++nProtected;
	/* Column 3: lend											*/
	SEXP lend_vector;
	PROTECT(lend_vector=allocVector(INTSXP,nRows));
	++nProtected;
	/* Column 4: rstart											*/
	SEXP rstart_vector;
	PROTECT(rstart_vector=allocVector(INTSXP,nRows));
	++nProtected;
	/* Column 5: rend											*/
	SEXP rend_vector;
	PROTECT(rend_vector=allocVector(INTSXP,nRows));
	++nProtected;
	/* Column 6: gap_len										*/
	SEXP gap_len_vector;
	PROTECT(gap_len_vector=allocVector(INTSXP,nRows));
	++nProtected;
	/* Column 7: nAligns										*/
	SEXP nAligns_vector;
	PROTECT(nAligns_vector=allocVector(INTSXP,nRows));
	++nProtected;
	/* Column 8: nProbes										*/
	SEXP nProbes_vector;
	PROTECT(nProbes_vector=allocVector(INTSXP,nRows));
	++nProtected;
	/* Column 9: nlstart										*/
	SEXP nlstart_vector;
	PROTECT(nlstart_vector=allocVector(INTSXP,nRows));
	++nProtected;
	/* Column 10: sum_left_cigar								*/
	SEXP qsm_vector;
	PROTECT(qsm_vector=allocVector(INTSXP,nRows));
	++nProtected;
	/* Column 11: lcl											*/
	SEXP nmcl_vector;
	PROTECT(nmcl_vector=allocVector(INTSXP,nRows));
	++nProtected;
	/* Column 12: mcl											*/
	SEXP mcl_vector;
	PROTECT(mcl_vector=allocVector(INTSXP,nRows));
	++nProtected;


	/* start position											*/
	site_list_element *el,*curr;
	curr=l->curr;
	site_list_set_curr_first(l);

	for(i=0;i<nRows;++i)
	{
		el=site_list_get_curr_pp(l);
		INTEGER(id_vector)                  [i]=i+1;
		INTEGER(refid_vector)        		[i]=l->refid;
		/* lend - max(left_cigar_len)+1							*/
		INTEGER(lstart_vector)        		[i]=el->lend-getByte(el->lcl,0)+1;
		INTEGER(lend_vector)	            [i]=el->lend;
		INTEGER(rstart_vector)	            [i]=el->rstart;
		INTEGER(rend_vector)  		        [i]=el->rstart+el->r_cigar_size-1;
		INTEGER(gap_len_vector)			    [i]=el->gap_len;
		INTEGER(nAligns_vector)		        [i]=el->nAligns;
		INTEGER(nProbes_vector)             [i]=el->nProbes;
		INTEGER(nlstart_vector)	            [i]=bitmask_nPos(el->lcl);
		INTEGER(qsm_vector)					[i]=bitmask_sumPos(el->lcl);
		INTEGER(nmcl_vector)				[i]=el->lcl;
		INTEGER(mcl_vector)					[i]=el->mcl;
	}
	/* Reset curr position										*/
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
	SET_VECTOR_ELT(dflist,10,qsm_vector);
	SET_VECTOR_ELT(dflist,11,nmcl_vector);
	SET_VECTOR_ELT(dflist,12,mcl_vector);

	/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
	 * Column Names
	 */
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
	SET_STRING_ELT(col_names,11,mkChar("lcl"));
	SET_STRING_ELT(col_names,12,mkChar("mcl"));
	setAttrib(dflist,R_NamesSymbol,col_names);

	/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
	 * Row Names
	 */
	SEXP row_names;
    PROTECT(row_names=allocVector(STRSXP,nRows));
    ++nProtected;

    int buf_size=128;
    char *buf=(char*) calloc(buf_size,sizeof(char));
    for(i=0;i<nRows;++i)
    {
    	sprintf(buf,"%i",i+1);
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

SEXP gap_site_list_get_nAlignGaps(SEXP pGapList)
{
	if(TYPEOF(pGapList)!=EXTPTRSXP)
		error("[gap_site_list_get_nAlignGaps] No external pointer!");

	site_list *l=(site_list*)(R_ExternalPtrAddr(pGapList));
	SEXP ans;
	PROTECT(ans=allocVector(INTSXP,1));
	INTEGER(ans)[0]=(l->nAlignGaps);
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

/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 * gap_site_ll
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

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
		error("[gap_site_ll_fetch] Reader must not be NULL pointer!\n");
	if(index==NULL)
		error("[gap_site_ll_fetch] Index must not be NULL pointer!\n");

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
	l->nAlignGaps=ls->nAlignGaps+rs->nAlignGaps;

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

	// Column 3: nAlignGaps
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
		INTEGER(ngap_vector)  [i]=l->nAlignGaps;
	}

	SET_VECTOR_ELT(dflist, 0,refid_vector);
	SET_VECTOR_ELT(dflist, 1,size_vector);
	SET_VECTOR_ELT(dflist, 2,nalign_vector);
	SET_VECTOR_ELT(dflist, 3,ngap_vector);

	/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
	 * Column Names
	 */
	SEXP col_names;
	PROTECT(col_names=allocVector(STRSXP,nCols));
	++nProtected;

	SET_STRING_ELT(col_names, 0,mkChar("ID"));
	SET_STRING_ELT(col_names, 1,mkChar("size"));
	SET_STRING_ELT(col_names, 2,mkChar("nAligns"));
	SET_STRING_ELT(col_names, 3,mkChar("nAlignGaps"));
	setAttrib(dflist,R_NamesSymbol,col_names);

	/* Row Names											*/
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

	/* create data.frame											*/
	unsigned nProtected=0;
	unsigned nCols=13;
	unsigned nRefs=sll->size;
	unsigned nRows=sum_ll_sizes(sll);

	SEXP dflist;
	PROTECT(dflist=allocVector(VECSXP,nCols));
	++nProtected;

	/* Column 0: id													*/
	SEXP id_vector;
	PROTECT(id_vector=allocVector(INTSXP,nRows));
	++nProtected;
	/* Column 1: refid												*/
	SEXP refid_vector;
	PROTECT(refid_vector=allocVector(INTSXP,nRows));
	++nProtected;
	// Column 2: lstart
	SEXP lstart_vector;
	PROTECT(lstart_vector=allocVector(INTSXP,nRows));
	++nProtected;
	/* Column 3: lend												*/
	SEXP lend_vector;
	PROTECT(lend_vector=allocVector(INTSXP,nRows));
	++nProtected;
	/* Column 4: rstart												*/
	SEXP rstart_vector;
	PROTECT(rstart_vector=allocVector(INTSXP,nRows));
	++nProtected;
	/* Column 5: rend												*/
	SEXP rend_vector;
	PROTECT(rend_vector=allocVector(INTSXP,nRows));
	++nProtected;
	/* Column 6: gap_len											*/
	SEXP gap_len_vector;
	PROTECT(gap_len_vector=allocVector(INTSXP,nRows));
	++nProtected;
	/* Column 7: nAligns											*/
	SEXP nAligns_vector;
	PROTECT(nAligns_vector=allocVector(INTSXP,nRows));
	++nProtected;
	/* Column 8: nAligns											*/
	SEXP nProbes_vector;
	PROTECT(nProbes_vector=allocVector(INTSXP,nRows));
	++nProtected;
	/* Column 9: nlstart: number of different lstart positions on gap_site	*/
	SEXP nlstart_vector;
	PROTECT(nlstart_vector=allocVector(INTSXP,nRows));
	++nProtected;

	/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
	 *
	 * Column 10: qsm=Quadrupel mcl sum of minimal cigar size
	 * mcl contains minimum cigar size values (of flanking left and right M=match segment)
	 * Sum of Quadrupel of rightmost values (=4 largest values) in mcl
	 *
	 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
	SEXP qsm_vector;
	PROTECT(qsm_vector=allocVector(INTSXP,nRows));
	++nProtected;

	/* Column 11: nmcl												*/
	SEXP nmcl_vector;
	PROTECT(nmcl_vector=allocVector(INTSXP,nRows));
	++nProtected;

	/* Column 12: gqs												*/
	SEXP gqs_vector;
	PROTECT(gqs_vector=allocVector(INTSXP,nRows));
	++nProtected;

	/* site_ll_element *ell;										*/
	site_list *l;
	site_list_element *el;

	site_ll_set_curr_first(sll);
	unsigned i,j,k=0;

	/* Scaling factor: Maximum value for gqs is read-length			*/
	unsigned smcl_denom=bitmap_size*nQsm;

	for(i=0;i<nRefs;++i)
	{
		l=site_ll_get_curr_site_list_pp(sll);
		site_list_set_curr_first(l);
		for(j=0;j<l->size;++j,++k)
		{
			el=site_list_get_curr_pp(l);
			INTEGER(id_vector)                  [k]=k+1; 						/* row-id: should start with 1	*/
			INTEGER(refid_vector)        		[k]=l->refid+1;
			INTEGER(lstart_vector)        		[k]=el->lend-getByte(el->lcl,0)+1;
			INTEGER(lend_vector)	            [k]=el->lend;
			INTEGER(rstart_vector)	            [k]=el->rstart;
			INTEGER(rend_vector)  		        [k]=el->rstart+el->r_cigar_size-1;
			INTEGER(gap_len_vector)			    [k]=el->gap_len;
			INTEGER(nAligns_vector)		        [k]=el->nAligns;
			INTEGER(nProbes_vector)             [k]=el->nProbes;
			INTEGER(nlstart_vector)	            [k]=bitmask_nPos(el->lcl);
			INTEGER(qsm_vector)					[k]=getSmcl(el->mcl); 			/* smcl = sum of minimal cigar size	*/
			INTEGER(nmcl_vector)				[k]=bitmask_nPos(el->mcl);
			INTEGER(gqs_vector)                 [k]=INTEGER(nlstart_vector)[k]*INTEGER(qsm_vector)[k]*2*10/smcl_denom;
		}
	}

	/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
	 * Convert left_cigar_type to vector
	 */

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
	SET_VECTOR_ELT(dflist,10,qsm_vector);
	SET_VECTOR_ELT(dflist,11,nmcl_vector);
	SET_VECTOR_ELT(dflist,12,gqs_vector);

	/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
	 * Column Names
	 */
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
	SET_STRING_ELT(col_names,10,mkChar("qsm"));
	SET_STRING_ELT(col_names,11,mkChar("nmcl"));
	SET_STRING_ELT(col_names,12,mkChar("gqs"));
	setAttrib(dflist,R_NamesSymbol,col_names);

	/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
	 * Row Names
	 */
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

SEXP gap_site_ll_get_nAlignGaps(SEXP pGapList)
{
	if(TYPEOF(pGapList)!=EXTPTRSXP)
		error("gap_site_ll_get_nAlignGaps] No external pointer!");
	site_ll *l=(site_ll*) (R_ExternalPtrAddr(pGapList));
	SEXP ans;
	PROTECT(ans=allocVector(REALSXP,1));
	REAL(ans)[0]=(double)get_nAlignGaps(l);
	UNPROTECT(1);
	return ans;
}


/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 *
 *  BamRange
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

static void finalize_bam_range(SEXP ptr)
{
	if(TYPEOF(ptr)!=EXTPTRSXP)
		error("[finalize_bam_range] No external pointer!\n");
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
		error("[bam_range_fetch] pReader is No external pointer!\n");
	if(TYPEOF(pIndex)!=EXTPTRSXP)
		error("[bam_range_fetch] pIndex is No external pointer!\n");
	if(TYPEOF(pCoords)!=REALSXP)
		error("[bam_range_fetch] pCoords is no REAL!\n");
	if(LENGTH(pCoords)!=3)
		error("[bam_range_fetch] pCoords must contain three values (refid,begin,end)!\n");
	if(TYPEOF(pComplex)!=LGLSXP)
		error("[bam_range_fetch] pComplex must be logical!\n");

	/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
	 * Extract reader and index pointer
	 */
	samfile_t *reader=(samfile_t*)(R_ExternalPtrAddr(pReader));
	bam_index_t *index=(bam_index_t*)(R_ExternalPtrAddr(pIndex));
	if(reader==NULL)
		error("[bam_range_fetch] Reader must not be NULL pointer!\n");
	if(index==NULL)
		error("[bam_range_fetch] Index must not be NULL pointer!\n");


	/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
	 *  Get coordinates
	 */
	double *pi=REAL(pCoords);
	int refid=(int) pi[0];
	int begin=(int) pi[1];
	int end=(int) pi[2];

	if(refid<0 || refid >=(reader->header->n_targets))
		error("[bam_range_fetch] refid out of range!\n");
	if(begin<0 || begin>=end || end>(reader->header->target_len[refid]))
		error("[bam_range_fetch] Begin or end out of range!\n");

	/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
	 * Init align_list
	 */
	align_list *l=init_align_list();
	l->range_begin=begin;
	l->range_end=end;
	l->seqid=refid;
	l->complex= (LOGICAL(pComplex)[0]==TRUE) ? 1 : 0;
	
	/* seqname and LN												*/
	bam_header_t* header=reader->header;
	int nchar = strlen(header->target_name[refid]);
	l->refname=calloc(nchar+1,sizeof(char));
	strcpy(l->refname,header->target_name[refid]);
	l->seq_LN=header->target_len[refid];

	/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
	 * Fetch: Retrieve only complex aligns (nCigar>1) when pComplex is set:
	 */
	if(LOGICAL(pComplex)[0]==TRUE)
		bam_fetch(reader->x.bam, index, refid, begin, end, (void*)l, range_fetch_complex_func);
	else
		bam_fetch(reader->x.bam, index, refid, begin, end, (void*)l, range_fetch_func);

	// Wind back
    l->curr_el=NULL;

	/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
	 * Create bamRange S4 object
	 */
	SEXP bam_range;
    PROTECT(bam_range = NEW_OBJECT(MAKE_CLASS("bamRange")));

    SEXP list;
	PROTECT(list=R_MakeExternalPtr( (void*)(l),R_NilValue,R_NilValue));
	R_RegisterCFinalizer(list,finalize_bam_range);
	SEXP list_name;
	PROTECT(list_name=Rf_mkString("range"));
    bam_range=SET_SLOT(bam_range,list_name,list);


    /* Return													*/
	UNPROTECT(3);
	return bam_range;
}


static int bam_count_fetch_func(const bam1_t *align, void *data)
{
	unsigned long * pCount=(unsigned long*) data;
	/* Total number of aligns									*/
	++(pCount[9]);

	uint32_t i;
	uint32_t *cigar=bam1_cigar(align);

	/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
	 * Count number of aligns for each Cigar-OP type
	 * M=0, I=1, D=2, N=3, S=4, H=5, P=6, '='=7, X=8
	 */
	//
	for(i=0;i<align->core.n_cigar;++i)
		++(pCount[(cigar[i] & BAM_CIGAR_MASK)]);
	return 0;
}

SEXP bam_count(SEXP pReader,SEXP pIndex,SEXP pCoords)
{
	if(TYPEOF(pReader)!=EXTPTRSXP)
		error("[bam_count] pReader is No external pointer!\n");
	if(TYPEOF(pIndex)!=EXTPTRSXP)
		error("[bam_count] pIndex is No external pointer!\n");
	if(TYPEOF(pCoords)!=REALSXP)
		error("[bam_count] pCoords is no REAL!\n");
	if(LENGTH(pCoords)!=3)
		error("[bam_count] pCoords must contain three values (refid,begin,end)!\n");

	samfile_t *reader=(samfile_t*)(R_ExternalPtrAddr(pReader));
	bam_index_t *index=(bam_index_t*)(R_ExternalPtrAddr(pIndex));
	if(reader==NULL)
		error("[bam_count] Reader must not be NULL pointer!\n");
	if(index==NULL)
		error("[bam_count] Index must not be NULL pointer!\n");

	double *pi=REAL(pCoords);
	int refid=(int) pi[0];
	int begin=(int) pi[1];
	int end=(int) pi[2];

	if(refid<0 || refid >=(reader->header->n_targets))
		error("[bam_count] refid out of range!\n");
	if(begin<0 || begin>=end || end>(reader->header->target_len[refid]))
		error("[bam_count] Begin or end out of range!\n");


	// M=0, I=1, D=2,N=3,S=4,H=5,P=6,'='=7,X=8, all_aligns=9
	unsigned long count[10];
	int i;

	for(i=0;i<10;++i)
		count[i]=0;

	bam_fetch(reader->x.bam, index, refid, begin, end, (void*)&count, bam_count_fetch_func);

	SEXP ans;
	PROTECT(ans=allocVector(INTSXP,10));
	for(i=0;i<10;++i)
		INTEGER(ans)[i]=count[i];
	UNPROTECT(1);
	return ans;
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
		error("[bam_range_push_back] pRange is No external pointer!\n");
	if(TYPEOF(pAlign)!=EXTPTRSXP)
		error("[bam_range_push_back] pAlign is No external pointer!\n");
	write_current_align((align_list*)(R_ExternalPtrAddr(pRange)),(bam1_t*) (R_ExternalPtrAddr(pAlign)));
	return R_NilValue;
}

SEXP bam_range_insert_past_curr_align(SEXP pRange,SEXP pAlign)
{
	if(TYPEOF(pRange)!=EXTPTRSXP)
		error("[bam_range_push_back] pRange is No external pointer!\n");
	if(TYPEOF(pAlign)!=EXTPTRSXP)
		error("[bam_range_push_back] pAlign is No external pointer!\n");
	insert_past_curr_align((align_list*)(R_ExternalPtrAddr(pRange)),(bam1_t*) (R_ExternalPtrAddr(pAlign)));
	return R_NilValue;
}

SEXP bam_range_insert_pre_curr_align(SEXP pRange,SEXP pAlign)
{
	if(TYPEOF(pRange)!=EXTPTRSXP)
		error("[bam_range_push_back] pRange is No external pointer!\n");
	if(TYPEOF(pAlign)!=EXTPTRSXP)
		error("[bam_range_push_back] pAlign is No external pointer!\n");
	insert_pre_curr_align((align_list*)(R_ExternalPtrAddr(pRange)),(bam1_t*) (R_ExternalPtrAddr(pAlign)));
	return R_NilValue;
}

SEXP bam_range_mv_curr_align(SEXP pSrc, SEXP pTarget)
{
	/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
	 * Moves current align in src to end of target list
	 * and moves current align in src to next align
	 */
	if(TYPEOF(pSrc)!=EXTPTRSXP)
		error("[bam_range_mv_curr_align] pSrc is No external pointer!\n");
	if(TYPEOF(pTarget)!=EXTPTRSXP)
		error("[bam_range_mv_curr_align] pTarget is No external pointer!\n");
	align_list_mv_curr_elem((align_list*)(R_ExternalPtrAddr(pSrc)),(align_list*)(R_ExternalPtrAddr(pTarget)));
	return R_NilValue;
}

SEXP bam_range_get_align_df(SEXP pRange)
{
	if(TYPEOF(pRange)!=EXTPTRSXP)
		error("[bam_range_get_align_df] No external pointer!");

	align_list *l=(align_list*)(R_ExternalPtrAddr(pRange));
	/* Read Adress of current align->for reconstitution at the end			*/
	align_element *e=l->curr_el;
	wind_back(l);
	const bam1_t *align;

	/* create data.frame									*/
	int nProtected=0;
	int nCols=9;
	SEXP dflist;
	PROTECT(dflist=allocVector(VECSXP,nCols));
	++nProtected;
	int nRows=(l->size);
	int i,j,nuc_val;

	/* Column 0: refid										*/
	SEXP ref_vector;
	PROTECT(ref_vector=allocVector(INTSXP,nRows));
	++nProtected;

	/* Column 1: position									*/
	SEXP pos_vector;
	PROTECT(pos_vector=allocVector(INTSXP,nRows));
	++nProtected;

	/* Column 2: nCigar										*/
	SEXP nCigar_vector;
	PROTECT(nCigar_vector=allocVector(INTSXP,nRows));
	++nProtected;

	/* Column 3: cigar										*/
	SEXP cig_vector;
	PROTECT(cig_vector=allocVector(STRSXP,nRows));
	++nProtected;

	/* Column 4: flag										*/
	SEXP flag_vector;
	PROTECT(flag_vector=allocVector(INTSXP,nRows));
	++nProtected;

	/* Column 5: seq										*/
	SEXP seq_vector;
	PROTECT(seq_vector=allocVector(STRSXP,nRows));
	++nProtected;

	/* Column 6: qual										*/
	SEXP qual_vector;
	PROTECT(qual_vector=allocVector(STRSXP,nRows));
	++nProtected;

	/* Column 7: name										*/
	SEXP name_vector;
	PROTECT(name_vector=allocVector(STRSXP,nRows));
	++nProtected;

	/* Column 8: strand_reverse								*/
	SEXP strev_vector;
	PROTECT(strev_vector=allocVector(LGLSXP,nRows));
	++nProtected;

	/* seq+cigar											*/
	unsigned char *raw_seq;
	int32_t seq_len=l->max_seqlen+1;
	char *buf=(char*) calloc(buf_size>seq_len ? buf_size : seq_len,sizeof(char));

	uint8_t *quals;

	for(i=0;i<nRows;++i)
	{
		align=get_const_next_align(l);
		INTEGER(ref_vector)[i]=(align->core.tid);
		INTEGER(pos_vector)[i]=(align->core.pos);
		INTEGER(nCigar_vector)[i]=(align->core.n_cigar);

		/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
		 * Cigar String
		 */
		if(cigar2str(buf,align)==0)
			error("[bam_range_get_align_df] Cigar error!\n");

		SET_STRING_ELT(cig_vector,i,mkChar(buf));
		//Rprintf("%i\t%s\n",i,buf);
		clear_buf(buf,buf_size);


		INTEGER(flag_vector)[i]=(align->core.flag);
		seq_len=align->core.l_qseq;
		raw_seq=bam1_seq(align);
		for(j=0;j<seq_len;++j)
		{
			nuc_val=bam1_seqi(raw_seq,j);
			if(nuc_val<0 || nuc_val>15)
			{
				Rprintf("[bam_range_get_align_df] bam_seq out of range at align %u\n",i);
				nuc_val=15;
			}
			buf[j]=bam_nt16_rev_table[nuc_val];
		}
		buf[j]=0;

		SET_STRING_ELT(seq_vector,i,mkChar(buf));

		/* quals										*/
		quals=bam1_qual(align);
		for(j=0;j<seq_len;++j)
		{
			// Ensure that printed character is in range
			quals[j]= (quals[j]<0 ? 0 : quals[j]);
			quals[j]= (quals[j]>93 ? 93 : quals[j]);
			buf[j]=(char) (quals[j]+33);
		}

		for(j=0;j<seq_len;++j)
			buf[j]=(char) (quals[j]+33);
		buf[j]=0;
		SET_STRING_ELT(qual_vector,i,mkChar(buf));

		/* read name									*/
		SET_STRING_ELT(name_vector,i,mkChar(bam1_qname(align)));

		/* strand_info									*/
		LOGICAL(strev_vector)[i]=bam1_strand(align);
	}

	/* Restore curr_el pointer							*/
	l->curr_el=e;

	SET_VECTOR_ELT(dflist,0,ref_vector);
	SET_VECTOR_ELT(dflist,1,pos_vector);
	SET_VECTOR_ELT(dflist,2,nCigar_vector);
	SET_VECTOR_ELT(dflist,3,cig_vector);
	SET_VECTOR_ELT(dflist,4,flag_vector);
	SET_VECTOR_ELT(dflist,5,seq_vector);
	SET_VECTOR_ELT(dflist,6,qual_vector);
	SET_VECTOR_ELT(dflist,7,name_vector);
	SET_VECTOR_ELT(dflist,8,strev_vector);

	/* Column Names										*/
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
	SET_STRING_ELT(col_names,7,mkChar("name"));
	SET_STRING_ELT(col_names,8,mkChar("revstrand"));
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

SEXP bam_range_write(SEXP pWriter,SEXP pRange,SEXP pRefid)
{
	if(TYPEOF(pRange)!=EXTPTRSXP)
		error("[bam_range_write] pRange is No external pointer!");
	if(TYPEOF(pWriter)!=EXTPTRSXP)
		error("[bam_range_write] pWriter is No external pointer!");
	if(TYPEOF(pRefid)!=INTSXP)
		error("[bam_range_write] pRefid must be integer!");

	unsigned long bytes_written=0;
	unsigned long range_size, i;

	samfile_t *writer=(samfile_t*) R_ExternalPtrAddr(pWriter);
	align_list *l=(align_list*)(R_ExternalPtrAddr(pRange));

	range_size=l->size;

	/* For restoration										*/
	align_element *e=l->curr_el;
	wind_back(l);

	bam1_t *align;
	for(i=0;i<range_size;++i)
	{
		align=get_next_align(l);
		align->core.tid=INTEGER(pRefid)[0];
		bytes_written+=samwrite(writer,align);
		bam_destroy1(align);
	}

	/* restore curr_el										*/
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
		error("[bam_range_wind_back] pRange is No external pointer!\n");
	align_list *l=(align_list*)(R_ExternalPtrAddr(pRange));
	wind_back(l);
	return R_NilValue;
}

SEXP bam_range_get_size(SEXP pRange)
{
	if(TYPEOF(pRange)!=EXTPTRSXP)
		error("[bam_range_get_size] pRange is No external pointer!\n");
	align_list *l=(align_list*)(R_ExternalPtrAddr(pRange));
	SEXP ans;
	PROTECT(ans=allocVector(INTSXP,1));
	INTEGER(ans)[0]=(l->size);
	UNPROTECT(1);
	return ans;
}

SEXP bam_range_get_coords(SEXP pRange)
{
	if(TYPEOF(pRange)!=EXTPTRSXP)
		error("[bam_range_get_coords] pRange is No external pointer!\n");
	align_list *l=(align_list*)(R_ExternalPtrAddr(pRange));

	SEXP ans;
	PROTECT(ans=allocVector(INTSXP,3));
	INTEGER(ans)[0]=l->seqid;
	INTEGER(ans)[1]=l->range_begin;
	INTEGER(ans)[2]=l->range_end;

	/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
	 * Names
	 */
	SEXP names_vec;
	PROTECT(names_vec=allocVector(STRSXP,3));

	SET_STRING_ELT(names_vec,0,mkChar("seqid"));
	SET_STRING_ELT(names_vec,1,mkChar("begin"));
	SET_STRING_ELT(names_vec,2,mkChar("end"));
	setAttrib(ans,R_NamesSymbol,names_vec);

	UNPROTECT(2);
	return ans;
}

SEXP bam_range_get_params(SEXP pRange)
{
	if(TYPEOF(pRange)!=EXTPTRSXP)
		error("[bam_range_get_params] pRange is No external pointer!\n");
	align_list *l=(align_list*)(R_ExternalPtrAddr(pRange));

	SEXP ans;
	PROTECT(ans=allocVector(INTSXP,7));
	INTEGER(ans)[0]=l->seqid;
	INTEGER(ans)[1]=l->range_begin;
	INTEGER(ans)[2]=l->range_end;
	INTEGER(ans)[3]=l->complex;
	INTEGER(ans)[4]=l->seq_LN;
	INTEGER(ans)[5]=l->min_seqlen;
	INTEGER(ans)[6]=l->max_seqlen;

	/* Names												*/
	SEXP names_vec;
	PROTECT(names_vec=allocVector(STRSXP,7));
	SET_STRING_ELT(names_vec,0,mkChar("seqid"));
	SET_STRING_ELT(names_vec,1,mkChar("qrBegin"));
	SET_STRING_ELT(names_vec,2,mkChar("qrEnd"));
	SET_STRING_ELT(names_vec,3,mkChar("complex"));
	SET_STRING_ELT(names_vec,4,mkChar("rSeqLen"));
	SET_STRING_ELT(names_vec,5,mkChar("qSeqMinLen"));
	SET_STRING_ELT(names_vec,6,mkChar("qSeqMaxLen"));
	setAttrib(ans,R_NamesSymbol,names_vec);

	UNPROTECT(2);
	return ans;
}

SEXP bam_range_get_refname(SEXP pRange)
{
	if(TYPEOF(pRange)!=EXTPTRSXP)
		error("[bam_range_get_coords] pRange is No external pointer!\n");
	align_list *l=(align_list*)(R_ExternalPtrAddr(pRange));

	if(!l->refname)
		return R_NilValue;

	SEXP ans;
	PROTECT(ans=allocVector(STRSXP,1));
	SET_STRING_ELT(ans,0,mkChar(l->refname));
	UNPROTECT(1);
	return ans;
}

SEXP bam_range_get_align_range(SEXP pRange)
{
	if(TYPEOF(pRange)!=EXTPTRSXP)
		error("[bam_range_get_coords] pRange is No external pointer!\n");

	align_list *l=(align_list*)(R_ExternalPtrAddr(pRange));
	if(!l->size)
		return R_NilValue;

	const bam1_t * align;
	unsigned long min_pos,max_end,i;
	min_pos=0;
	--(min_pos); /* Create big number as start value			*/
	max_end=0;
	uint32_t calend;

	align_element *curr_el=l->curr_el;
	wind_back(l);
	for(i=0;i<l->size;++i)
	{
		align=get_const_next_align(l);
		if(align->core.pos<min_pos)
			min_pos=align->core.pos;
		calend=bam_calend(&align->core,bam1_cigar(align));
		if(calend>max_end)
			max_end=calend;
	}
	l->curr_el=curr_el;

	SEXP ans;

	PROTECT(ans=allocVector(INTSXP,2));

	/* Returns minimum of position= 0-based position of leftmost match nucleotide	*/
	INTEGER(ans)[0]=min_pos;

	/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
	 * calend adds cigar shifts to position and thereby returns the
	 * 0-based position of the first nucleotide after the last match
	 * Returns the 0-based position of the rightmost match nucleotide
	 */
	INTEGER(ans)[1]=max_end-1;

	// + + + + + + + + + + + + + + + + + + + + + + + + + +
	// Names
	SEXP names_vec;
	PROTECT(names_vec=allocVector(STRSXP,2));

	SET_STRING_ELT(names_vec,0,mkChar("min_pos"));
	SET_STRING_ELT(names_vec,1,mkChar("max_end"));
	setAttrib(ans,R_NamesSymbol,names_vec);

	UNPROTECT(2);
	return ans;
}


SEXP bam_range_push_back(SEXP pRange,SEXP pAlign)
{
	if(TYPEOF(pRange)!=EXTPTRSXP)
		error("[bam_range_push_back] pRange is No external pointer!\n");
	if(TYPEOF(pAlign)!=EXTPTRSXP)
		error("[bam_range_push_back] pAlign is No external pointer!\n");
	align_list_push_back((align_list*)(R_ExternalPtrAddr(pRange)),(bam1_t*) (R_ExternalPtrAddr(pAlign)));
	return R_NilValue;
}

SEXP bam_range_push_front(SEXP pRange,SEXP pAlign)
{
	if(TYPEOF(pRange)!=EXTPTRSXP)
		error("[bam_range_push_front] pRange is No external pointer!\n");
	if(TYPEOF(pAlign)!=EXTPTRSXP)
		error("[bam_range_push_front] pAlign is No external pointer!\n");
	align_list_push_front((align_list*)(R_ExternalPtrAddr(pRange)),(bam1_t*) (R_ExternalPtrAddr(pAlign)));
	return R_NilValue;
}

SEXP bam_range_pop_back(SEXP pRange)
{
	if(TYPEOF(pRange)!=EXTPTRSXP)
		error("[bam_range_pop_back] pRange is No external pointer!\n");
	align_list_pop_back((align_list*)(R_ExternalPtrAddr(pRange)));
	return R_NilValue;
}
SEXP bam_range_pop_front(SEXP pRange)
{
	if(TYPEOF(pRange)!=EXTPTRSXP)
		error("[bam_range_pop_front] pRange is No external pointer!\n");
	align_list_pop_front((align_list*)(R_ExternalPtrAddr(pRange)));
	return R_NilValue;
}

SEXP bam_range_write_fastq(SEXP pRange,SEXP pFilename,SEXP pAppend)
{
	if(TYPEOF(pRange)!=EXTPTRSXP)
		error("[bam_range_write_fastq] No external pointer!");
	if(TYPEOF(pFilename)!=STRSXP)
		error("[bam_range_write_fastq] Filename must be a string!\n");
	if(TYPEOF(pAppend)!=LGLSXP)
		error("[bam_range_write_fastq] pAppend must be logical!\n");

	align_list *l=(align_list*)(R_ExternalPtrAddr(pRange));

	/* Read Adress of current align->for reconstitution at the end		*/
	align_element *e=l->curr_el;
	wind_back(l);
	const bam1_t *align;

	/* Open file for ouptut												*/
	_Bool append =*(LOGICAL(AS_LOGICAL(pAppend)));
	gzFile gz;
	if(append)
		gz =gzopen(CHAR(STRING_ELT(pFilename,0)),"ba");
	else
		gz =gzopen(CHAR(STRING_ELT(pFilename,0)),"bw");

	if(!gz)
		error("[bam_range_write_fastq] File could not be opened!");

	/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
	 * create buffer for sequence and qualities
	 */
	unsigned char *raw_seq;
	int32_t seq_len;
	int buf_size=2048;
	char *buf=(char*) calloc(buf_size,sizeof(char));
	uint8_t *quals;
	unsigned i=0,j;
	int nRows=(l->size);


	for(i=0;i<nRows;++i)
	{
		align=get_const_next_align(l);

		/* read sequence									*/
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

		/* write header line and sequence					*/
		gzprintf(gz,"@%s\n%s\n+\n",bam1_qname(align),buf);

		/* write quality string								*/
		quals=bam1_qual(align);
		for(j=0;j<seq_len;++j)
			buf[j]=(char) (quals[j]+33);
		buf[j]=0;
		gzprintf(gz,"%s\n",buf);

	}
	/* Reset curr_el pointer								*/
	l->curr_el=e;
	gzclose(gz);
	free(buf);

	//Rprintf("[bam_range_write_fastq] %u records written to file '%s'.\n",i,CHAR(STRING_ELT(pFilename,0)));
	return R_NilValue;
}

SEXP bam_range_write_fastq_index(SEXP pRange,SEXP pFilename,SEXP pWhichCopy,SEXP pAppend)
{
	if(TYPEOF(pRange)!=EXTPTRSXP)
		error("[bam_range_write_fastq_idx] No external pointer!");
	if(TYPEOF(pFilename)!=STRSXP)
		error("[bam_range_write_fastq_idx] Filename must be a string!\n");
	if(TYPEOF(pWhichCopy)!=INTSXP)
		error("[bam_range_write_fastq_idx] pWhichCopy must be integer!\n");
	if(TYPEOF(pAppend)!=LGLSXP)
		error("[bam_range_write_fastq_idx] pAppend must be logical!\n");

	/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
	 * The function reads aligns from range, checks pWichCopy and writes out to fastq file
	 */

	align_list *l=(align_list*)(R_ExternalPtrAddr(pRange));

	/* Read Adress of current align->for reconstitution at the end				*/
	align_element *e=l->curr_el;
	wind_back(l);
	const bam1_t *align;

	/* Open file for ouptut									*/
	_Bool append =*(LOGICAL(AS_LOGICAL(pAppend)));
	gzFile gz;
	if(append)
		gz=gzopen(CHAR(STRING_ELT(pFilename,0)),"ba");
	else
		gz=gzopen(CHAR(STRING_ELT(pFilename,0)),"bw");

	if(!gz)
		error("[bam_range_write_fastq_idx] File could not be opened!");

	/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
	 * create buffer for sequence and qualities
	 */
	unsigned char *raw_seq;
	int32_t seq_len;
	int buf_size=2048;
	char *buf=(char*) calloc(buf_size,sizeof(char));
	uint8_t *quals;
	unsigned j;

	unsigned i,nWritten=0;
	unsigned nCheck=INTEGER(pWhichCopy)[LENGTH(pWhichCopy)-1];

	nCheck= ((l->size)<nCheck) ? (l->size) : nCheck;
	//Rprintf("[bam_range_write_fastq_idx] Checking %u records.\n",nCheck);
	for(i=1;i<=nCheck;++i)
	{
		align=get_const_next_align(l);
		if(i==(INTEGER(pWhichCopy)[nWritten]))
		{
			/* read sequence									*/
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

			/* write header line and sequence					*/
			gzprintf(gz,"@%s\n%s\n+\n",bam1_qname(align),buf);

			/* write quality string								*/
			quals=bam1_qual(align);
			for(j=0;j<seq_len;++j)
				buf[j]=(char) (quals[j]+33);
			buf[j]=0;
			gzprintf(gz,"%s\n",buf);
			++nWritten;
		}
	}
	/* Reset curr_el pointer									*/
	gzclose(gz);
	free(buf);
	l->curr_el=e;
	//Rprintf("[bam_range_write_fastq_idx] %u records written to file '%s'.\n",nWritten,CHAR(STRING_ELT(pFilename,0)));

	SEXP ans;
	PROTECT(ans=allocVector(INTSXP,1));
	INTEGER(ans)[0]=nWritten;
	UNPROTECT(1);
	return ans;
}

SEXP bam_range_get_seqlen(SEXP pRange)
{
	if(TYPEOF(pRange)!=EXTPTRSXP)
		error("[bam_range_get_seqlen] No external pointer!");

	align_list *l=(align_list*)(R_ExternalPtrAddr(pRange));
	SEXP res;
	PROTECT(res=allocVector(INTSXP,2));
	if(l->size==0)
	{
		INTEGER(res)[0]=0;
		INTEGER(res)[1]=0;
	}
	else
	{
		INTEGER(res)[0]=l->min_seqlen;
		INTEGER(res)[1]=l->max_seqlen;
	}

	// + + + + + + + + + + + + + + + + + + + + + + + + + +
	// Names
	SEXP names_vec;
	PROTECT(names_vec=allocVector(STRSXP,2));

	SET_STRING_ELT(names_vec,0,mkChar("min"));
	SET_STRING_ELT(names_vec,1,mkChar("max"));
	setAttrib(res,R_NamesSymbol,names_vec);

	UNPROTECT(2);
	return res;
}

SEXP bam_range_get_qual_df(SEXP pRange)
{
	/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
	 * Count phred quality values of all aligns in range
	 * into data.frame
	 */

	if(TYPEOF(pRange)!=EXTPTRSXP)
		error("[bam_range_get_qual_df] No external pointer!");

	align_list *l=(align_list*)(R_ExternalPtrAddr(pRange));
	unsigned i,j,nCols=l->max_seqlen;

	/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
	 * SAM Format (p.3)	: QUAL=ASCII of Phred + 33
	 * Sanger FASTQ		: ASCII	33 - 126
	 * 					: phred  0 -  93 (94 values)
	 */
	unsigned nRows=94;
	unsigned char *quals;
	SEXP dflist, column_vector;
	PROTECT(dflist=allocVector(VECSXP,nCols));

	/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
	 * Column vectors: One for each position in sequence
	 */
	for(i=0;i<nCols;++i)
	{
		PROTECT(column_vector=allocVector(INTSXP,nRows));
		SET_VECTOR_ELT(dflist,i,column_vector);
		memset(INTEGER(column_vector),0,sizeof(int)*nRows);
	}

	align_element *e=l->curr_el;
	wind_back(l);
	for(i=0;i<l->size;++i)
	{
		const bam1_t *align=get_const_next_align(l);

		/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
		 * Count qual values for each
		 * position in 1 to seqlen
		 */
		quals=bam1_qual(align);
		for(j=0;j<(align->core.l_qseq);++j)
		{
			column_vector=VECTOR_ELT(dflist,j);
			// Counted values are truncated at nRows
			// because of static array size
			if(quals[j]<nRows)
				++(INTEGER(column_vector)[quals[j]]);
			else
				++(INTEGER(column_vector)[nRows-1]);
		}
	}
	// Restore curr_el pointer
	l->curr_el=e;

	/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
	 * Column Names
	 */
	SEXP col_names;
	PROTECT(col_names=allocVector(STRSXP,nCols));
	char *buf=R_alloc(buf_size,sizeof(char));
	// 1-based (as indexing of sequences in R)
	for(i=0;i<nCols;++i)
	{
    	sprintf(buf,"%i",i+1);
    	SET_STRING_ELT(col_names,i,mkChar(buf));
	}
	setAttrib(dflist,R_NamesSymbol,col_names);

	/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
	 *  Row Names
	 */
	SEXP row_names;
    PROTECT(row_names=allocVector(STRSXP,nRows));
    for(i=0;i<nRows;++i)
    {
    	sprintf(buf,"%i",i);
    	SET_STRING_ELT(row_names,i,mkChar(buf));
    }
    setAttrib(dflist,R_RowNamesSymbol,row_names);

	/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
	 * Return
	 */
	setAttrib(dflist,R_ClassSymbol,mkString("data.frame"));
	UNPROTECT(3+nCols);
	return dflist;
}

SEXP bam_range_get_align_depth(SEXP pRange,SEXP pGap)
{
	/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
	 * Count align depth for given range
	 * into REAL vector
	 */
	if(TYPEOF(pRange)!=EXTPTRSXP)
		error("[bam_range_get_align_depth] No external pointer!");
	if(TYPEOF(pGap)!=LGLSXP)
		error("[bam_range_get_align_depth] pGap must be LGLSXP!");


	unsigned nProtected=0;
	align_list *l=(align_list*)(R_ExternalPtrAddr(pRange));
	if(l->range_end==0)
		error("[bam_range_get_align_depth] range_end=0!");
	if(l->range_begin>=l->range_end)
		error("[bam_range_get_align_depth] end>=begin!");

	unsigned long begin,end,range_len,i,j;
	begin=l->range_begin;
	end=l->range_end;
	range_len=end-begin+1;

	unsigned long *c = (unsigned long*) calloc(range_len,sizeof(unsigned long));
	align_element *curr_el=l->curr_el;
	wind_back(l);

	if(LOGICAL(pGap)[0]==TRUE)
	{
		for(i=0;i<l->size;++i)
			count_align_gap_depth(c,begin,end,get_const_next_align(l));
	}
	else
	{
		for(i=0;i<l->size;++i)
			count_align_depth(c,begin,end,get_const_next_align(l));
	}
	/* Reset current pointer in align_list		*/
	l->curr_el=curr_el;

	/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
	 * alignDepth vector
	 */
	SEXP res;
	PROTECT(res=allocVector(INTSXP,range_len));
	++nProtected;
	for(i=0;i<range_len;++i)
	{
		INTEGER(res)[i]=c[i];
	}
	free(c);

	/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
	 * Position
	 * Directly passing res-INTEGER pointer:
	 * iteration for 0-initialization
	 * necessary
	 */

	SEXP pos_vector;
	PROTECT(pos_vector=allocVector(INTSXP,range_len));
	++nProtected;
    for(i=0,j=begin+1;i<range_len;++i,++j)
    	INTEGER(pos_vector)[i]=j;


	/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
	 * Create S4 object
	 */

    SEXP alignDepth,depth_name,pos_name,params,params_name,refname,refname_name;
    PROTECT(alignDepth=NEW_OBJECT(MAKE_CLASS("alignDepth")));
	++nProtected;

    PROTECT(depth_name=Rf_mkString("depth"));
	++nProtected;
    alignDepth=SET_SLOT(alignDepth,depth_name,res);

    PROTECT(pos_name=Rf_mkString("pos"));
	++nProtected;
    alignDepth=SET_SLOT(alignDepth,pos_name,pos_vector);

    PROTECT(params=allocVector(INTSXP,8));
	++nProtected;
    INTEGER(params)[0]=l->seqid;
    INTEGER(params)[1]=l->range_begin;
    INTEGER(params)[2]=l->range_end;
    INTEGER(params)[3]=l->complex;
    INTEGER(params)[4]=l->seq_LN;
    INTEGER(params)[5]=l->min_seqlen;
    INTEGER(params)[6]=l->max_seqlen;
    INTEGER(params)[7]=LOGICAL(pGap)[0]==TRUE ? 1:0;

	/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
	 *  Names for params
	 */

	SEXP names_vec;
	PROTECT(names_vec=allocVector(STRSXP,8));
	++nProtected;

	SET_STRING_ELT(names_vec,0,mkChar("seqid"));
	SET_STRING_ELT(names_vec,1,mkChar("qrBegin"));
	SET_STRING_ELT(names_vec,2,mkChar("qrEnd"));
	SET_STRING_ELT(names_vec,3,mkChar("complex"));
	SET_STRING_ELT(names_vec,4,mkChar("rSeqLen"));
	SET_STRING_ELT(names_vec,5,mkChar("qSeqMinLen"));
	SET_STRING_ELT(names_vec,6,mkChar("qSeqMaxLen"));
	SET_STRING_ELT(names_vec,7,mkChar("gap"));
	setAttrib(params,R_NamesSymbol,names_vec);

    PROTECT(params_name=Rf_mkString("params"));
	++nProtected;
    alignDepth=SET_SLOT(alignDepth,params_name,params);

    PROTECT(refname=allocVector(STRSXP,1));
    ++nProtected;
    if(!l->refname)
    	SET_STRING_ELT(refname,0,mkChar(""));
    else
    	SET_STRING_ELT(refname,0,mkChar(l->refname));

    PROTECT(refname_name=allocVector(STRSXP,1));
    ++nProtected;
    SET_STRING_ELT(refname_name,0,mkChar("refname"));
    alignDepth=SET_SLOT(alignDepth,refname_name,refname);

	UNPROTECT(nProtected);
	return alignDepth;
}

SEXP bam_range_count_nucs(SEXP pRange)
{
	if(TYPEOF(pRange)!=EXTPTRSXP)
		error("[bam_range_count_nucs] No external pointer!");

	align_list *l=(align_list*)(R_ExternalPtrAddr(pRange));

	/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
	 * Prepare counter
	 */
	const unsigned nCount=5;
	SEXP pCount=PROTECT(allocVector(INTSXP,nCount));
	int *count=INTEGER(pCount);
	memset(count,0,nCount*sizeof(int));

	/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
	 * Prepare list
	 */
	int32_t seq_len,i,j;
	align_element *e=l->curr_el;
	const bam1_t *align;
	unsigned char *raw_seq;

	/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
	 * Count nucleotides on align sequences
	 */
	wind_back(l);
	for(i=0;i<l->size;++i)
	{
		align=get_const_next_align(l);
		seq_len=align->core.l_qseq;
		raw_seq=bam1_seq(align);
		for(j=0;j<seq_len;++j)
			++count[bam_nt16_dna_table[bam1_seqi(raw_seq,j)]];
	}
	/* Restore curr_el pointer								*/
	l->curr_el=e;

	/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
	 * Names for pCount
	 */
	SEXP names_vec;
	PROTECT(names_vec=allocVector(STRSXP,nCount));

	SET_STRING_ELT(names_vec,0,mkChar("A"));
	SET_STRING_ELT(names_vec,1,mkChar("C"));
	SET_STRING_ELT(names_vec,2,mkChar("G"));
	SET_STRING_ELT(names_vec,3,mkChar("T"));
	SET_STRING_ELT(names_vec,4,mkChar("N"));
	setAttrib(pCount,R_NamesSymbol,names_vec);

	UNPROTECT(2);
	return pCount;
}

SEXP bam_range_idx_copy(SEXP pRange, SEXP pIndex)
{
	if(TYPEOF(pRange)!=EXTPTRSXP)
		error("[bam_range_idx_copy] No external pointer!");
	if(TYPEOF(pIndex)!=INTSXP)
		error("[bam_range_idx_copy] Index must be Integer!");

	// Expects:
	//   i) idx must be sorted in ascending order
	//  ii) idx[0] >0
	unsigned * idx = (unsigned*)INTEGER(pIndex);
	unsigned n_idx = (unsigned) LENGTH(pIndex);

	align_list *src, *res;
	src = (align_list*) (R_ExternalPtrAddr(pRange));

	if(src==NULL)
		error("[bam_range_idx_copy] align_list pointer=NULL!");

	res = align_list_idx_copy(src, idx, n_idx);


	/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
	 * Create bamRange S4 object
	 */
	SEXP bam_range;
    PROTECT(bam_range = NEW_OBJECT(MAKE_CLASS("bamRange")));

    SEXP list;
	PROTECT(list = R_MakeExternalPtr( (void*)(res), R_NilValue, R_NilValue));
	R_RegisterCFinalizer(list, finalize_bam_range);
	SEXP list_name;
	PROTECT(list_name = Rf_mkString("range"));
    bam_range = SET_SLOT(bam_range, list_name,list);

    /* Return													*/
	UNPROTECT(3);
	return bam_range;
}



/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 * bam_header
 */


static void finalize_bam_header(SEXP ptr)
{
	if(TYPEOF(ptr)!=EXTPTRSXP)
		error("[finalize_bam_header] No external pointer!");
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
		error("[init_bam_header] Header Text must be a string.\n");

	bam_header_t *h=(bam_header_t*)calloc(1, sizeof(bam_header_t));

	/* Copy header text														*/
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
		error("[bam_header_get_header_text] No external pointer!");

	bam_header_t* header=(bam_header_t*)(R_ExternalPtrAddr(pHeader));
	SEXP ans;
	PROTECT(ans=Rf_allocVector(STRSXP,1));
	SET_STRING_ELT(ans,0,mkChar(header->text));
	UNPROTECT(1);
	return ans;
}

/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 * bam_writer
 */

static void finalize_bam_writer(SEXP ptr)
{
	if(TYPEOF(ptr)!=EXTPTRSXP)
		error("[finalize_bam_writer] No external pointer!");
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
		error("[bam_writer_open] pHeader No external pointer!\n");
	if(TYPEOF(pFilename)!=STRSXP)
		error("[bam_writer_open] pFilename no string!\n");

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
		error("[bam_writer_open] pReader No external pointer!\n");
	if(TYPEOF(pFilename)!=STRSXP)
		error("[bam_writer_open] pFilename no string!\n");

	samfile_t *reader=(samfile_t*) (R_ExternalPtrAddr(pReader));
	samfile_t *writer=samopen(CHAR(STRING_ELT(pFilename,0)),"wb",reader->header);

	SEXP ptr;
	PROTECT(ptr=R_MakeExternalPtr( (void*) (writer),R_NilValue,R_NilValue));
	R_RegisterCFinalizer(ptr,finalize_bam_writer);
	UNPROTECT(1);
	return ptr;
}

SEXP bam_writer_save_align(SEXP pWriter,SEXP pAlign,SEXP pRefid)
{
	if(TYPEOF(pWriter)!=EXTPTRSXP)
		error("[bam_writer_save_align] No external pointer!\n");
	if(TYPEOF(pAlign)!=EXTPTRSXP)
		error("[bam_writer_save_align] No external pointer!\n");
	if(TYPEOF(pRefid)!=INTSXP)
		error("[bam_writer_save_align] pRefid must be integer!\n");

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
		error("[bam_writer_close] No exteranl pointer!\n");

	samfile_t *writer= (samfile_t*) (R_ExternalPtrAddr(pWriter));
	samclose(writer);
	R_SetExternalPtrAddr(pWriter,NULL);
	Rprintf("[bamWriter] closed.\n");
	return R_NilValue;
}

/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 * bam_align
 */


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
			error("[bam_align_getCigar_df] Cigar_type not in defined range!");

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
		error("[bam_align_get_segment_sequence] No external pointer!");

	/* Extract char* sequence with samtools						*/
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
		error("[bam_align_get_qualities] No external pointer!");
	bam1_t *align= (bam1_t*) (R_ExternalPtrAddr(pAlign));

	SEXP ans;
	PROTECT(ans=allocVector(STRSXP,1));
	SET_STRING_ELT(ans,0,mkChar((char*)bam1_qual(align)));
	UNPROTECT(1);
	return ans;
}

SEXP bam_align_get_qual_values(SEXP pAlign)
{
	if(TYPEOF(pAlign)!=EXTPTRSXP)
		error("[bam_align_get_qual_values] No external pointer!");
	bam1_t *align= (bam1_t*) (R_ExternalPtrAddr(pAlign));

	SEXP ans;
	int32_t seq_len=align->core.l_qseq;
	char *qual=(char*)bam1_qual(align);

	PROTECT(ans=allocVector(INTSXP,seq_len));
	int *ansval=INTEGER(ans);
	memset(ansval,0,seq_len*sizeof(int));

	unsigned i;
	for(i=0;i<seq_len;++i)
		INTEGER(ans)[i]=(int)qual[i];
	UNPROTECT(1);
	return ans;
}

SEXP bam_align_count_nucs(SEXP pAlign)
{
	if(TYPEOF(pAlign)!=EXTPTRSXP)
		error("[bam_align_count_nucs] No external pointer!");

	/* Prepare counter											*/
	const unsigned nCount=5;
	SEXP pCount=PROTECT(allocVector(INTSXP,nCount));
	int *count=INTEGER(pCount);
	memset(count,0,nCount*sizeof(int));

	/* Extract char* sequence with samtools						*/
	bam1_t *align= (bam1_t*) (R_ExternalPtrAddr(pAlign));
	int32_t seq_len=align->core.l_qseq;
	unsigned char *raw_seq=bam1_seq(align);
	int32_t i;
	for(i=0;i<seq_len;++i)
		++count[bam_nt16_dna_table[bam1_seqi(raw_seq,i)]];

	/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
	 * Names for pCount
	 */
	SEXP names_vec;
	PROTECT(names_vec=allocVector(STRSXP,nCount));

	SET_STRING_ELT(names_vec,0,mkChar("A"));
	SET_STRING_ELT(names_vec,1,mkChar("C"));
	SET_STRING_ELT(names_vec,2,mkChar("G"));
	SET_STRING_ELT(names_vec,3,mkChar("T"));
	SET_STRING_ELT(names_vec,4,mkChar("N"));
	setAttrib(pCount,R_NamesSymbol,names_vec);
	UNPROTECT(2);
	return pCount;
}

/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 * Alignment flags
 */

SEXP bam_align_is_paired(SEXP pAlign)
{
	if(TYPEOF(pAlign)!=EXTPTRSXP)
		error("[bam_align_is_paired] No external pointer!");
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
		error("[bam_align_mapped_in_proper_pair] No external pointer!");
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
		error("[bam_align_is_unmapped] No external pointer!");
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
		error("[bam_align_mate_is_unmapped] No external pointer!");
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
		error("[bam_align_strand_reverse] No external pointer!");
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
		error("[bam_align_strand_reverse] No external pointer!");

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
		error("[bam_align_is_first_in_pair] No external pointer!");
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
		error("[bam_align_is_second_in_pair] No external pointer!");
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
		error("[bam_align_is_secondary_align] No external pointer!");
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
		error("[bam_align_fail_qc] No external pointer!");

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
		error("[bam_align_is_pcr_or_optical_dup] No external pointer!");

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
		error("[bam_align_get_flag] No external pointer!");
	bam1_t *align=(bam1_t*)(R_ExternalPtrAddr(pAlign));

	SEXP ans;
	PROTECT(ans=allocVector(INTSXP,1));
	INTEGER(ans)[0]=align->core.flag;
	UNPROTECT(1);
	return ans;
}

/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 *  Writing accessors
 */

SEXP bam_align_set_refid(SEXP pAlign,SEXP pRefid)
{
	if(TYPEOF(pAlign)!=EXTPTRSXP)
		error("[bam_align_set_refid] No external pointer!");
	if(TYPEOF(pRefid)!=INTSXP)
		error("[bam_align_set_refid] pRefid must be integer!");
	bam1_t *align= (bam1_t*) (R_ExternalPtrAddr(pAlign));
	align->core.tid=INTEGER(pRefid)[0];

	return R_NilValue;
}

SEXP bam_align_set_is_paired(SEXP pAlign, SEXP val)
{
	if(TYPEOF(pAlign)!=EXTPTRSXP)
		error("[bam_align_set_is_paired] No external pointer!");
	if(TYPEOF(val)!=LGLSXP)
		error("[bam_align_set_is_paired] No bool value!");

	bam1_t *align=(bam1_t*)(R_ExternalPtrAddr(pAlign));
	set_flag(align,*(LOGICAL(AS_LOGICAL(val))),BAM_FPAIRED);
	return R_NilValue;
}

SEXP bam_align_set_mapped_in_proper_pair(SEXP pAlign, SEXP val)
{
	if(TYPEOF(pAlign)!=EXTPTRSXP)
		error("[bam_align_set_mapped_in_proper_pair] No external pointer!");
	if(TYPEOF(val)!=LGLSXP)
		error("[bam_align_set_mapped_in_proper_pair] No bool value!");

	bam1_t *align=(bam1_t*)(R_ExternalPtrAddr(pAlign));
	set_flag(align,*(LOGICAL(AS_LOGICAL(val))),BAM_FPROPER_PAIR);
	return R_NilValue;
}

SEXP bam_align_set_is_unmapped(SEXP pAlign, SEXP val)
{
	if(TYPEOF(pAlign)!=EXTPTRSXP)
		error("[bam_align_set_is_unmapped] No external pointer!");
	if(TYPEOF(val)!=LGLSXP)
		error("[bam_align_set_is_unmapped] No bool value!");

	bam1_t *align=(bam1_t*)(R_ExternalPtrAddr(pAlign));
	set_flag(align,*(LOGICAL(AS_LOGICAL(val))),BAM_FUNMAP);
	return R_NilValue;
}

SEXP bam_align_set_mate_is_unmapped(SEXP pAlign, SEXP val)
{
	if(TYPEOF(pAlign)!=EXTPTRSXP)
		error("[bam_align_set_mate_is_unmapped] No external pointer!");
	if(TYPEOF(val)!=LGLSXP)
		error("[bam_align_set_mate_is_unmapped] No bool value!");

	bam1_t *align=(bam1_t*)(R_ExternalPtrAddr(pAlign));
	set_flag(align,*(LOGICAL(AS_LOGICAL(val))),BAM_FMUNMAP);
	return R_NilValue;
}

SEXP bam_align_set_strand_reverse(SEXP pAlign, SEXP val)
{
	if(TYPEOF(pAlign)!=EXTPTRSXP)
		error("[bam_align_set_strand_reverse] No external pointer!");
	if(TYPEOF(val)!=LGLSXP)
		error("[bam_align_set_strand_reverse] No bool value!");

	bam1_t *align=(bam1_t*)(R_ExternalPtrAddr(pAlign));
	set_flag(align,*(LOGICAL(AS_LOGICAL(val))),16);
	return R_NilValue;
}

SEXP bam_align_set_mate_strand_reverse(SEXP pAlign, SEXP val)
{
	if(TYPEOF(pAlign)!=EXTPTRSXP)
		error("[bam_align_set_mate_strand_reverse] No external pointer!");
	if(TYPEOF(val)!=LGLSXP)
		error("[bam_align_set_mate_strand_reverse] No bool value!");

	bam1_t *align=(bam1_t*)(R_ExternalPtrAddr(pAlign));
	set_flag(align,*(LOGICAL(AS_LOGICAL(val))),32);
	return R_NilValue;
}

SEXP bam_align_set_is_first_in_pair(SEXP pAlign, SEXP val)
{
	if(TYPEOF(pAlign)!=EXTPTRSXP)
		error("[bam_align_set_is_first_in_pair] No external pointer!");
	if(TYPEOF(val)!=LGLSXP)
		error("[bam_align_set_is_first_in_pair] No bool value!");
	bam1_t *align=(bam1_t*)(R_ExternalPtrAddr(pAlign));
	set_flag(align,*(LOGICAL(AS_LOGICAL(val))),BAM_FREAD1);
	return R_NilValue;
}

SEXP bam_align_set_is_second_in_pair(SEXP pAlign, SEXP val)
{
	if(TYPEOF(pAlign)!=EXTPTRSXP)
		error("[bam_align_set_is_second_in_pair] No external pointer!");
	if(TYPEOF(val)!=LGLSXP)
		error("[bam_align_set_is_second_in_pair] No bool value!");

	bam1_t *align=(bam1_t*)(R_ExternalPtrAddr(pAlign));
	set_flag(align,*(LOGICAL(AS_LOGICAL(val))),BAM_FREAD2);
	return R_NilValue;
}

SEXP bam_align_set_is_secondary_align(SEXP pAlign, SEXP val)
{
	if(TYPEOF(pAlign)!=EXTPTRSXP)
		error("[bam_align_set_is_secondary_align] No external pointer!");
	if(TYPEOF(val)!=LGLSXP)
		error("[bam_align_set_is_secondary_align] No bool value!");

	bam1_t *align=(bam1_t*)(R_ExternalPtrAddr(pAlign));
	set_flag(align,*(LOGICAL(AS_LOGICAL(val))),BAM_FSECONDARY);
	return R_NilValue;
}

SEXP bam_align_set_fail_qc(SEXP pAlign, SEXP val)
{
	if(TYPEOF(pAlign)!=EXTPTRSXP)
		error("[bam_align_set_fail_qc] No external pointer!");
	if(TYPEOF(val)!=LGLSXP)
		error("[bam_align_set_fail_qc] No bool value!");

	bam1_t *align=(bam1_t*)(R_ExternalPtrAddr(pAlign));
	set_flag(align,*(LOGICAL(AS_LOGICAL(val))),BAM_FQCFAIL);
	return R_NilValue;
}

SEXP bam_align_set_is_pcr_or_optical_dup(SEXP pAlign, SEXP val)
{
	if(TYPEOF(pAlign)!=EXTPTRSXP)
		error("[bam_align_set_is_pcr_or_optical_dup] No external pointer!");
	if(TYPEOF(val)!=LGLSXP)
		error("[bam_align_set_is_pcr_or_optical_dup] No bool value!");

	bam1_t *align=(bam1_t*)(R_ExternalPtrAddr(pAlign));
	set_flag(align,*(LOGICAL(AS_LOGICAL(val))),BAM_FDUP);
	return R_NilValue;
}

SEXP bam_align_set_flag(SEXP pAlign, SEXP val)
{
	if(TYPEOF(pAlign)!=EXTPTRSXP)
		error("[bam_align_set_flag] No external pointer!");
	if(TYPEOF(val)!=INTSXP)
		error("[bam_align_set_flag] No integer value!");

	bam1_t *align=(bam1_t*)(R_ExternalPtrAddr(pAlign));
	align->core.flag=*INTEGER(val);
	return R_NilValue;
}

SEXP bam_align_create(SEXP pStrVals, SEXP pIntVals)
{

	/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
	 * String-values:
	 * 1) query-name
	 * 2) query sequence
	 * 3) quality string
	 * 4) CIGAR string
	 *
	 * Integer-values:
	 * 1) refid
	 * 2) position
	 * 3) flag
	 * 4) align quality
	 * 5) mate refid
	 * 6) mate position
	 * 7) insert size
	 */

	if(TYPEOF(pStrVals)!=STRSXP)
		error("[bam_align_create] pStrVals must be string value!");
	if(LENGTH(pStrVals)!=4)
		error("[bam_align_create] pStrVals must have length 4!");
	if(TYPEOF(pIntVals)!=INTSXP)
		error("[bam_align_create] pIntVals must be integer value!");
	if(LENGTH(pIntVals)!=7)
		error("[bam_align_create] pIntVals must have length 7!");


	const char* qname=CHAR(STRING_ELT(pStrVals,0));
	const char* seq=  CHAR(STRING_ELT(pStrVals,1));
	const char* qual= CHAR(STRING_ELT(pStrVals,2));
	const char* cigar=CHAR(STRING_ELT(pStrVals,3));

	int32_t refid       	=INTEGER(pIntVals)[0];
	int32_t position		=INTEGER(pIntVals)[1];
	long flag				=INTEGER(pIntVals)[2];
	uint32_t align_qual		=INTEGER(pIntVals)[3];
	// mate
	int32_t mate_refid		=INTEGER(pIntVals)[4];
	int32_t mate_position	=INTEGER(pIntVals)[5];

	int32_t insert_size		=INTEGER(pIntVals)[6];


	int doff=0,i;
	bam1_t *align=bam_init1();
	bam1_core_t *c = &align->core;
	uint8_t *p = 0;

#ifndef BAM1_ADD_CIGAR
	uint32_t *cigar_data;
#endif


	/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
	 * query name
	 */
	c->l_qname = strlen(qname) + 1;
	memcpy(alloc_align_data(align, doff + c->l_qname) + doff, qname, c->l_qname);
	doff += c->l_qname;


	/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
	 * Set atomic fields
	 */
	c->flag = flag;
	c->tid  = refid;
	c->pos  = position;
	c->qual = align_qual;

	// mate values
	c->mtid=mate_refid;
	c->mpos=mate_position;
	c->isize=insert_size;


	/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
	 * cigar
	 */
	{
		char *s, *t;
		int i, op;
		long x;
		c->n_cigar = 0;

		/* cigar[0]='*' checks for empty cigar				*/
		if (cigar[0] != '*')
		{
			/* check incoming string s and count n_cigar	*/
			for (s = (char*) cigar; *s; ++s)
			{
				/* *s must either be alpha or digit			*/
				if ((isalpha(*s)) || (*s=='='))
					++c->n_cigar;
				else if (!isdigit(*s))
					error("[bam_align_create] CIGAR character must be alphabetic or numeric or '='!");
			}

			align->data = alloc_align_data(align, doff + c->n_cigar * 4);

#ifdef BAM1_ADD_CIGAR
			align->cigar=calloc(c->n_cigar,sizeof(uint32_t));
#else
			cigar_data=calloc(c->n_cigar,sizeof(uint32_t));
#endif


			for (i=0, s= (char*)cigar; i != c->n_cigar; ++i)
			{
				/* convert string to integer							*/
				/* t will then point to first non-digit position behind	*/
				x = strtol(s, &t, 10);
				op = toupper(*t);
				if (op == 'M') op = BAM_CMATCH;
				else if (op == 'I') op = BAM_CINS;
				else if (op == 'D') op = BAM_CDEL;
				else if (op == 'N') op = BAM_CREF_SKIP;
				else if (op == 'S') op = BAM_CSOFT_CLIP;
				else if (op == 'H') op = BAM_CHARD_CLIP;
				else if (op == 'P') op = BAM_CPAD;
				else if (op == '=') op = BAM_CEQUAL;
				else if (op == 'X') op = BAM_CDIFF;
				else error("[bam_align_create] invalid CIGAR operation!");

				/* goto next								*/
				s = t + 1;

#ifdef BAM1_ADD_CIGAR
				bam1_cigar(align)[i]=(uint32_t) ((x << BAM_CIGAR_SHIFT) | ((long)op));
#else
				cigar_data[i] = (x << BAM_CIGAR_SHIFT) | op;
#endif


			}
			if(*s)
				error("[bam_align_create] Unmatched CIGAR operation");
			c->bin = bam_reg2bin(c->pos, bam_calend(c, bam1_cigar(align)));
			doff += c->n_cigar * 4;
		}
		else
		{
			/* s[0]='*' (empty cigar)						*/
			if (!(c->flag&BAM_FUNMAP)) {
				Rprintf("[bam_align_create] Parse warning: mapped sequence without CIGAR\n");
				c->flag |= BAM_FUNMAP;
			}
			c->bin = bam_reg2bin(c->pos, c->pos + 1);
		}
	}


	/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
	 * Sequence and qualities
	 */
	c->l_qseq = strlen(seq);

	if(strcmp(qual,"*"))
	{
		if (c->n_cigar && c->l_qseq != (int32_t)bam_cigar2qlen(c, bam1_cigar(align)))
			error("[bam_align_set_seq] cigar does not match sequence length.");
		p = (uint8_t*)alloc_align_data(align, doff + c->l_qseq + (c->l_qseq+1)/2) + doff;
		memset(p, 0, (c->l_qseq+1)/2);
		for (i = 0; i < c->l_qseq; ++i)
			p[i/2] |= bam_nt16_table[(int)seq[i]] << 4*(1-i%2);
	}
	else
		c->l_qseq = 0;

	/* qualities											*/
	if (strcmp(qual, "*") && c->l_qseq != strlen(qual))
				error("[bam_align_set_seq] Sequence and quality are inconsistent");
	p += (c->l_qseq+1)/2;
	if (strcmp(qual, "*") == 0)
		for (i = 0; i < c->l_qseq; ++i)
			p[i] = 0xff;
	else
		memcpy(p,qual,c->l_qseq+1);
	doff += c->l_qseq + (c->l_qseq+1)/2;


	/* No auxiliary data									*/
	align->l_aux=0;

	align->data_len=doff;

	SEXP ptr;
	PROTECT(ptr=R_MakeExternalPtr((void*)(align),R_NilValue,R_NilValue));
	R_RegisterCFinalizer(ptr,finalize_bam_align);
	UNPROTECT(1);
	return ptr;
}



/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 *
 * Miscellaneous functions
 *
 */


SEXP copy_fastq_records(SEXP pInfile,SEXP pOutfile,SEXP pWhichCopy,SEXP pAppend)
{
	if(TYPEOF(pInfile)!=STRSXP)
		error("[copy_fastq_records] pInfile must be a string!");
	if(TYPEOF(pOutfile)!=STRSXP)
		error("[copy_fastq_records] pOutfile must be a string!\n");
	if(TYPEOF(pWhichCopy)!=INTSXP)
		error("[copy_fastq_records] pWhichCopy must be integer!\n");
	if(TYPEOF(pAppend)!=LGLSXP)
		error("[copy_fastq_records] pAppend must be logical!\n");


	FILE *fin=fopen(CHAR(STRING_ELT(pInfile,0)),"r");
	if(fin==NULL)
		error("[copy_fastq_records] Infile does not exist!");

	/* Open file for ouptut									*/
	_Bool append =*(LOGICAL(AS_LOGICAL(pAppend)));
	FILE *fout;
	if(append)
		fout =fopen(CHAR(STRING_ELT(pOutfile,0)),"a+");
	else
		fout =fopen(CHAR(STRING_ELT(pOutfile,0)),"w");

	/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
	 * create buffer for sequence and qualities
	 */
	unsigned buf_size=2048;
	char *head,*seq,*qual,*buf;
	head=(char*) calloc(buf_size,sizeof(char));
	seq =(char*) calloc(buf_size,sizeof(char));
	qual=(char*) calloc(buf_size,sizeof(char));
	buf =(char*) calloc(buf_size,sizeof(char));

	unsigned nCheck=INTEGER(pWhichCopy)[LENGTH(pWhichCopy)-1];
	Rprintf("[copy_fastq_records] Checking %u records\n",nCheck);
	unsigned nWritten=0;
	rewind(fin);

	int c;
	unsigned i=0,ma=0,mp=0;
	while(!feof(fin) && i<nCheck)
	{
			c=fgetc(fin);
			if(c=='@')
			{
				/* Number of identified fastq items						*/
				++i;
				/* read head and seq									*/
				if(!fgets(head,buf_size,fin))
					break;
				if(!fgets(seq,buf_size,fin))
					break;
				c=fgetc(fin);
				if(c!='+')
				{
					++mp;
					/* skip following text until next @					*/
					while(!feof(fin) && c!='@')
						c=fgetc(fin);
					if(!feof(fin))
						ungetc('@',fin);
				}
				/* This should eat the '\n' after the '+'				*/
				c=fgetc(fin);
				if(!fgets(qual,buf_size,fin))
					break;

				/* check and eventually write item to outfile			*/
				if(i==(INTEGER(pWhichCopy)[nWritten]))
				{
					fprintf(fout,"@%s%s+\n%s",head,seq,qual);
					++nWritten;
				}
			}
			else
			{
				if(!feof(fin))
					++ma;

				/* skip following text until next @						*/
				while(!feof(fin) && c!='@')
					c=fgetc(fin);
				if(!feof(fin))
					ungetc('@',fin);
			}
	}
	Rprintf("[copy_fastq_records] Nr of '@' found: %u\tNr of missing '+': %u\tNr of missing '@': %u\n",i,mp,ma);

	fclose(fin);
	fclose(fout);
	free(head);
	free(seq);
	free(buf);
	free(qual);
	Rprintf("[copy_fastq_records] %u records written to file '%s'.\n",nWritten,CHAR(STRING_ELT(pOutfile,0)));

	SEXP ans;
	PROTECT(ans=allocVector(INTSXP,1));
	INTEGER(ans)[0]=nWritten;
	UNPROTECT(1);
	return ans;
}


SEXP count_fastq(SEXP pInfile,SEXP pMaxCol)
{
	/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
	 * Check incoming args
	 */

	if(TYPEOF(pInfile)!=STRSXP)
		error("[count_fastq] pInfile must be a string!");
	if(TYPEOF(pMaxCol)!=INTSXP)
		error("[count_fastq] pMaxCol must be Integer!");

	FILE *fin=fopen(CHAR(STRING_ELT(pInfile,0)),"r");
	if(fin==NULL)
		error("[count_fastq] Infile does not exist!");

	if(INTEGER(pMaxCol)[0]<1)
		error("[count_fastq] pMaxCol must be positive (>0)!");

	unsigned nProtected=0;


	/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
	 * Prepare vector for counting of nucleotide frequencies.
	 */

	const unsigned alpha_size  =nFastq+5;		/* Total number of tracked values							*/
	const unsigned iRead_pos   =nFastq;			/* Number of fastq items									*/
	const unsigned iMin_seqlen =nFastq+1;		/* Length of shortest sequence								*/
	const unsigned iMax_seqlen =nFastq+2;		/* Length of longest  sequence								*/
	const unsigned iMissedLines=nFastq+3;		/* Number of skipped lines when searched for next '@'		*/
	const unsigned iSeqQualLenUneq=nFastq+4;	/* Number of reads where Sequence and Qual differ in length	*/
												/* (or nCols is too small).									*/
	SEXP pAlpha=PROTECT(allocVector(INTSXP,alpha_size));
	++nProtected;
	int *alpha=INTEGER(pAlpha);
	memset(alpha,0,sizeof(unsigned)*alpha_size);
	/* Set to large value									*/
	alpha[iMin_seqlen]=30000;

	/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
	 * Names for pAlpha
	 */
	SEXP names_vec;
	PROTECT(names_vec=allocVector(STRSXP,alpha_size));
	++nProtected;

	SET_STRING_ELT(names_vec, 0,mkChar("A"));
	SET_STRING_ELT(names_vec, 1,mkChar("C"));
	SET_STRING_ELT(names_vec, 2,mkChar("G"));
	SET_STRING_ELT(names_vec, 3,mkChar("T"));
	SET_STRING_ELT(names_vec, 4,mkChar("M"));
	SET_STRING_ELT(names_vec, 5,mkChar("R"));
	SET_STRING_ELT(names_vec, 6,mkChar("W"));
	SET_STRING_ELT(names_vec, 7,mkChar("S"));
	SET_STRING_ELT(names_vec, 8,mkChar("Y"));
	SET_STRING_ELT(names_vec, 9,mkChar("K"));
	SET_STRING_ELT(names_vec,10,mkChar("V"));
	SET_STRING_ELT(names_vec,11,mkChar("H"));
	SET_STRING_ELT(names_vec,12,mkChar("D"));
	SET_STRING_ELT(names_vec,13,mkChar("B"));
	SET_STRING_ELT(names_vec,14,mkChar("N"));
	SET_STRING_ELT(names_vec,15,mkChar("+"));
	SET_STRING_ELT(names_vec,16,mkChar("-"));
	SET_STRING_ELT(names_vec,17,mkChar("="));
	SET_STRING_ELT(names_vec,18,mkChar("other"));
	SET_STRING_ELT(names_vec,19,mkChar("nReads"));
	SET_STRING_ELT(names_vec,20,mkChar("minSeqlen"));
	SET_STRING_ELT(names_vec,21,mkChar("maxSeqlen"));
	SET_STRING_ELT(names_vec,22,mkChar("missedLines"));
	SET_STRING_ELT(names_vec,23,mkChar("seqQualLenUneq"));
	setAttrib(pAlpha,R_NamesSymbol,names_vec);

	/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
	 *  Prepare data frame for counting of quality values
	 */
	int i,c,seqlen,qual_value,nRows,nCols;
	nCols=INTEGER(pMaxCol)[0];							/* Should be >= max seq-length				*/
	nRows=128;											/* Maximum phred quality score				*/

	SEXP dflist;
	PROTECT(dflist=allocVector(VECSXP,nCols));
	++nProtected;

	SEXP col_vec;
	for(i=0;i<nCols;++i)
	{
		col_vec=PROTECT(allocVector(INTSXP,nRows));
		++nProtected;
		memset(INTEGER(col_vec),0,nRows*sizeof(int));
		SET_VECTOR_ELT(dflist,i,col_vec);
	}

	/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
	 * Row names
	 */
	SEXP row_names;
    PROTECT(row_names=allocVector(STRSXP,nRows));
    ++nProtected;

    char *char_buf=(char*)calloc(buf_size,sizeof(char));
    for(i=0;i<nRows;++i)
    {
    	sprintf(char_buf,"%i",i);
    	SET_STRING_ELT(row_names,i,mkChar(char_buf));
    }
    setAttrib(dflist,R_RowNamesSymbol,row_names);

	/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
	 * Column Names
	 */
	SEXP col_names;
	PROTECT(col_names=allocVector(STRSXP,nCols));
	++nProtected;
    for(i=0;i<nCols;++i)
    {
    	sprintf(char_buf,"%i",i);
    	SET_STRING_ELT(col_names,i,mkChar(char_buf));
    }
	setAttrib(dflist,R_NamesSymbol,col_names);
	setAttrib(dflist,R_ClassSymbol,mkString("data.frame"));
	free(char_buf);

	/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
	 * Create buffer for sequence and qualities
	 */
	unsigned buf_size=2048;
	char *head,*seq,*qual,*iter;
	head=(char*) calloc(buf_size,sizeof(char));
	seq =(char*) calloc(buf_size,sizeof(char));
	qual=(char*) calloc(buf_size,sizeof(char));

	/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
	 * Collect data from fastq file
	 */
	while(!feof(fin))
	{
		c=fgetc(fin);
		if(c=='@')
		{
			/* Number of identified fastq items				*/
			++alpha[iRead_pos];
			/* read head and seq							*/
			if(!fgets(head,buf_size,fin))
				break;
			if(fgets(seq,buf_size,fin))
				break;

			/* count nucleotides in sequence				*/
			iter=seq;
			seqlen=0;
			while((*iter!='\0')&(*iter!='\n'))
			{
				++alpha[FASTQ_LETTERS[(unsigned char)*iter]];
				++iter;
				++seqlen;
			}

			if(seqlen>alpha[iMax_seqlen])
				alpha[iMax_seqlen]=seqlen;
			if(seqlen<alpha[iMin_seqlen])
				alpha[iMin_seqlen]=seqlen;

			/* Read '+' and Quality string					*/
			c=fgetc(fin);
			/* This should eat the '\n' after the '+'		*/
			c=fgetc(fin);
			if(!fgets(qual,buf_size,fin))
				break;

			iter=qual;
			i=0;
			while(((*iter)!='\n')&((*iter)!='\0')&(i<nCols))
			{
				col_vec=VECTOR_ELT(dflist,i);

				/* Calculate phred quality from fastq ASCII encoding	*/
				qual_value=(*iter)-33;
				if(qual_value<0)
					++(INTEGER(col_vec)[0]);
				else if(qual_value<nRows)
					++(INTEGER(col_vec)[qual_value]);
				else
					++(INTEGER(col_vec)[nRows-1]);

				++iter;
				++i;
			}
			/* Sequence and Quality have different length				*/
			/* (or nCols<seqlen)										*/
			if(i!=seqlen)
				++alpha[iSeqQualLenUneq];
		}
		else
		{
			// skip following text until next @
			while(!feof(fin) && c!='@')
			{
				++alpha[iMissedLines];
				c=fgetc(fin);
			}
			if(!feof(fin))
				ungetc('@',fin);
		}
	}
	fclose(fin);
	free(head);
	free(seq);
	free(qual);

	/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
	 * Prepare list as return container
	 */
	unsigned list_size=2;
	SEXP res_list=PROTECT(allocVector(VECSXP,list_size));
	++nProtected;

	SET_VECTOR_ELT(res_list,0,pAlpha);
	SET_VECTOR_ELT(res_list,1,dflist);

	SEXP list_names=PROTECT(allocVector(STRSXP,list_size));
	++nProtected;
	SET_STRING_ELT(list_names,0,mkChar("nuc_counts"));
	SET_STRING_ELT(list_names,1,mkChar("qual_counts"));
	setAttrib(res_list,R_NamesSymbol,list_names);
	setAttrib(res_list,R_ClassSymbol,mkString("list"));

	UNPROTECT(nProtected);
	return res_list;
}

/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 * Expects a vector of quantiles (pQuant) and
 * data.frame (pDf) where each column contains relative quantities (sums up to 1)
 *
 * Steps down each column and writes values for each quantile
 * into output data.frame
 */

SEXP get_col_quantiles(SEXP pQuant, SEXP pDf)
{
	/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
	 * pQuant is expected to contain unique
	 * and ascending sorted values
	 */
	if(TYPEOF(pQuant)!=REALSXP)
		error("[get_col_quantiles] pQuant must be REAL!");

	/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
	 * pDf is expected to be a data.frame
	 * where all columns contain relative values
	 * (between 0 and 1), so direct comparison to pQuant works
	 * and all columns sum up to 1.
	 */

	if(TYPEOF(pDf)!=VECSXP)
		error("[get_col_quantiles] pDf must be VECSXP!");

	/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
	 * Nr of Rows for writing = Nr of quantile values
	 */
	unsigned nwRows=LENGTH(pQuant);
	/* Nr of Rows for reading = Nr of rows in pDf			*/
	unsigned nrRows=LENGTH(VECTOR_ELT(pDf,0));
	/* Nr of Columns: equal in both data.frames				*/
	unsigned nCols=LENGTH(pDf);
	double val, *q=REAL(pQuant);
	unsigned i,j,k, nProtected=0;

	/* read and write qual-vectors							*/
	SEXP rQvec,wQvec;
	SEXP dflist;
	PROTECT(dflist=allocVector(VECSXP,nCols));
	++nProtected;

	for(i=0;i<nCols;++i)
	{
		PROTECT(wQvec=allocVector(INTSXP,nwRows));
		++nProtected;
		SET_VECTOR_ELT(dflist,i,wQvec);

		rQvec=VECTOR_ELT(pDf,i);
		if(TYPEOF(rQvec)!=REALSXP)
			error("[get_col_quantiles] All columns in pDf must be REAL!");

		val=0;
		k=0; /* which quantile we look for					*/
		for(j=0;(j<nrRows) & (k<nwRows);++j)
		{
			val+=REAL(rQvec)[j];
			if(val>q[k])
			{
				/* * * * * * * * * * * * * * * * * * *
				 * The k-th entry gets the row index
				 * where the cumulative value exeeds
				 * the k-th quantile
				 */
				INTEGER(wQvec)[k]=j;
				++k;
			}
		}
		/* Fill the rest of quantiles with max values		*/
		for(;k<nwRows;++k)
			INTEGER(wQvec)[k]=nrRows;
	}

	/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
	 * Column names
	 */
	SEXP col_names;
	PROTECT(col_names=allocVector(STRSXP,nCols));
	++nProtected;
	SEXP in_col_names=getAttrib(pDf,R_NamesSymbol);

	char *buf=(char*) calloc(buf_size,sizeof(char));
	for(i=0;i<nCols;++i)
	{
    	sprintf(buf,"%i",i);
    	SET_STRING_ELT(col_names,i,mkChar(CHAR(STRING_ELT(in_col_names,i))));
	}
	setAttrib(dflist,R_NamesSymbol,col_names);

	/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
	 * Row names
	 */
	SEXP row_names;
    PROTECT(row_names=allocVector(STRSXP,nwRows));
    ++nProtected;

    for(i=0;i<nwRows;++i)
    {
    	sprintf(buf,"q_%i",(int)(q[i]*100));
    	SET_STRING_ELT(row_names,i,mkChar(buf));
    }
    free(buf);
    setAttrib(dflist,R_RowNamesSymbol,row_names);
	setAttrib(dflist,R_ClassSymbol,mkString("data.frame"));

	UNPROTECT(nProtected);
	return(dflist);
}

SEXP count_text_lines(SEXP pInfile)
{
	if(TYPEOF(pInfile)!=STRSXP)
		error("[count_text_lines] pInfile must be a string!");

	FILE *fin=fopen(CHAR(STRING_ELT(pInfile,0)),"r");
	if(fin==NULL)
		error("[copy_fastq_records] Infile does not exist!");

	unsigned lines=0;
	char * buf = (char*) calloc(buf_size,sizeof(char));

	while(!feof(fin))
	{
		if(!fgets(buf,buf_size,fin))
			break;
		++lines;
	}
	fclose(fin);
	free(buf);
	--lines;

	SEXP res;
	PROTECT(res=allocVector(INTSXP,1));
	INTEGER(res)[0]=lines;
	UNPROTECT(1);
	return res;
}

/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 *
 * Definitions for R_registerRoutines
 *
 */


void R_init_rbamtools(DllInfo *info)
{
	R_CallMethodDef cmd[] ={
			{ "is_nil_externalptr", 				(DL_FUNC) &is_nil_externalptr,					1},

			/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
			 * bamHeader
			 */

			{ "init_bam_header",					(DL_FUNC) &init_bam_header,						1},
			{ "bam_header_get_header_text",			(DL_FUNC) &bam_header_get_header_text,			1},

			/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
			 * BamWriter
			 */
			{ "bam_writer_open",					(DL_FUNC) &bam_writer_open,						2},
			{ "bam_reader_open_writer",				(DL_FUNC) &bam_reader_open_writer,				2},
			{ "bam_writer_save_align",				(DL_FUNC) &bam_writer_save_align,				3},
			{ "bam_writer_close",					(DL_FUNC) &bam_writer_close,					1},

			/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
			 * BamReader
			 */
			{ "bam_reader_open",					(DL_FUNC) &bam_reader_open,						1},
			{ "bam_reader_close",					(DL_FUNC) &bam_reader_close,					1},
			{ "bam_reader_get_header_text",			(DL_FUNC) &bam_reader_get_header_text,			1},
			{ "bam_reader_get_ref_count",   		(DL_FUNC) &bam_reader_get_ref_count,    		1},
			{ "bam_reader_get_target_name",			(DL_FUNC) &bam_reader_get_target_name,  		2},
			{ "bam_reader_get_ref_data",			(DL_FUNC) &bam_reader_get_ref_data,				1},

			{ "bam_reader_create_index",			(DL_FUNC) &bam_reader_create_index,				2},
			{ "bam_reader_load_index",				(DL_FUNC) &bam_reader_load_index,				1},
			{ "bam_reader_unload_index",			(DL_FUNC) &bam_reader_unload_index,				1},
			{ "bam_reader_get_next_align",  		(DL_FUNC) &bam_reader_get_next_align,   		1},
			{ "bam_reader_save_aligns",				(DL_FUNC) &bam_reader_save_aligns,  			2},
			{ "bam_reader_sort_file",				(DL_FUNC) &bam_reader_sort_file,				4},

			{ "bam_reader_get_header",				(DL_FUNC) &bam_reader_get_header,				1},
			{ "bam_reader_tell",					(DL_FUNC) &bam_reader_tell,						1},
			{ "bam_reader_seek",					(DL_FUNC) &bam_reader_seek,						2},
			{ "bam_reader_write_fastq",  			(DL_FUNC) &bam_reader_write_fastq,   			3},
			{ "bam_reader_write_fastq_index",		(DL_FUNC) &bam_reader_write_fastq_index,  		4},
			{ "bam_count",							(DL_FUNC) &bam_count,							3},

			/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
			 * gap_list
			 */
			{ "create_gap_list",					(DL_FUNC) &create_gap_list,						0},
			{ "gap_list_fetch",						(DL_FUNC) &gap_list_fetch,						3},
			{ "gap_list_get_df",					(DL_FUNC) &gap_list_get_df,						1},
			{ "gap_list_get_size",  				(DL_FUNC) &gap_list_get_size,   				1},
			{ "gap_list_get_nAligns",				(DL_FUNC) &gap_list_get_nAligns,  				1},
			{ "gap_list_get_nAlignGaps",			(DL_FUNC) &gap_list_get_nAlignGaps,				1},

			/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
			 * gap_site_list
			 */
			{ "create_gap_site_list",				(DL_FUNC) &create_gap_site_list,				0},
			{ "gap_site_list_fetch",				(DL_FUNC) &gap_site_list_fetch,					3},
			{ "gap_site_list_get_df",				(DL_FUNC) &gap_site_list_get_df,				1},
			{ "gap_site_list_get_ref_id",  			(DL_FUNC) &gap_site_list_get_ref_id,   			1},
			{ "gap_site_list_get_size",				(DL_FUNC) &gap_site_list_get_size,  			1},
			{ "gap_site_list_get_nAligns",			(DL_FUNC) &gap_site_list_get_nAligns,			1},
			{ "gap_site_list_get_nAlignGaps",		(DL_FUNC) &gap_site_list_get_nAlignGaps,		1},
			{ "gap_site_list_merge",				(DL_FUNC) &gap_site_list_merge,					3},
			{ "gap_site_list_copy",					(DL_FUNC) &gap_site_list_copy,					1},
			{ "bitmask_r_zip",  					(DL_FUNC) &bitmask_r_zip,   					2},

			/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
			 * gap_site_list_list
			 */
			{ "gap_site_ll_init",					(DL_FUNC) &gap_site_ll_init,					0},
			{ "gap_site_ll_fetch",					(DL_FUNC) &gap_site_ll_fetch,					5},
			{ "gap_site_ll_get_df",					(DL_FUNC) &gap_site_ll_get_df,					2},
			{ "gap_site_ll_get_size",   			(DL_FUNC) &gap_site_ll_get_size,    			1},
			{ "gap_site_ll_get_nAligns",			(DL_FUNC) &gap_site_ll_get_nAligns,  			1},
			{ "gap_site_ll_get_nAlignGaps",			(DL_FUNC) &gap_site_ll_get_nAlignGaps,			1},
			{ "gap_site_ll_add_curr_pp",			(DL_FUNC) &gap_site_ll_add_curr_pp,				3},
			{ "gap_site_ll_add_merge_pp",			(DL_FUNC) &gap_site_ll_add_merge_pp,			4},
			{ "gap_site_ll_reset_refid",			(DL_FUNC) &gap_site_ll_reset_refid,				1},
			{ "gap_site_ll_get_summary_df",  		(DL_FUNC) &gap_site_ll_get_summary_df,   		1},
			{ "gap_site_ll_set_curr_first",			(DL_FUNC) &gap_site_ll_set_curr_first,  		1},

			/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
			 * bam_range
			 */
			{ "bam_range_init",						(DL_FUNC) &bam_range_init,						0},
			{ "bam_range_fetch",					(DL_FUNC) &bam_range_fetch,						4},
			{ "bam_range_get_size",					(DL_FUNC) &bam_range_get_size,					1},
			{ "bam_range_get_coords",   			(DL_FUNC) &bam_range_get_coords,    			1},
			{ "bam_range_get_params",				(DL_FUNC) &bam_range_get_params,  				1},
			{ "bam_range_get_refname",				(DL_FUNC) &bam_range_get_refname,				1},

			{ "bam_range_get_align_range",			(DL_FUNC) &bam_range_get_align_range,			1},
			{ "bam_range_get_next_align",			(DL_FUNC) &bam_range_get_next_align,			1},
			{ "bam_range_get_prev_align",			(DL_FUNC) &bam_range_get_prev_align,			1},
			{ "bam_range_step_next_align",			(DL_FUNC) &bam_range_step_next_align,  			1},
			{ "bam_range_step_prev_align",			(DL_FUNC) &bam_reader_sort_file,				1},

			{ "bam_range_get_align_df",				(DL_FUNC) &bam_range_get_align_df,				1},
			{ "bam_range_write",					(DL_FUNC) &bam_range_write,						3},
			{ "bam_range_wind_back",				(DL_FUNC) &bam_range_wind_back,					1},
			{ "bam_range_push_back",  				(DL_FUNC) &bam_range_push_back,   				2},
			{ "bam_range_pop_back",					(DL_FUNC) &bam_range_pop_back,  				1},
			{ "bam_range_push_front",				(DL_FUNC) &bam_range_push_front,				2},

			{ "bam_range_pop_front",				(DL_FUNC) &bam_range_pop_front,					1},
			{ "bam_range_write_current_align",		(DL_FUNC) &bam_range_write_current_align,		2},
			{ "bam_range_insert_past_curr_align",	(DL_FUNC) &bam_range_insert_past_curr_align,	2},
			{ "bam_range_insert_pre_curr_align",	(DL_FUNC) &bam_range_insert_pre_curr_align,    	2},
			{ "bam_range_mv_curr_align",			(DL_FUNC) &bam_range_mv_curr_align,  			2},
			{ "bam_range_write_fastq",				(DL_FUNC) &bam_range_write_fastq,				3},

			{ "bam_range_write_fastq_index",		(DL_FUNC) &bam_range_write_fastq_index,			4},
			{ "bam_range_get_seqlen",				(DL_FUNC) &bam_range_get_seqlen,				1},
			{ "bam_range_get_qual_df",				(DL_FUNC) &bam_range_get_qual_df,				1},
			{ "bam_range_get_align_depth",  		(DL_FUNC) &bam_range_get_align_depth,   		2},
			{ "bam_range_count_nucs",               (DL_FUNC) &bam_range_count_nucs,                1},
			{ "bam_range_idx_copy",               	(DL_FUNC) &bam_range_idx_copy,                	2},

			/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
			 * BamAlignment
			 */
			{ "bam_align_get_name",					(DL_FUNC) &bam_align_get_name,					1},
			{ "bam_align_get_refid",				(DL_FUNC) &bam_align_get_refid,					1},
			{ "bam_align_get_position",				(DL_FUNC) &bam_align_get_position,				1},
			{ "bam_align_get_nCigar",   			(DL_FUNC) &bam_align_get_nCigar,    			1},
			{ "bam_align_get_cigar_df",				(DL_FUNC) &bam_align_get_cigar_df,  			1},
			{ "bam_align_get_mate_refid",			(DL_FUNC) &bam_align_get_mate_refid,			1},

			{ "bam_align_get_mate_position",		(DL_FUNC) &bam_align_get_mate_position,			1},
			{ "bam_align_get_insert_size",			(DL_FUNC) &bam_align_get_insert_size,			1},
			{ "bam_align_get_map_quality",			(DL_FUNC) &bam_align_get_map_quality,			1},
			{ "bam_align_get_segment_sequence",  	(DL_FUNC) &bam_align_get_segment_sequence,   	1},
			{ "bam_align_get_qualities",			(DL_FUNC) &bam_align_get_qualities,  			1},
			{ "bam_align_get_qual_values",			(DL_FUNC) &bam_align_get_qual_values,			1},
			{ "bam_align_count_nucs",				(DL_FUNC) &bam_align_count_nucs,                1},


			/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
			 * alignment flags
			 */
			{ "bam_align_is_paired",				(DL_FUNC) &bam_align_is_paired,					1},
			{ "bam_align_mapped_in_proper_pair",	(DL_FUNC) &bam_align_mapped_in_proper_pair,		1},
			{ "bam_align_is_unmapped",				(DL_FUNC) &bam_align_is_unmapped,				1},
			{ "bam_align_mate_is_unmapped",			(DL_FUNC) &bam_align_mate_is_unmapped,			1},
			{ "bam_align_strand_reverse",  			(DL_FUNC) &bam_align_strand_reverse,   			1},
			{ "bam_align_mate_strand_reverse",		(DL_FUNC) &bam_align_mate_strand_reverse,  		1},
			{ "bam_align_is_first_in_pair",			(DL_FUNC) &bam_align_is_first_in_pair,			1},
			{ "bam_align_is_second_in_pair",		(DL_FUNC) &bam_align_is_second_in_pair,			1},
			{ "bam_align_is_secondary_align",		(DL_FUNC) &bam_align_is_secondary_align,		1},
			{ "bam_align_fail_qc",  				(DL_FUNC) &bam_align_fail_qc,   				1},
			{ "bam_align_is_pcr_or_optical_dup",	(DL_FUNC) &bam_align_is_pcr_or_optical_dup,  	1},
			{ "bam_align_get_flag",					(DL_FUNC) &bam_align_get_flag,					1},

			{ "bam_align_set_refid",				(DL_FUNC) &bam_align_set_refid,					2},
			{ "bam_align_set_is_paired",			(DL_FUNC) &bam_align_set_is_paired,				2},
			{ "bam_align_set_mapped_in_proper_pair",(DL_FUNC) &bam_align_set_mapped_in_proper_pair,	2},
			{ "bam_align_set_is_unmapped",   		(DL_FUNC) &bam_align_set_is_unmapped,    		2},
			{ "bam_align_set_mate_is_unmapped",		(DL_FUNC) &bam_align_set_mate_is_unmapped,  	2},
			{ "bam_align_set_strand_reverse",		(DL_FUNC) &bam_align_set_strand_reverse,		2},
			{ "bam_align_set_mate_strand_reverse",	(DL_FUNC) &bam_align_set_mate_strand_reverse,	2},
			{ "bam_align_set_is_first_in_pair",		(DL_FUNC) &bam_align_set_is_first_in_pair,		2},
			{ "bam_align_set_is_second_in_pair",	(DL_FUNC) &bam_align_set_is_second_in_pair,		2},
			{ "bam_align_set_is_secondary_align",  	(DL_FUNC) &bam_align_set_is_secondary_align,   	2},
			{ "bam_align_set_fail_qc",				(DL_FUNC) &bam_align_set_fail_qc,				2},
			{ "bam_align_set_is_pcr_or_optical_dup",(DL_FUNC) &bam_align_set_is_pcr_or_optical_dup,	2},
			{ "bam_align_set_flag",					(DL_FUNC) &bam_align_set_flag,					2},
			{ "bam_align_create",  					(DL_FUNC) &bam_align_create,   					2},

			/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
			 * Miscellaneous functions
			 */
			{ "copy_fastq_records",					(DL_FUNC) &copy_fastq_records,					4},
			{ "count_fastq",						(DL_FUNC) &count_fastq,							2},
			{ "get_col_quantiles",					(DL_FUNC) &get_col_quantiles,					2},
			{ "count_text_lines",					(DL_FUNC) &count_text_lines,					1},

			{NULL, NULL, 0}
	};
	/*
	 * 			{ "",	(DL_FUNC) &,	}
	 */
	R_registerRoutines(info, NULL, cmd, NULL, NULL);
}


#endif
