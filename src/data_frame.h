/*
 * data_frame.h
 *
 *  Created on: 10.02.2015
 *      Author: wolfgang
 */

#ifndef DATA_FRAME_H_
#define DATA_FRAME_H_

#include <Rdefines.h>
#include <sstream>
//#include <cstdlib>

class data_frame
{
public:
	data_frame(unsigned nrow, unsigned ncol): nrow_(nrow), ncol_(ncol), next_column(0)
	{
		PROTECT(dflist=allocVector(VECSXP, ncol_));
		PROTECT(col_names=allocVector(STRSXP, ncol_));
		PROTECT(row_names=allocVector(STRSXP, nrow_));

		fill_row_names();
		setAttrib(dflist,R_ClassSymbol,mkString("data.frame"));
	    setAttrib(dflist,R_RowNamesSymbol,row_names);
		setAttrib(dflist,R_NamesSymbol,col_names);
	}

	~data_frame() { UNPROTECT(3 + next_column); }

	int addColumn(SEXP pCol, const char* colname)
	{
		if( ((unsigned) length(pCol)) != nrow_)
			error("[data_frame] Length of column vector must be equal to number of rows!");

		if(next_column < ncol_)
		{
			SEXP p;
			PROTECT(p=pCol);
			SET_VECTOR_ELT(dflist, next_column, p);
			SET_STRING_ELT(col_names, next_column, mkChar(colname));
			return ++next_column;
		}
		return -1;
	}

	operator SEXP() const { return dflist; }

	unsigned nrow() const { return nrow_; }
	unsigned ncol() const { return ncol_; }

	int * addIntColumn(const char * colname)
	{
		if(next_column < ncol_)
		{
			SEXP p = PROTECT(allocVector(INTSXP, nrow_));
			SET_VECTOR_ELT(dflist, next_column, p);
			SET_STRING_ELT(col_names, next_column, mkChar(colname));
			++next_column;
			return INTEGER(p);
		}
		return 0;
	}

	double * addRealColumn(const char * colname)
	{
		if(next_column < ncol_)
		{
			SEXP p = PROTECT(allocVector(REALSXP, nrow_));
			SET_VECTOR_ELT(dflist, next_column, p);
			SET_STRING_ELT(col_names, next_column, mkChar(colname));
			++next_column;
			return REAL(p);
		}
		return 0;
	}

private:
	unsigned nrow_;
	unsigned ncol_;
	unsigned next_column;

	SEXP dflist;
	SEXP col_names;
	SEXP row_names;

private:
	void fill_row_names()
	{
		unsigned i;
		std::stringstream sst;
	    for(i=0; i < nrow_; ++i)
	    {
	    	sst << i + 1;
	    	SET_STRING_ELT(row_names, i, mkChar(sst.str().c_str()));
	    	sst.str(std::string());
	    	sst.clear();
	    }
	}
};


#endif /* DATA_FRAME_H_ */
