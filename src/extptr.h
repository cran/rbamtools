/*
 * exptr.hpp
 *
 *  Created on: 30.01.2015
 *      Author: kaisers
 */

#ifndef EXPTR_HPP_
#define EXPTR_HPP_

#include<Rdefines.h>
#include<Rinternals.h>
#include<memory>

template<typename T>
static void _finalizer(SEXP ext)
{
	if(TYPEOF(ext)==EXTPTRSXP)
	{
		shared_ptr<T> *p = (shared_ptr<T>*) R_ExternalPtrAddr(ext);
		//if(p)
		//	Rprintf("[_finalizer]: ptr %p\n", p);
		delete p;
	}
}



template<typename T>
class extptr
{
public:

	// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - //
	// Constructors
	// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - //
	explicit extptr(): sp_(new shared_ptr<T>(make_shared<T>())){
		//Rprintf("[extptr] extptr() id: %3i\n", (*sp_)->getValue());
	}
	// Copies incoming object
	extptr(T &p): sp_(new shared_ptr<T>(make_shared<T>())) { **sp_ = p;	}
	//extptr(shared_ptr<T> &p): sp_(new shared_ptr<T>(make_shared<T>())) { *sp_= p; }

	explicit extptr(SEXP pPtr)//: sp_(new shared_ptr<T>(make_shared<T>())) // ToDo: Throw exception ??
	{
		if(TYPEOF(pPtr) != EXTPTRSXP)
			error("[extptr] No external pointer!");

		if(!R_ExternalPtrAddr(pPtr))
			error("[extptr] Received Nil pointer!");

		shared_ptr<T> * sp = (shared_ptr<T> *) (R_ExternalPtrAddr(pPtr));
		sp_ = new shared_ptr<T>(*sp);

	}

	virtual ~extptr()
	{
		//Rprintf("[extptr] ~exptr ptr : %p\n", sp_);
		delete sp_;
	}

	// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - //
	// Operator overloading
	// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - //
    extptr& operator=(const extptr& rhs)
    {
    	if(this != &rhs)
    		*sp_ = *(rhs.sp_); // Share pointer
   	    return *this ;
    }

    T& operator*() const { return *(*sp_); }
    T* operator->() const { return sp_->get(); }

	operator SEXP() const
	{
		shared_ptr<T> * p = new shared_ptr<T>;
		// Shared copy
		*p = *sp_;
		//Rprintf("[extptr] operator SEXP: %p\n", p);

	    SEXP ext = PROTECT(R_MakeExternalPtr(p, R_NilValue, R_NilValue));
	    R_RegisterCFinalizerEx(ext, _finalizer<T>, TRUE);
	    UNPROTECT(1);
	    return ext;
	}

private:
	shared_ptr<T> * sp_;
};



template<typename T>
SEXP to_string(const extptr<T> &p)
{
	SEXP pRes  = PROTECT(allocVector(STRSXP, 1));
	SET_STRING_ELT(pRes, 0, mkChar(string(*p).c_str()));
	UNPROTECT(1);
	return pRes;
}


template<typename T>
class atmptr
{
	explicit atmptr(unsigned size): protected_(false), pRob(R_NilValue), p_(0) { }
	virtual ~atmptr() { if(protected_) UNPROTECT(pRob); }

    T& operator*() const { return *p_; }
    T* operator->() const { return p_; }
    T& operator [] (int i) { return p_[i]; }
    int length() const { return length(pRob); }

	operator SEXP() const { return pRob; }

private:
	bool protected_;
	SEXP pRob;
	T * p_;
};

template<>
class atmptr<int>
{
public:
	explicit atmptr(unsigned size): protected_(true), pRob(R_NilValue), p_(0)
	{
		pRob = PROTECT(allocVector(INTSXP, size));
		p_ = INTEGER(pRob);
	}

	explicit atmptr(SEXP p): protected_(false)
	{
		if(TYPEOF(p) != INTSXP)
			error("atmptr<int> SEXP must be of type INTSXP");

		pRob=p;
		p_ = INTEGER(pRob);
	}

	virtual ~atmptr()
	{
		if(protected_)
			UNPROTECT(1);
	}

    int& operator*() const { return *p_; }
    int* operator->() const { return p_; }
    int& operator [] (int i) { return p_[i]; }
    int size() const { return length(pRob); }

	operator SEXP() const { return pRob; }

private:
	bool protected_;
	SEXP pRob;
	int * p_;
};


template<>
class atmptr<char>
{
public:
	explicit atmptr(unsigned size): protected_(true), pRob(R_NilValue)
	{
		pRob = PROTECT(allocVector(STRSXP, size));
	}

	explicit atmptr(SEXP p): protected_(false)
	{
		if(TYPEOF(p) != STRSXP)
			error("atmptr<char> SEXP must be of type STRSXP");

		pRob=p;
	}

	explicit atmptr(const vector<string> &v): protected_(true)
	{
		int i, n = v.size();
		pRob = PROTECT(allocVector(STRSXP, n));
		for(i = 0; i < n; ++i)
			SET_STRING_ELT(pRob, i, mkChar(v[i].c_str()));
	}

	virtual ~atmptr()
	{
		if(protected_)
			UNPROTECT(1);
	}

	void get(unsigned i, string &s)
	{
		s.clear();
		if(i < (unsigned) LENGTH(pRob))
			s = CHAR(STRING_ELT(pRob, i));
	}

	void set(unsigned i, string &s)
	{
		if(i < (unsigned) LENGTH(pRob))
			SET_STRING_ELT(pRob, i, mkChar(s.c_str()));
	}

    int length() const { return LENGTH(pRob); }

	operator SEXP() const { return pRob; }

private:
	bool protected_;
	SEXP pRob;
};



SEXP to_string(const string & s)
{
	SEXP pRes = PROTECT(allocVector(STRSXP, 1));
	SET_STRING_ELT(pRes, 0, mkChar(s.c_str()));
	UNPROTECT(1);
	return pRes;
}


#endif /* EXPTR_HPP_ */
