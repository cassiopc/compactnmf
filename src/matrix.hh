/*
## matrix.hh (package compactnmf) Copyright (C) 2012 C. P. de Campos (cassiopc@acm.org)

## This program is free software: you can redistribute it and/or modify
## it under the terms of the GNU General Public License as published by
## the Free Software Foundation, either version 3 of the License, or
## (at your option) any later version.
## This program is distributed in the hope that it will be useful,
## but WITHOUT ANY WARRANTY; without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
## GNU General Public License for more details.
## A copy of the GNU General Public License is available at
## <http://www.gnu.org/licenses/>.

## I appreciate very much if any derivate publication cite my work.
## Please send me an email to ask for the proper citation to be
## included in a work that uses the code available here.
*/

#ifndef MATRIX_HH
#define MATRIX_HH

#ifdef USE_CBLAS
extern "C" {
#include <gsl/gsl_cblas.h>
}
#endif
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include "withr.hh"
#include "general.hh"

typedef long tint;

#include<iostream>
#include<string>
#include<map>
#include<vector>
#include<algorithm>
#include<cmath>

// template matrix of elements with 1 or 2 dimensions
template <class TT=double>
class matrix {
  // dimensions of the matrix
  tint n;
  tint m;
  // desaloc is a flag to indicate if the data shall be desaloccated when the matrix is destroyed
  bool desaloc;
  // data actually contains the data of the matrix
  TT *data;
  // column names
  std::vector<std::string> colnames;
  std::map<std::string,tint> mapcolnames;
  // line names
  std::vector<std::string> linenames;
  std::map<std::string,tint> maplinenames;

	// function to aloc the space necessary to keep n lines and m columns
	void aloc(int tozero=0, double zero=0);
	// copy the content of the other matrix to the current matrix, but does not
	// allocate memory for that. So the memory must be already available.
	void copy(const matrix &other);
	// auxiliary pointer and function to sort the columns of the matrix  
	static matrix<TT> *sortcolstemp;
	static int colcmp(const void *a, const void *b);

public:
#ifndef CPURE
	static SEXP makeobj(std::vector<std::string> const &names, std::vector<matrix<double> > const &mats) {
		int n = names.size();
		SEXP rl = PROTECT(allocVector(VECSXP,n));
		SEXP nm = PROTECT(allocVector(STRSXP,n));
		for(int k=0; k<n; k++) {
			SEXP w;
			PROTECT(w = allocMatrix(REALSXP, mats[k].N(), mats[k].M()));
			SEXP dimv = getAttrib(w, R_DimSymbol);
			for(int i=0; i<mats[k].N(); i++)
				for(int j=0; j<mats[k].M(); j++)
					getp(w,dimv,i,j) = mats[k](i,j);
			
			SET_VECTOR_ELT(rl, k, w);
			SET_STRING_ELT(nm, k, mkChar(names[k].c_str()));
		}
		setAttrib(rl, R_NamesSymbol, nm);
		UNPROTECT(2+n);
		return rl;
	}
	
	static SEXP mat2sexp(const matrix &mat) {
		DEBMACRO(Rprintf("mat2sexp %d %d\n", mat.N(), mat.M());)
			SEXP v;
		PROTECT(v = allocMatrix(REALSXP, mat.N(), mat.M()));
		SEXP dimv = getAttrib(v, R_DimSymbol);
		for(int i=0; i<mat.N(); i++)
			for(int j=0; j<mat.M(); j++)
				getp(v,dimv,i,j) = mat(i,j);
		
		UNPROTECT(1);
		DEBMACRO(Rprintf("mat2sexp end\n");)
			return v;
	}
#endif

	tint getcolumn(const std::string &a) const {
	  std::map<std::string,tint>::iterator i = mapcolnames.find(a);
		if(i == mapcolnames.end()) return -1;
		return i->second;
	}
	tint getline(const std::string &a) const {
	  std::map<std::string,tint>::iterator i = maplinenames.find(a);
		if(i == maplinenames.end()) return -1;
		return i->second;
	}

	std::string getname(tint i=-1, tint j=-1) {
		if(i>=0) {
			return linenames[i];
		} else if(j<0) {
			error("Invalid call to getname");
		}
		return colnames[j];
	}
	void setname(tint i=-1, tint j=-1, const std::string &a="") {
		if(i>=0) {
			linenames[i] = a;
			maplinenames.insert(std::make_pair<std::string,tint>(a,i));
		}
		if(j>=0) {
			colnames[j] = a;
			mapcolnames.insert(std::make_pair<std::string,tint>(a,j));
		}
	}

  static double kbs;
  static double megas();

  double euclidean(const matrix<TT> &other) const;
  double correlation(const matrix<TT> &other) const;
  double l1distance(const matrix<TT> &other) const;
  void maxpos(int &a, int &b) const;

#ifndef CPURE
  // set the data of the current matrix to be the data that came from R
  // note that no extra data is allocated, but the same data space
  // of R is used. This avoid duplication of a lot of data, but changes
  // will affect the data that came from R (this is a feature!)
  void set(SEXP s) {
	  SEXP dim = getAttrib(s, R_DimSymbol);
	  if(length(dim) > 0) {
		  n = INTEGER(dim)[0];
		  if(length(dim)>1)
			  m = INTEGER(dim)[1];
		  else m = 1;
		  desaloc = false;
		  data = REAL(s);
	  } else {
		  n=m=0;
		  data=0;
	  }
  }
  // gets a matrix in R format and returns a matrix
  static matrix fromsexp(SEXP t) {
    matrix mat(t);
    return mat;
  }
  // constructor from an R matrix
  // SEXP is already all we need, so we just point data to the current data of SEXP
  // and set the sizes accordingly. The only care here is to not desalloc memory
  // when destroying the object, because it was not allocated by me
  matrix(SEXP s): desaloc(false) {
    set(s);
  }
#endif

  void fill(TT zero=0);
  tint njumps() const;
	matrix<int> jumps() const;
  bool isvector() const;
  bool iscolumn() const { return size()>0 && m==1; }
  matrix<int> toint() const;
  matrix<long> tolong() const;
  matrix<char> tochar() const;
  matrix<double> todouble() const;
  matrix distinct() const;
	matrix<int> tocat_bycol() const;

/*
	void sort(int (*compar)(const void *, const void *) const = &(matrix::cmp));
	int cmp(const void *a, const void *b) const {
		TT d = *(TT *)a - *(TT *)b;
		return d<0? -1: (d>0? 1: 0);
	}
*/

	void sort() {
		tint n=size();
		std::vector<TT> v(n);
		for(tint i=0; i<n; i++) v[i] = data[i];
		std::sort(v.begin(),v.end());
		for(tint i=0; i<n; i++) data[i] = v[i];
	}

//	void sortint();

  // read-only access to n
  inline tint N() const { return n; }
  // read-only access to m
  inline tint M() const { return m; }

  // returns an integer array of the indices of the columns ordered by
  // their elements. Matrices with multiple lines will have the elements
  // ordered using all lines, comparing first line, then second, and so on.
  // TThe matrix is not altered and only an array of indices corresponding
  // to the desired ordering is returned.
  matrix<double> sortcols();

  static matrix tomatrix(int *v, tint n);
  static matrix tomatrix(long *v, tint n);
  static matrix tomatrix(double *v, tint n);

  // return true iff this matrix doesnt contain the same data as the other
  // and of course have the same size
  bool operator!=(const matrix &other) const;
  bool operator==(const matrix &other) const;
  matrix operator>(TT val) const;
  bool operator<(const matrix &other) const;

  double entropy() const;
  double condentropy() const;
  double chisquared() const;
	double mutualinfo(const matrix &other) const;

//	matrix filter(

  // writes the binary data of the matrix to the file, that must be
  // already open and with writing permission.  
  void tofile(FILE *fp) const;
  // reads from the file fp binary data to populate the current matrix.
  // It desaloc the current data e alocs space for the new data.
  void fromfile(FILE *fp);
  // writes an array of matrices to a file. fn must contain the filename,
  // lis contains nmat elements, each one is a matrix. TThe number of
  // matrices being written is also put in the file.
  static void tofile(const char *fn, const matrix *lis, int nmat);
  // reads an array of matrices from a file. fn is the filename, and
  // nmat is used to return the number of elements that is read.
  static matrix<TT> *fromfile(const char *fn, int *nmat);
    
  // default constructor creates an empty matrix  
  matrix(tint n=0, tint m=1): n(n), m(m), desaloc(true) {
    aloc();
  }

  // default constructor creates an empty matrix  
  matrix(tint n, tint m, double zero): n(n), m(m), desaloc(true) {
    aloc(1,zero);
  }

  // copy constructor
  matrix(const matrix &other): n(other.n), m(other.m), desaloc(true) {
    aloc();
    copy(other);      
  }
  inline bool empty() const { return size()==0 || !data; }
  
	void print(std::string s=std::string(""), tint x=-1, tint y=-1, bool typeint=false, FILE *saida=stderr) const;
	void tocsv(std::string *s=0, FILE *saida=stderr) const;
  static matrix<TT> fromcsv(FILE *ent=stdin, tint n=0, tint m=0, bool lineheader=false, bool colheader=false);

	TT operation(int type, tint i=-1, tint j=-1) const;
	matrix<TT> apply(int dim, int type) const;
    
  // returns the minimum element of column j
  TT min(tint i=-1, tint j=-1) const;
  // returns the maximum element of column j
  TT max(tint i=-1, tint j=-1) const;
  TT findmax(tint i=-1, tint j=-1) const;
  // returns the sum of the matrix elements
  TT sum(tint i=-1, tint j=-1) const;
  // returns the absolute maximum value of the matrix elements
  double absmax(tint i=-1, tint j=-1) const;
  // returns the mean of the matrix elements
  double mean(tint i=-1, tint j=-1) const;
  // returns the mean of squared elements of the matrix
  double meansquare(tint i=-1, tint j=-1) const;
  // returns the variance of the elements in the matrix
  double var(tint i=-1, tint j=-1) const {
    double d = mean(i,j);
    return meansquare(i,j) - d*d;
  }
  // returns the standard deviation of the elements
  double sd(tint i=-1, tint j=-1) const {
	  return std::sqrt(var(i,j));
  }

  // operator to access the elements of the matrix. This is the
  // updatable version
  inline TT &operator()(tint i, tint j=0) {
#ifndef OPTIMIZE
    if(n==1 && j==0 && i>=0 && i<m) return data[i];
    if(i<0 || j<0 || i>=n || j>=m || !data)
      error("error in () i=%ld j=%ld n=%ld m=%ld\n", i,j,n,m);
#endif
    return data[i+j*n];
  }
  // operator to access the elements of the matrix. This is the
  // non-updatable version to be used in const matrices
  inline const TT &operator()(tint i, tint j=0) const {
#ifndef OPTIMIZE
    if(n==1 && j==0 && i>=0 && i<m) return data[i];
    if(i<0 || j<0 || i>=n || j>=m || !data)
      error("error in () i=%ld j=%ld n=%ld m=%ld\n", i,j,n,m);
#endif
    return data[i+j*n];
  }
  // access the elements of the matrix as integers.
  inline const int integer(tint i, tint j=0) const {
#ifndef OPTIMIZE
    if(i<0 || j<0 || i>=n || j>=m || !data)
      error("error in integer i=%ld j=%ld n=%ld m=%ld\n", i,j,n,m);
#endif
    TT d = data[i+j*n];
    return (int) (d + (d<0? -0.5: 0.5));
  }
  
  // returns true iff m1 has the same dimensions of this matrix  
  inline bool matchsize(matrix const &m1) const {
    return n == m1.n && m == m1.m;
  }
  inline tint size() const { return n*m; }

  // subtract matrix other from the current matrix and keep the result
  // on the current matrix
  matrix &operator-=(const matrix &other);
  matrix operator-(const matrix &other) const;
  // negation of matrix
  matrix operator-() const;
  // multiply element-wise this and matrix other
  matrix &operator*=(const matrix &other);
  matrix &operator*=(const TT &other);
  matrix operator*(const matrix &other) const;
  matrix operator*(const TT &other) const;
  matrix operator/(const TT &other) const;
	// matrix product
	matrix multiply(const matrix &other) const;
	matrix inverse() const;

	matrix &operator/=(const matrix &other);
	matrix operator/(const matrix &other) const;
	matrix sqrt() const;
  matrix abs() const;
  matrix loge() const;
	void applyloge();

	matrix &basicop(matrix &tmp, const matrix &other, int type, TT x=0) const;	

	
  bool pattern(const matrix &other, TT p1, TT p2) const;

  // sum element-wise this and matrix other
  matrix operator+(const TT &other) const;
  matrix operator+(const matrix &other) const;
  // add matrix other from the current matrix and keep the result
  // on the current matrix
  matrix &operator+=(const TT &other);
  matrix &operator+=(const matrix &other);
  // returns a submatrix with the columns from j to jj-1
  // if jj is not defined, returns a single column j.
  matrix column(tint j, tint jj=-1) const {
    if(jj<0) jj=j+1;
    if(j==jj) return matrix();
#ifndef OPTIMIZE
    if(j<0 || j>=jj || jj>m) error("column j=%ld to %ld n=%ld m=%ld\n", j, jj, n, m);
#endif
    return submat(0,n,j,jj);
  }

  // update the column of the current matrix with the given column
  void updatecolumn(tint j, const matrix &col);
  // update the line of the current matrix with the given line
  void updateline(tint i, const matrix &lin);
  // returns a transpose matrix
  matrix T() const;
  // returns a submatrix with the lines from j to jj-1
  // if jj is not defined, returns a single line j.
  matrix line(tint j, tint jj=-1) const {
    if(jj<0) jj=j+1;
    if(j==jj) return matrix();
#ifndef OPTIMIZE
    if(j<0 || j>=jj || jj>n) error("line i=%ld to %ld n=%ld m=%ld\n", j, jj, n, m);
#endif
    return submat(j,jj,0,m);
  }

  // concatenate other at the end of the current matrix (end with the respect to the columns),
  // that is, columns of other are included to the right of the current matrix  
  matrix &columncat(const matrix &other);

  // concatenate other at the end of the current matrix (end with the respect to the lines),
  // that is, lines of other are included to the bottom of the current matrix  
  matrix &linecat(const matrix &other);
  // copy operator. It desaloc the current matrix and alloc again with the correct
  // sizes to fit other  
  matrix &operator=(const matrix &other);  
  matrix &realoc(tint newn, tint newm=-1);

  // returns a new matrix containing a sub-matrix of the current matrix
  // from line inin to line fimn (not including fimn), 
  // and from column inim to column fimm (not including fimm)
  matrix submat(tint inin, tint fimn, tint inim, tint fimm) const;

  // desallocation function to destroy the matrix. It needs to delete data
  // in case the flag desaloc is active.
  ~matrix() {
    clear();
  }
  void clear();
  void categorize(matrix &allint, matrix &categ);
  matrix getjoint(const matrix &bb, matrix &bret) const;

  tint nnonzero() const;
  bool haspositive() const;
  bool hasnegative() const;

  static matrix singleton(TT x);
};

template <class TT> matrix<TT> *matrix<TT>::sortcolstemp = 0;
template <class TT> double matrix<TT>::kbs=0.;
template <class TT> double matrix<TT>::megas() { return matrix<double>::kbs/1024.; }

template <class TT>
void matrix<TT>::clear() {
  if(data && desaloc) {
    matrix<double>::kbs -= (double) size()*sizeof(TT)/1024.;
    delete [] data;
  }
  data = 0;
  n = m = 0;
}

// function to aloc the space necessary to keep n lines and m columns
template <class TT>
void matrix<TT>::aloc(int tozero, double zero) {
  if(n<=0 || m<=0) { n=m=0; data=0; }
  else {
	  data = new TT[size()];
	  matrix<double>::kbs += (double) size()*sizeof(TT)/1024.;
	  if(tozero) fill(zero);
	  colnames.resize(m,"");
	  linenames.resize(n,"");
  }
}

// copy the content of the other matrix to the current matrix, but does not
// allocate memory for that. So the memory must be already available.
template <class TT>
void matrix<TT>::copy(const matrix &other) {
  matrix &tmp = *this;
  myassert(matchsize(other));
  for(tint i=0; i<n; i++)
    for(tint j=0; j<m; j++)
      tmp(i,j) = other(i,j);
}

template <class TT>
int matrix<TT>::colcmp(const void *a, const void *b) {
  tint aa = *(tint *) a;
  tint bb = *(tint *) b;
  tint n = sortcolstemp->n;
  for(tint i=0; i<n; i++) {
    TT da = (*sortcolstemp)(i,aa), db = (*sortcolstemp)(i,bb);
    if(fabs(da-db)<eps) continue;
    if(da<db) return -1;
    if(db<da) return 1;
  }
  return 0;
}

template <class TT>
void matrix<TT>::fill(TT zero) {
  for(tint i=0; i<size(); i++) data[i]=zero;
}

template <class TT>
bool matrix<TT>::isvector() const {
  return ((n==1 || m==1) && !empty());
}

template <class TT>
tint matrix<TT>::njumps() const {
	return jumps().size();
}
template <class TT>
matrix<int> matrix<TT>::jumps() const {
  if(!isvector()) throw std::string("jumps only works for vectors");
  tint nj=0;
  const matrix<TT> &temp = *this;
  tint mm = size();
  matrix<int> resp(mm);
  TT d = temp(0);
  for(tint i=0; i<mm; i++) {
    if(fabs(temp(i)-d)>eps) {
      d = temp(i);
	  resp(nj++) = i;
    }
  }
  resp.realoc(nj);
  return resp;
}

template <class TT>
matrix<double> matrix<TT>::todouble() const {
  const matrix<TT> &thi = *this;
  matrix<double> tmp(n,m);
  for(tint i=0; i<n; i++)
    for(tint j=0; j<m; j++)
      tmp(i,j) = (double) thi(i,j);
  return tmp;
}
template <class TT>
matrix<int> matrix<TT>::toint() const {
  const matrix<TT> &thi = *this;
  matrix<int> tmp(n,m);
  for(tint i=0; i<n; i++)
    for(tint j=0; j<m; j++)
      tmp(i,j) = (int) (thi(i,j) < 0 ? (thi(i,j)-0.5): (thi(i,j)+0.5));
  return tmp;
}
template <class TT>
matrix<long> matrix<TT>::tolong() const {
  const matrix<TT> &thi = *this;
  matrix<long> tmp(n,m);
  for(tint i=0; i<n; i++)
    for(tint j=0; j<m; j++)
      tmp(i,j) = (long) (thi(i,j) < 0 ? (thi(i,j)-0.5): (thi(i,j)+0.5));
  return tmp;
}
template <class TT>
matrix<char> matrix<TT>::tochar() const {
  const matrix<TT> &thi = *this;
  matrix<char> tmp(n,m);
  for(tint i=0; i<n; i++)
    for(tint j=0; j<m; j++)
      tmp(i,j) = (char) (thi(i,j) < 0 ? (thi(i,j)-0.5): (thi(i,j)+0.5));
  return tmp;
}

/*
template <class TT>
void matrix<TT>::sort(int (*compar)(const void *, const void *)) {
  qsort(data, size(), sizeof(TT), compar);
}
template <class TT>
void matrix<TT>::sortint() {
  qsort(data, size(), sizeof(TT), intcmp);
}
*/

// returns an integer array of the indices of the columns ordered by
// their elements. Matrices with multiple lines will have the elements
// ordered using all lines, comparing first line, then second, and so on.
// TThe matrix is not altered and only an array of indices corresponding
// to the desired ordering is returned.
template <class TT>
matrix<double> matrix<TT>::sortcols() {
  sortcolstemp = this;
  tint *ind = new tint[m];
  for(tint i=0; i<m; i++) ind[i]=i;
  qsort(ind, m, sizeof(tint), colcmp);
  matrix<double> t = matrix<TT>::tomatrix(ind,m);
  delete [] ind;
  return t;
}

template <class TT>
matrix<TT> matrix<TT>::tomatrix(int *v, tint n) {
  matrix<TT> t(n);
  for(tint i=0; i<n; i++) t(i)=(TT) v[i];
  return t;
}
template <class TT>
matrix<TT> matrix<TT>::tomatrix(long *v, tint n) {
  matrix<TT> t(n);
  for(tint i=0; i<n; i++) t(i)=(TT) v[i];
  return t;
}
template <class TT>
matrix<TT> matrix<TT>::tomatrix(double *v, tint n) {
  matrix<TT> t(n);
  for(tint i=0; i<n; i++) t(i)=(TT) v[i];
  return t;
}

// return true iff this matrix doesnt contain the same data as the other
// and of course have the same size
template <class TT>
bool matrix<TT>::operator!=(const matrix &other) const {
  const matrix &tmp = *this;
  if(!matchsize(other)) return true;
  for(tint i=0; i<n; i++)
    for(tint j=0; j<m; j++)
      if(fabs(tmp(i,j)-other(i,j))>eps) return true;
  return false;
}
template <class TT>
bool matrix<TT>::operator==(const matrix &other) const {
  return !(*this != other);
}

template <class TT>
matrix<TT> matrix<TT>::operator>(TT val) const {
  matrix<TT> resp(n,m,0);
  const matrix &tmp=*this;
  for(tint i=0; i<n; i++) {
    for(tint j=0; j<m; j++)
      if(tmp(i,j)>val) resp(i,j)=1;
  }
  return resp;
}

template <class TT>
bool matrix<TT>::operator<(const matrix &other) const {
  const matrix &tmp = *this;
  if(tmp.N() < other.N()) return true;
  if(tmp.N() > other.N()) return false;
  if(tmp.M() < other.M()) return true;
  if(tmp.M() > other.M()) return false;
  
  for(tint i=0; i<n; i++)
    for(tint j=0; j<m; j++)
      if(tmp(i,j)<(double) other(i,j)-eps) return true;
      else if(tmp(i,j)>(double) other(i,j)+eps) return false;
  return false;
}

template <class TT>
void matrix<TT>::categorize(matrix &allint, matrix &categ) {
  matrix &r = *this;
  for(int j=0; j<r.M(); j++) {
    allint(j) = 1;
    categ(j) = 0;
    // verify if the column is composed only by integers
    for(int i=0; i<r.N(); i++) {
      if(fabs(r(i,j) - r.integer(i,j))>eps) {
	allint(j) = 0;
	break;
      }
    }
    // if column is all integer, verify if there are not so many of them
    if(allint(j)) {
      //      fprintf(stderr,"col %d all int\n", j);

      // the idea is to sort the numbers and then count the number of jumps
      // this gives us the number of distinct numbers
      matrix cc = r.column(j);
      //      cc.T().print();
      matrix c = cc.distinct();
      //      c.T().print();

      // if there are less than 10 distinct numbers, we define the feature column
      // as categorical
      if(c.size()<10) {
	//	fprintf(stderr, "col %d has %d categories:", j, q);

	// now we just ensure that the categories range from 0 to q-1, in case
	// the original numbers were not continuous integer numbers
	for(int ii=0; ii<r.N(); ii++) {
	  int i;
	  for(i=0; i<c.size(); i++)
	    if(r.integer(ii,j)==c.integer(i)) break;
	  myassert(i<c.size());
	  r(ii,j) = i;
	}
	categ(j) = c.size();
	//	fprintf(stderr, " (%.1lf)\n", categ(j));
      }
    }
  }
}

template <class TT>
matrix<TT> matrix<TT>::getjoint(const matrix &bb, matrix &bret) const {
  myassert(matchsize(bb));
  myassert(iscolumn() && bb.iscolumn());
  matrix a(*this), b(bb);
  matrix cata(1), inta(1), catb(1), intb(1);
  a.categorize(inta, cata);
  b.categorize(intb, catb);
  myassert(cata.integer(0)>0);
  myassert(catb.integer(0)>0);
  matrix p(cata.integer(0), catb.integer(0));
  bool doit = bret.empty();
  if(doit) {
    bret.realoc(catb.integer(0));
    for(int i=0; i<catb.integer(0); i++)
      bret(i) = 0;
  }
  for(int i=0; i<p.N(); i++) {
    for(int j=0; j<p.M(); j++) p(i,j)=0;
  }
  for(int i=0; i<a.N(); i++) {
    p(a.integer(i), b.integer(i))++;
    if(doit)
      bret(b.integer(i))++;
  }
  return p;
}


template <class TT>
double matrix<TT>::entropy() const {
  const matrix &tmp = *this;
  double ent = 0;
  tint ene=0;
  for(tint i=0; i<n; i++)
    for(tint j=0; j<m; j++) ene += tmp(i,j);
  
  for(tint i=0; i<n; i++)
    for(tint j=0; j<m; j++) {
      double p = tmp(i,j)/ene;
      if(p>eps)
	ent -= p * log(p);
    }
  return ent/log(2);
}

template <class TT>
double matrix<TT>::condentropy() const {
  const matrix &tmp = *this;
  matrix p(tmp);
  tint ene=0;
  for(tint i=0; i<n; i++)
    for(tint j=0; j<m; j++) ene += tmp(i,j);
  
  for(tint i=0; i<p.N(); i++)
    for(tint j=0; j<p.M(); j++) p(i,j) /= ene;
  double *pcond = new double[n];
  for(tint i=0; i<p.N(); i++) {
    pcond[i] = 0;
    for(tint j=0; j<p.M(); j++)
      pcond[i] += p(i,j);
  }
  double mi = 0;
  for(tint i=0; i<p.N(); i++)
    for(tint j=0; j<p.M(); j++)
      if(p(i,j)>eps && pcond[i]>eps)
	mi += p(i,j) * log(pcond[i]/p(i,j));
  
  delete [] pcond;
  return mi/log(2);
}

template <class TT>
double matrix<TT>::chisquared() const {
  const matrix &tmp = *this;
  tint ene=0;
  double *prow = new double[n];
  for(tint i=0; i<n; i++) {
    prow[i] = 0;
    for(tint j=0; j<m; j++)
      prow[i] += tmp(i,j);
    ene += prow[i];
  }
  double *pcol = new double[n];
  for(tint j=0; j<m; j++) {
    pcol[j] = 0;
    for(tint i=0; i<n; i++)
      pcol[j] += tmp(i,j);
  }
  
  double v = 0;
  for(tint i=0; i<n; i++)
    for(tint j=0; j<m; j++)
      if(tmp(i,j)>eps) {
	double d = (prow[i] * pcol[j]) / ene - tmp(i,j);
	v += (d*d)/(tmp(i,j));
      }
  
  delete [] prow;
  delete [] pcol;
  return v/ene;
}

// writes the binary data of the matrix to the file, that must be
// already open and with writing permission.  
template <class TT>
void matrix<TT>::tofile(FILE *fp) const {
  myassert(fwrite(this, sizeof(matrix), 1, fp)>0);
  if(size() > 0)
	  myassert(fwrite(data, sizeof(TT), size(), fp)>0);
}

// reads from the file fp binary data to populate the current matrix.
// It desaloc the current data e alocs space for the new data.
template <class TT>
void matrix<TT>::fromfile(FILE *fp) {
  if(data && desaloc) {
    matrix<double>::kbs -= (double) size()*sizeof(TT)/1024.;
    delete [] data;
  }
  myassert(fread(this, sizeof(matrix), 1, fp)>0);
  aloc();
  if(size() > 0)
	  myassert(fread(data, sizeof(TT), size(), fp)>0);
}

// writes an array of matrices to a file. fn must contain the filename,
// lis contains nmat elements, each one is a matrix. TThe number of
// matrices being written is also put in the file.
template <class TT>
void matrix<TT>::tofile(const char *fn, const matrix<TT> *lis, int nmat) {
  FILE *fp = fopen(fn,"w");
  myassert(fp);
  myassert(fwrite(&nmat, sizeof(int), 1, fp)>0);
  for(tint i=0; i<nmat; i++)
    lis[i].tofile(fp);
  fclose(fp);
}

// reads an array of matrices from a file. fn is the filename, and
// nmat is used to return the number of elements that is read.
template <class TT>
matrix<TT> *matrix<TT>::fromfile(const char *fn, int *nmat) {
  FILE *fp = fopen(fn,"r");
  myassert(fp);
  myassert(fread(nmat, sizeof(int), 1, fp)>0);
  matrix *lis = new matrix[*nmat];
  for(tint i=0; i<*nmat; i++)
    lis[i].fromfile(fp);
  fclose(fp);
  return lis;
}

template <class TT>
void matrix<TT>::tocsv(std::string *s, FILE *saida) const {
	for(tint i=0; i<m; i++) 
		if(s)
			fprintf(saida,"\"%s\"%s", s[i].c_str(), (i<m-1?",":""));
		else
			fprintf(saida,"v%ld%s", i+1, (i<m-1?",":""));
	fprintf(saida,"\n");

	for(tint i=0; i<n; i++) {  
		for(tint j=0; j<m; j++)
			fprintf(saida,"%lf%s", (double) (*this)(i,j), (j<m-1?",":""));
		fprintf(saida,"\n");
	}
}
template <class TT>
matrix<TT> matrix<TT>::fromcsv(FILE *ent, tint n, tint m, bool lineheader, bool colheader) {
	const int s = 10000000;
	char *tt = 0;
	char *lin = new char[s+2];
	if(m==0 || n==0) {
	  myassert(fgets(lin, s, ent)!=NULL);
//	  std::cerr << "lin = " << lin << std::endl;
	  m=0;
	  lineheader=false;
	  colheader=false;
	  char *tok = std::strtok(lin,",");
	  while(tok) {
	    tt = 0;
		double d = strtod(tok, &tt);
//		std::cerr << "tok=" << tok << " " << (unsigned long) tok << "tt=" << (unsigned long) tt << std::endl;
		if(tt==tok) lineheader=true;
	    if(std::strchr(tok,'\"') != NULL || std::strchr(tok,'\'') != NULL) lineheader=true;
	    m++;
	    tok = std::strtok(NULL,",");
	  }
	  n=1;
	  while(fgets(lin, s, ent)!=NULL) {
	    tok = std::strtok(lin,",");
	    myassert(tok);
	    tt = 0;
        double d = strtod(tok, &tt); 
//		std::cerr << "tok=" << tok << " " << (unsigned long) tok << "tt=" << (unsigned long) tt << std::endl;
		if(tt==tok) colheader=true;
	    if(std::strchr(tok,'\"') != NULL || std::strchr(tok,'\'') != NULL) colheader=true;
	    n++;
	  }
	  fseek(ent,0L,SEEK_SET);
	  if(colheader) m--;
	  if(lineheader) n--;
	}
	matrix<> mat(n,m);
	std::cerr << "n=" << n << " m=" << m << " lineheader=" << lineheader << " colheader=" << colheader << std::endl;
	//	myassert(fgets(lin, s, ent)!=NULL);
	if(lineheader) {
		myassert(fgets(lin, s, ent)!=NULL);
//		std::cerr << "line header = " << lin << std::endl;
	}
	for(tint i=0; i<n; i++) {
		myassert(fgets(lin, s, ent)!=NULL);
		char *tok = std::strtok(lin,",");
		if(colheader) {
			tok = std::strtok(NULL,",");
//			std::cerr << "col header = " << tok << std::endl;
		}
		for(tint j=0; j<m; j++) {
			myassert(tok);
			while(*tok == ' ' || *tok == '\"' || *tok == '\'') tok++;
			char *t;
			for(t=tok; *t && *t != '\'' && *t != '\"' && *t != ' '; t++);
			*t = 0;
			mat(i,j) = atof(tok);
			tok = std::strtok(NULL,",");
		}
	}
	delete [] lin;
	return mat;
}

#ifndef CPURE
template <class TT>
void matrix<TT>::print(std::string s, tint x, tint y, bool typeint, FILE *saida) const {
	if(x<0 || x>n) x=n;
	if(y<0 || y>m) y=m;
	Rprintf("matrix::print(%s,%ld,%ld)\n",s.c_str(),n,m);
	for(tint i=0; i<x; i++) {
		for(tint j=0; j<y; j++)
			if(typeint)
				Rprintf("%ld ", (long int) (*this).integer(i,j));
			else
				Rprintf("%lf ", (double) (*this)(i,j));
		Rprintf("\n");
  }
}
#else
template <class TT>
void matrix<TT>::print(std::string s, tint x, tint y, bool typeint, FILE *saida) const {
  if(x<0 || x>n) x=n;
  if(y<0 || y>m) y=m;
  fprintf(saida,"matrix::print(%s,%ld,%ld)\n",s.c_str(),n,m);
  for(tint i=0; i<x; i++) {
    for(tint j=0; j<y; j++)
      if(typeint)
	fprintf(saida,"%ld ", (long int) (*this).integer(i,j));
      else
	fprintf(saida,"%lf ", (double) (*this)(i,j));
    fprintf(saida,"\n");
  }
}
#endif

// returns the minimum element of column j, or line i, or all
template <class TT>
TT matrix<TT>::min(tint i, tint j) const {
  return operation(2, i, j);
}

// returns the maximum element of column j, or line i, or all
template <class TT>
TT matrix<TT>::max(tint i, tint j) const {
  return operation(1, i, j);
}

template <class TT>
TT matrix<TT>::findmax(tint i, tint j) const {
	myassert(i<-eps || j<-eps);
  return operation(4, i, j);
}


// returns the sum of the matrix elements of column j, or line i, or all
template <class TT>
TT matrix<TT>::sum(tint i, tint j) const {
  return operation(0, i, j);
}

// dim == 1 : column-wise, 
// dim == 2 : line-wise
enum {
	COLUMNWISE = 1,
	LINEWISE = 2
};
// type == 1: sum 
// type == 2: which.max
enum {
	SUM = 1,
	WHICHMAX = 2,
	VALMAX = 3,
};
template <class TT>
matrix<TT> matrix<TT>::apply(int dim, int type) const {
	matrix<TT> r;
	switch(dim) {
		case COLUMNWISE:
			r.realoc(M());
			for(tint i=0; i<M(); i++)
				switch(type) {
					case SUM:
						r(i) = sum(-1,i);
						break;
					case WHICHMAX:
						r(i) = findmax(-1,i);
						break;
					case VALMAX:
						r(i) = r(findmax(-1,i),i);
						break;
				}
			break;
		case LINEWISE:
			r.realoc(1,N());
			for(tint i=0; i<N(); i++)
				switch(type) {
					case SUM:
						r(i) = sum(i,-1);
						break;
					case WHICHMAX:
						r(i) = findmax(i,-1);
						break;
				}
			break;
		default: myassert(1==0);
	}
	return r;
}

template <class TT>
TT matrix<TT>::operation(int type, tint iz, tint jz) const {
  myassert(iz<n && jz<m);
  const matrix &tmp = *this;
  if(tmp.empty()) return 0.;

  tint inii=0, fimi=n-1;
  tint inij=0, fimj=m-1;
  if(iz>=0) inii=fimi=iz;
  if(jz>=0) inij=fimj=jz;

  TT s=0;
  tint pos=-1;
  TT d;
  switch(type) {
	  case 1:
	  case 2:
		  s = tmp(inii,inij);
		  break;
	  case 4:
		  s = -1;
		  break;
  }
  for(tint i=inii; i<=fimi; i++)
	  for(tint j=inij; j<=fimj; j++)
		  switch(type) {
			  case 0:
				  s += tmp(i,j); 
				  break;
			  case 1:
				  if(tmp(i,j)>s) s = tmp(i,j);
				  break;
			  case 2:
				  if(tmp(i,j)<s) s = tmp(i,j);
				  break;
			  case 3:
				  d = tmp(i,j);
				  s += d*d;
				  break;
			  case 4:
				  if(s<-eps || tmp(i,j)>s) {
					  s = tmp(i,j);
					  pos = (iz>=0? j: i);
				  }
				  break;
			case 5:
			  if(fabs(tmp(i,j)) > s) s = tmp(i,j);
			  break;
		  }
  if(type != 4)
	  return s;
  return pos;
}

// returns the absolute maximum value of the matrix elements of column j, or line i, or all
template <class TT>
double matrix<TT>::absmax(tint i, tint j) const {
  if(size()>0) return operation(5,i,j);
  else return 0.;
}


// returns the mean of the matrix elements of column j, or line i, or all
template <class TT>
double matrix<TT>::mean(tint i, tint j) const {
  if(size()>0) return sum(i,j)/size();
  else return 0.;
}

// returns the mean of squared elements of the matrix
template <class TT>
double matrix<TT>::meansquare(tint i, tint j) const {
  if(size()>0) return operation(3,i,j)/size();
  else return 0.;
}

// subtract matrix other from the current matrix and keep the result
// on the current matrix
template <class TT>
matrix<TT> &matrix<TT>::operator-=(const matrix &other) {
	matrix &tmp = *this;
	return basicop(tmp,other,2);
}

template <class TT>
matrix<TT> matrix<TT>::operator-(const matrix &other) const {
	matrix<TT> tmp(*this);
	return basicop(tmp,other,2);
}

// negation of matrix
template <class TT>
matrix<TT> matrix<TT>::operator-() const {
	matrix<TT> tmp(*this);
	return basicop(tmp,tmp,5);
}

// multiply element-wise this and matrix other
template <class TT>
matrix<TT> matrix<TT>::operator*(const matrix &other) const {
	matrix<TT> tmp(*this);
	return basicop(tmp,other,3);
}
template <class TT>
matrix<TT> &matrix<TT>::operator*=(const matrix &other) {
	matrix<TT> &tmp(*this);
	return basicop(tmp,other,3);
}

template <class TT>
matrix<TT> &matrix<TT>::basicop(matrix &tmp, const matrix &other, int type, TT x) const {
	myassert(matchsize(other) || (m==other.M() && other.N()==1) || (n==other.N() && other.M()==1) || (other.N()==1 && other.M()==1));
	tint k = 0, kk = 0;
	if(type == 7)
		myassert(fabs(x)>eps);
	for(tint i=0; i<n; i++)
		for(tint j=0; j<m; j++) {
			if(type < 6) {
				if(other.N()==1) k=0;
				else k = i;
				if(other.M()==1) kk=0;
				else kk = j;
			}
			switch(type) {
				case 1:		  
					tmp(i,j) += other(k,kk);
					break;
				case 2:
					tmp(i,j) -= other(k,kk);
					break;
				case 3:
					tmp(i,j) *= other(k,kk);
					break;
				case 4:
					myassert(fabs(other(k,kk))>eps);
					tmp(i,j) /= other(k,kk);
					break;
				case 5:
					tmp(i,j) = -other(k,kk);
					break;
				case 6:
					tmp(i,j) *= x;
					break;
				case 7:
					tmp(i,j) /= x;
					break;
				case 8:
					tmp(i,j) += x;
					break;
			}
		}
	return tmp;
}

// take square root element-wise
template <class TT>
matrix<TT> matrix<TT>::sqrt() const {
	const matrix<TT> &th(*this);
	matrix<TT> tmp(N(),M());
	for(tint i=0; i<n; i++)
		for(tint j=0; j<m; j++)
			tmp(i,j) = std::sqrt(th(i,j));
	return tmp;
}

// absolute element-wise value
template <class TT>
matrix<TT> matrix<TT>::abs() const {
  const matrix<TT> &th(*this);
  matrix<TT> tmp(N(),M());
  for(tint i=0; i<n; i++)
	for(tint j=0; j<m; j++)
	  tmp(i,j) = fabs(th(i,j));
  return tmp;
}
// log element-wise value
template <class TT>
matrix<TT> matrix<TT>::loge() const {
  const matrix<TT> &th(*this);
  matrix<TT> tmp(N(),M());
  for(tint i=0; i<n; i++)
	for(tint j=0; j<m; j++) {
	  myassert(th(i,j)>0.);
	  tmp(i,j) = log(th(i,j));
	}
  return tmp;
}
template <class TT>
void matrix<TT>::applyloge() {
  matrix<TT> &th(*this);
  for(tint i=0; i<n; i++)
	for(tint j=0; j<m; j++) {
	  myassert(th(i,j)>0.);
	  th(i,j) = log(th(i,j));
	}
}

// multiply element-wise this and the constant
template <class TT>
matrix<TT> matrix<TT>::operator*(const TT &other) const {
	matrix<TT> tmp(*this);
	return basicop(tmp,tmp, 6, other);
}
template <class TT>
matrix<TT> &matrix<TT>::operator*=(const TT &other) {
	matrix<TT> &tmp(*this);
	return basicop(tmp,tmp, 6, other);
}
template <class TT>
matrix<TT> matrix<TT>::operator/(const TT &other) const {
	matrix<TT> tmp(*this);
	return basicop(tmp,tmp, 7, other);
}

// multiply this and matrix other (product of the matrices)
template <class TT>
matrix<TT> matrix<TT>::multiply(const matrix &other) const {
	const matrix<TT> &th(*this);
	matrix<TT> tmp(N(),other.M());
	myassert(M() == other.N());

// 	cblas_dgemm(ML_BLAS_STORAGE,CblasNoTrans,CBlasNoTrans,
//             A.rows(),B.cols(),A.cols(),
//             1.0,A.raw_data(),A.data_stride(),
//             B.raw_data(),B.data_stride(),
//             0.0,C.raw_data(),C.data_stride());

//	cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, tmp.N(), tmp.M(), M(), 1.0, th.data, M(), other.data, other.M(), 0.0, tmp.data, tmp.M());

#ifdef USE_CBLAS
	cblas_dgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, tmp.N(), tmp.M(), M(), 1.0, th.data, N(), other.data, other.N(), 0.0, tmp.data, tmp.N());
#else
	for(tint i=0; i<tmp.N(); i++)
		for(tint j=0; j<tmp.M(); j++) {
		    tmp(i,j) = 0.;
			for(tint k=0; k<M(); k++)
				tmp(i,j) += th(i,k) * other(k,j);
		}
#endif
	return tmp;
}

// inverse
template <class TT>
matrix<TT> matrix<TT>::inverse() const {
	const matrix<TT> &th(*this);

	tint m = M();
	if(m!=N()) throw "Only Square Matrix can have inverse!";
	matrix SI(m,2*m);

	for(tint i=0; i<m;i++) {
		for(tint j=0; j<2*m; j++) {
			if(j<m) {
				SI.data[i*2*m+j] = th(i,j);
			} else {
				SI.data[i*2*m+j] = 0;
			}
			if(j==i+m) SI.data[i*2*m+j] =1.0;
		}
	}
	for(tint ii=0; ii<m-1; ii++) {
		for(tint kk=ii+1; kk < m; kk++) {
			for(tint jj=2*m-1; jj>=ii; jj--) {
				SI.data[kk*2*m+jj] -= SI.data[ii*2*m+jj] * 
					(SI.data[kk*2*m+ii]/SI.data[ii*2*m+ii]);
			}
		}
	}      
	for(tint iii=m-1; iii>0; iii--) {
		for(tint kkk=iii-1; kkk >=0; kkk--){
			for(tint jjj=2*m-1; jjj>=iii; jjj--){
				SI.data[kkk*2*m+jjj] -= SI.data[iii*2*m+jjj] * 
					(SI.data[kkk*2*m+iii]/SI.data[iii*2*m+iii]);
			}
		}
	}      
	
	for(tint iiii=0; iiii<m; iiii++) {
		for(tint jjjj=2*m-1; jjjj >=0; jjjj--) {
			SI.data[iiii*2*m+jjjj] =SI.data[iiii*2*m+jjjj]/SI.data[iiii*2*m+iiii];                     
		}
	}      

	matrix I(m,m);
	for(int s = 0; s<m; s++)
		for(int r=0; r<m; r++) {
			I.data[s*m+r] = SI.data[s*2*m+r+m];
		}
	return I;
}

// divide element-wise this and matrix other
template <class TT>
matrix<TT> matrix<TT>::operator/(const matrix &other) const {
	matrix tmp(*this);
	return basicop(tmp,other,4);
}
template <class TT>
matrix<TT> &matrix<TT>::operator/=(const matrix &other) {
	matrix &tmp(*this);
	return basicop(tmp,other,4);
}

template <class TT>
bool matrix<TT>::pattern(const matrix &other, TT p1, TT p2) const {
  tint e = size();
  for(tint i=0; i<e; i++)
    if(data[i] == p1 && other.data[i] == p2) return true;
  return false;
}

// sum element-wise this and matrix other
template <class TT>
matrix<TT> matrix<TT>::operator+(const TT &other) const {
	matrix<TT> tmp(*this);
	return basicop(tmp,tmp, 8, other);
}
template <class TT>
matrix<TT> &matrix<TT>::operator+=(const TT &other) {
	matrix<TT> &tmp(*this);
	return basicop(tmp,tmp, 8, other);
}

template <class TT>
matrix<TT> matrix<TT>::operator+(const matrix &other) const {
	matrix<TT> tmp(*this);
	return basicop(tmp,other,1);
}

// add matrix other from the current matrix and keep the result
// on the current matrix
template <class TT>
matrix<TT> &matrix<TT>::operator+=(const matrix &other) {
	matrix<TT> &tmp = *this;
	return basicop(tmp,other,1);
}

// update the column of the current matrix with the given column
template <class TT>
void matrix<TT>::updatecolumn(tint j, const matrix &col) {
  myassert(j>=0 && j<m);
  myassert(col.N() == n);
  matrix<TT> &tmp = *this;
  for(tint i=0; i<n; i++) tmp(i,j) = col(i);
}
// update the line of the current matrix with the given line
template <class TT>
void matrix<TT>::updateline(tint i, const matrix &lin) {
  myassert(i>=0 && i<n);
  myassert(lin.M() == m);
  matrix<TT> &tmp = *this;
  for(tint j=0; j<m; j++) tmp(i,j) = lin(0,j);
}

// returns a transpose matrix
template <class TT>
matrix<TT> matrix<TT>::T() const {
  matrix<TT> t(m,n);
  for(tint i=0; i<n; i++) 
    for(tint j=0; j<m; j++)
      t(j,i) = (*this)(i,j);
  return t;
}

  // concatenate other at the end of the current matrix (end with the respect to the columns),
  // that is, columns of other are included to the right of the current matrix  
template <class TT>
matrix<TT> &matrix<TT>::columncat(const matrix &other) {
  if(empty()) *this = other;
  else {
    if(!other.empty()) {
      matrix<TT> &tmp = *this;
      tint oldm = m;
      if(other.n != n) error("columncat number of lines must agree (%ld %ld) <-> (%ld %ld)\n",n,m,other.n,other.m);
      realoc(n,m+other.m);
      for(tint i=0; i<other.n; i++)
	for(tint j=0; j<other.m; j++)
	  tmp(i,j+oldm) = other(i,j);
    }
  }
  return *this;
}

// concatenate other at the end of the current matrix (end with the respect to the lines),
// that is, lines of other are included to the bottom of the current matrix  
template <class TT>
matrix<TT> &matrix<TT>::linecat(const matrix &other) {
  if(empty()) *this = other;
  else {
    if(!other.empty()) {
      matrix<TT> &tmp = *this;
      tint oldn = n;
      if(other.m != m) error("linecat number of columns must agree (%ld %ld) <-> (%ld %ld)\n",n,m,other.n,other.m);
      realoc(n+other.n,m);
      for(tint i=0; i<other.n; i++)
	for(tint j=0; j<other.m; j++)
	  tmp(i+oldn,j) = other(i,j);
    }
  }
  return *this;
}

// copy operator. It desaloc the current matrix and alloc again with the correct
// sizes to fit other  
template <class TT>
matrix<TT> &matrix<TT>::operator=(const matrix<TT> &other) {
  if(n != other.n || m != other.m) {
    if(data && desaloc) {
      matrix<double>::kbs -= (double) size()*sizeof(TT)/1024.;
      delete [] data;
    }
    n = other.n;
    m = other.m;
    aloc();
  }
  copy(other);
  return *this;
}

 
template <class TT>
matrix<TT> &matrix<TT>::realoc(tint newn, tint newm) {
  if(newm<=0) {
    if(newn>0 && newm<=0) newm=1;
	else newm=m;
  }
  if(newn==n && newm==m) return *this;
  if(newn<=0) {
	  clear();
	  return *this;
  }
  myassert(newn>0 && newm>0);
  TT *newdata = new TT[newn*newm];
  matrix<double>::kbs += (double) newn * newm *sizeof(TT)/1024.;
  for(tint i=0; i<mini(newn,n); i++) 
    for(tint j=0; j<mini(newm,m); j++)
      newdata[i + newn * j] = (*this)(i,j);
  myassert(desaloc);
  if(data && size()>0) {
    matrix<double>::kbs -= (double) size()*sizeof(TT)/1024.;
    delete [] data;
  }
  data = newdata;
  n = newn;
  m = newm;
  return *this;
}

template <class TT>
double matrix<TT>::euclidean(const matrix<TT> &other) const {
  const matrix<TT> &tmp = *this;
  myassert(matchsize(other));
  double d=0;
  for(int i=0; i<n; i++)
    for(int j=0; j<m; j++) {
      double e = tmp(i,j) - other(i,j);
      d += e*e;
    }
  return std::sqrt(d);
}

template <class TT>
double matrix<TT>::correlation(const matrix<TT> &other) const {
  const matrix<TT> &tmp = *this;
  myassert(matchsize(other));
  double sum=0, sumo=0, sumq=0, sumoq=0, summ=0;
  for(int i=0; i<n; i++)
    for(int j=0; j<m; j++) {
      double d = tmp(i,j), e = other(i,j);
      sum += d;
      sumo += e;
      sumq += d*d;
      sumoq += e*e;
      summ += d*e;
    }
  return (n*m*summ - sum*sumo)/(std::sqrt(n*m*sumq - sum*sum)*std::sqrt(n*m*sumoq - sumo*sumo));
}

template <class TT>
void matrix<TT>::maxpos(int &a, int &b) const {
  myassert(n>0 && m>0);
  const matrix &tmp = *this;
  a=0;
  b=0;
  for(int i=0; i<n; i++)
    for(int j=0; j<m; j++)
      if(tmp(i,j)>tmp(a,b)) { a=i; b=j; }
}

template <class TT>
double matrix<TT>::l1distance(const matrix<TT> &other) const {
  const matrix<TT> &tmp = *this;
  myassert(matchsize(other));
  double d=0;
  for(int i=0; i<n; i++)
    for(int j=0; j<m; j++) {
      double e = tmp(i,j) - other(i,j);
      d += fabs(e);
    }
  return d;
}

template <class TT>
matrix<int> matrix<TT>::tocat_bycol() const {
	const matrix<TT> &tmp = *this;
	matrix<int> resp(N(),M(),-1);
	for(tint j=0; j<M(); j++) {
		matrix<TT> dis = column(j).distinct();
		for(tint i=0; i<N(); i++) {
			tint k;
			for(k=0; k<dis.size(); k++)
				if(tmp(i,j)==dis(k)) {
					resp(i,j)=k;
					break;
				}
			myassert(k<dis.size());
		}
	}
	return resp;
}

template <class TT>
matrix<TT> matrix<TT>::distinct() const {
	matrix<TT> tmp = *this;
	matrix<TT> tmp2(n*m,1);
	tmp.sort();
	int k=1;
	tmp2(0) = tmp(0,0);
	for(int i=0; i<n*m-1; i++) {
		if(fabs(tmp.data[i]-tmp.data[i+1])>eps) {
			tmp2(k)=tmp.data[i+1];
			k++;
		}
	}
	return tmp2.submat(0,k,0,1);
}

// returns a new matrix containing a sub-matrix of the current matrix
// from line inin to line fimn (not including fimn), 
// and from column inim to column fimm (not including fimm)
template <class TT>
matrix<TT> matrix<TT>::submat(tint inin, tint fimn, tint inim, tint fimm) const {
  matrix<TT> t(fimn-inin,fimm-inim);
  const matrix<TT> &tmp = *this;
  for(tint i=0; i<t.n; i++)
    for(tint j=0; j<t.m; j++)
      t(i,j) = tmp(i+inin,j+inim);
  return t;
}

template <class TT>
tint matrix<TT>::nnonzero() const {
  const matrix<TT> &tmp = *this;
  tint nz=0;
  for(tint i=0; i<n; i++) for(tint j=0; j<m; j++) if(fabs(tmp(i,j))>eps) nz++;
  return nz;
}
template <class TT>
bool matrix<TT>::haspositive() const {
  const matrix<TT> &tmp = *this;
  for(tint i=0; i<n; i++) for(tint j=0; j<m; j++) if(tmp(i,j)>eps) return true;
  return false;
}
template <class TT>
bool matrix<TT>::hasnegative() const {
  const matrix<TT> &tmp = *this;
  for(tint i=0; i<n; i++) for(tint j=0; j<m; j++) if(tmp(i,j)<-eps) return true;
  return false;
}

template <class TT>
matrix<TT> matrix<TT>::singleton(TT x) {
  matrix<TT> tmp(1);
  tmp(0) = x;
  return tmp;
}

template <class TT>
double matrix<TT>::mutualinfo(const matrix<TT> &other) const {
	matrix<TT> tc;
    matrix<TT> p = getjoint(other, tc);
    return tc.entropy() - p.condentropy();
}

bool issimilar(const matrix<double> &c1, const matrix<double> &c2, double mean, double var, int type, double thres=2.);

#ifndef CPURE
extern "C" {
	SEXP similar(SEXP c1, SEXP c2, SEXP mean, SEXP var, SEXP type, SEXP thres);
}
#endif

#endif
