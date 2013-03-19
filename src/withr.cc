/*
## withr.cc (package compactnmf) 
## Copyright (C) 2012 C. P. de Campos (cassiopc@acm.org)

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

#include "withr.hh"
#include<cstdio>
#include<cstdlib>
#include<cstring>

#ifdef CPURE
#include <cxxabi.h>
#include <execinfo.h>

// this function dumps the stack in case of a premature abortion of the execution
// char *ass is a string contaning the assertion that has been violated
// char *fname is the name of the source file where the assertion was violated
// int line is the number of the line of the assertion
void backtrace(char *ass, char *fname, int line) {
  void *addresses[20];
  char **strings;
  int size = backtrace(addresses, 20);
  strings = backtrace_symbols(addresses, size);
  fprintf(stderr,"Assertion %s failed at file %s line %d\nStack frames: %d\n", ass, fname, line, size);
  for(int i = 0; i < size; i++) {
    fprintf(stderr,"[%d: %lX] ", i, (unsigned long)addresses[i]);
    fprintf(stderr,"%s -- ", strings[i]); 
    char *function = new char[1000];
    char *begin = 0, *end = 0;
    // find the parentheses and address offset surrounding the mangled name
    for (char *j = strings[i]; *j; ++j) {
      if (*j == '(') {
	begin = j;
      }
      else if (*j == '+') {
	end = j;
      }
    }
    if (begin && end) {
      begin++;
      *end = 0;
      // found our mangled name, now in [begin, end)
      
      int status;
      size_t sz=1000;
      char *ret = abi::__cxa_demangle(begin, function, &sz, &status);
      if (ret) {
	// return value may be a realloc() of the input
	function = ret;
      }
      else {
	// demangling failed, just pretend it's a C function with no args
	std::strncpy(function, begin, sz);
	std::strncat(function, "()", sz);
	function[sz-1] = 0;
      }
      fprintf(stderr, "%s (%u,%d)", function,(unsigned)sz,status);
    }
    fprintf(stderr, "\n");
    //    ad.lookup(addresses[i]);
    delete [] function;
  }
  free(strings);
}
#else
int geti(SEXP x, int p) {
  return INTEGER(x)[p];
}
  
double &getd(SEXP x, int p) {
  return REAL(x)[p];
}
  
double& getp(SEXP t, SEXP d, int x, int y) {
  return REAL(t)[(x) + INTEGER(d)[0] * (y)];
}
  
double getpc(SEXP t, SEXP d, int x, int y) {
  const char *s = CHAR(STRING_PTR(t)[(x) + INTEGER(d)[0] * (y)]);
  return atof((char *)s);
}

std::string getstr(SEXP t) {
  const char *s = CHAR(STRING_PTR(t)[0]);
  return std::string(s);
}

char *getstr(SEXP t, SEXP d, int x, int y) {
	return (char *) (STRING_PTR(t)[(x) + INTEGER(d)[0] * (y)]);
}

/*
SEXP str2num(SEXP t) {
  SEXP dim = getAttrib(t, R_DimSymbol);
  int n = INTEGER(dim)[0]; 
  SEXP v;
  PROTECT(v = allocMatrix(REALSXP, n, 4));
  SEXP dimv = getAttrib(v, R_DimSymbol);
  for(int i=0; i<n; i++) {
    double d = getpc(t, dim, i, 1);
    if(fabs(d)<eps) d=23;
    getp(v, dimv, i, 0) = (int) d;
    getp(v, dimv, i, 1) = getpc(t, dim, i, 2);
    getp(v, dimv, i, 2) = getpc(t, dim, i, 3);
    getp(v, dimv, i, 3) = getpc(t, dim, i, 4);
  }
  UNPROTECT(1);
  return v;
}
*/

double sumsquare(SEXP a, SEXP b) {
  SEXP dim = getAttrib(b, R_DimSymbol);
  int n = INTEGER(dim)[0];
  double s = 0;
  for(int i=0; i<n; ++i) {
    double d = getd(a, i) - getd(b, i);
    s += d*d;
  }
  return s;
}
  
#endif
