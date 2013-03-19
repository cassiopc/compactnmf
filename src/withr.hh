/*
## withr.hh (package compactnmf) 
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

#ifndef WITHR_HH
#define WITHR_HH

#ifdef CPURE
#define Rprintf(...) (fprintf(stderr,__VA_ARGS__), fflush(stderr))
#define error(...) ((void) (fprintf(stderr,__VA_ARGS__), backtrace((char *) "",(char *) __FILE__,__LINE__), throw 1, 0))
#define myassert(e) ((void) ((e)? 0: ((void) backtrace((char *) #e,(char *) __FILE__,__LINE__), throw 0, 0)))

// this function dumps the stack in case of a premature abortion of the execution
// char *ass is a string contaning the assertion that has been violated
// char *fname is the name of the source file where the assertion was violated
// int line is the number of the line of the assertion
void backtrace(char *ass, char *fname, int line);
#else
#include <R.h> 
#include <Rinternals.h> 
#include <Rmath.h>
#include <cassert>
#define myassert(e) ((void) ((e)? 0: ((void) error((char *) #e), 0)))
#include <string>
#include <vector>

int geti(SEXP x, int p=0);
double &getd(SEXP x, int p=0);
double& getp(SEXP t, SEXP d, int x, int y);
double getpc(SEXP t, SEXP d, int x, int y);
std::string getstr(SEXP t);
char *getstr(SEXP t, SEXP d, int x, int y);

SEXP str2num(SEXP t);

#endif
#endif
