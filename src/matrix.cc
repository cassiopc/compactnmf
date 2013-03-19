/*
## matrix.cc (package compactnmf) Copyright (C) 2012 C. P. de Campos (cassiopc@acm.org)

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
#include "matrix.hh"

#ifndef CPURE
extern "C" {
	SEXP similar(SEXP c1, SEXP c2, SEXP mean, SEXP var, SEXP type, SEXP thres) {
		return matrix<>::mat2sexp(matrix<>::singleton(issimilar(matrix<>(c1),matrix<>(c2),getd(mean),getd(var),geti(type),getd(thres))));
	}
}
#endif

bool issimilar(const matrix<double> &c1, const matrix<double> &c2, double mean, double var, int type, double thres) {
  bool d=false;
  matrix<double> tc, p;
  switch(type) {
  case 1: 
    d = c1.euclidean(c2) < mean - thres*sqrt(var);
    break;
  case 2:
    d = fabs(c1.correlation(c2)) > mean + thres*sqrt(var);
    break;
  case 3:
    d = c1.l1distance(c2) < mean - thres*sqrt(var);
    break;
  case 4:
    p = c1.getjoint(c2, tc);
    d = p.chisquared() > mean + thres*sqrt(var);
    break;
  default:
    p = c1.getjoint(c2, tc);
    d = tc.entropy() - p.condentropy() > mean + thres*sqrt(var);
  }
  return d;
}

