/*
## nmf.hh (package compactnmf) Copyright (C) 2012 C. P. de Campos (cassiopc@acm.org)

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

#ifndef NMF_HH
#define NMF_HH

#include "matrix.hh"

enum {
	EUCLIDEAN = 1,
	DIVERGENCE,
	MODIFIED,
	LOCAL,
	NNSC,
	NNSC2,
	KMEANS
};

void nmfCpp(matrix<> &v, int r, int verbose, int method, double lambda, double mu,
			double gamma, int reorder, matrix<> &w0, matrix<> &h0, matrix<> &w, matrix<> &h);
void nmfconsensusCpp(matrix<> &a, int kstart, int kend, int nloop, 
		     matrix<>*&opt_w, matrix<>*&opt_h, matrix<>*&consensus, matrix<>*&opt_conn, matrix<>**&ww, matrix<>**&hh, 
		     int verbose, int intmethod, std::string conn_method,
		     double lambda, double mu, double gamma, char *inputfile, bool transp);
void nmfconsensusCppList(matrix<> &a, int kstart, int kend, int iniloop, int fimloop, std::vector<int> looplist, 
			 matrix<>*&opt_w, matrix<>*&opt_h, matrix<>*&consensus, matrix<>*&opt_conn, matrix<>**&ww, matrix<>**&hh, 
			 int verbose, int intmethod, std::string conn_method,
			 double lambda, double mu, double gamma, char *inputfile, bool nopt, bool transp);
void nmfconsensusCppReadFiles(matrix<> &a, int kstart, int kend, int iniloop, int fimloop, std::vector<int> looplist, 
			      matrix<>**&w, matrix<>**&h, int verbose, char *inputfile, bool transp);

#endif
