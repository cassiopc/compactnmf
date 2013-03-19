/*
## nmf-standalone.cc (package compactnmf) 
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

#include <cstdio>
#include <cassert>
#include "nmf.hh"

matrix<> gentest(int n, int m, int r) {
	matrix<> h0(n,r);
	matrix<> w0(r,m);
	for(int i=0; i<w0.N(); i++)
		for(int j=0; j<w0.M(); j++)
			w0(i,j) = randnum(1e6+1)/1.e6; 
	for(int i=0; i<h0.N(); i++)
		for(int j=0; j<h0.M(); j++)
			h0(i,j) = randnum(1e6+1)/1.e6; 
	return h0.multiply(w0);
}

void gentests() {
	int nsizes = 6;
	int sizes[] = { 10, 50, 100, 500, 1000, 5000 };
	int nranks = 3;
	int ranksizes[] = { 2, 4, 8 };
	char *name = new char[100];
	for(int i=nsizes-1; i<nsizes; i++) {
		for(int r=0; r<nranks; r++) {
			for(int j=100; j<201; j+=100) {
				sprintf(name,"test-%d-%d-%d.csv",sizes[i],j,ranksizes[r]);
				FILE *fp = fopen(name,"w");
				myassert(fp);
				matrix<> mat = gentest(sizes[i],j,ranksizes[r]);
				mat.tocsv(0, fp);
				fclose(fp);
			}
		}
	}
	delete [] name;
}

//#define TIPO1
#ifdef TIPO1
int main() {
	gentests();
	return 0;
}
#else
matrix <> origmat, comparetoorigmat, filtermat;
double origsum, klorig;
int nmfold=0;
int main(int argc, char **argv) {
	char sspace[300], dirname[200];
	if(argc < 3) {
		fprintf(stderr, "NMF solver. Copyright (c) 2008-2012 C. P. de Campos (cassiopc@acm.org). All rights reserved.\n\
Usage: %s <inputfile> <rankini> [<rankfim> [<#loops> [<runnum> [<directoryname> [<originalmat> <filtering>]]]]]\n\
\n\
inputfile: complete path of the input file\n\
rankini: initial rank to be processed\n\
         if < 0, then AMPL code is generated from the matrix.\n\
rankfim: final rank to be processed (by default equals to rankini)\n\
runnum: if greater than zero, then it corresponds to the number of the NMF run to be executed (ignored if >0 and #loops > 1).\n\
        if <= -3, then the NMF cross-validation is executed over the already done NMFs.\n\
        if == -2, then the matrix is transposed before analyzing the results.\n\
        if == -2 or -1, then the NMF runs (previously done) are analyzed and put together in the consensus matrix.\n\
#loops: number of times to run the NMF.\n\
directoryname: if the files are to be read and generated in a different directory than the current, specify where\n\
originalmat/filtering: csv files with the original matrix and with a list of the columns that were merged to obtain\n\
        the inputfile as it is. This is used only if one wants to compare the factorization with the original matrix\n\
        when runmun <= -3. The number of lines in filtering must be equal to the number of data/numerical lines in inputfile\n\
        (if you check the files, filtering may have one line less, because inputfile may have a header line\n\
", argv[0]);
		return 1;
	}

	FILE *inputfile = fopen(argv[1],"r");
	myassert(inputfile);
	initrand();

	Rprintf("reading file %s\n",argv[1]);
	matrix<> mat = matrix<>::fromcsv(inputfile);
	//	mat.print("mat");

	int rini = atoi(argv[2]);

	if(rini < 0) {
	  printf("param k := %d;\n", -rini);
	  printf("param N := %ld;\n", mat.N());
	  printf("param M := %ld;\n", mat.M());
	  printf("set matlins :=");
	  for(int i=0; i<mat.N(); ++i) printf(" %d", i+1);

	  printf(";\nset matcols :=");
	  for(int i=0; i<mat.M(); ++i) printf(" %d", i+1);
	  //	  printf("set matlins := 1..%ld;\nset matcols := 1..%ld;\n",mat.N(),mat.M());

	  printf(";\nparam mat:\n");
	  for(int i=0; i<mat.M(); ++i) printf("%d ", i+1);
	  printf(":=\n");
	  for(int i=0; i<mat.N(); ++i) {
	    printf("\n%d", i+1);
	    for(int j=0; j<mat.M(); ++j)
	      printf(" %.1lf", mat(i,j));
	  }
	  printf(" ;\n");
	  /*
	  printf("\nparam mat :=");
	  for(int i=0; i<mat.N(); ++i) {
	    for(int j=0; j<mat.M(); ++j) {
	      printf("\n[%d,%d]", i+1, j+1);
	      printf(" %.1lf", mat(i,j));
	    }
	  }
	  printf(" ;\n");
	  */
	  return 0;
	}
	  

	int rend = rini;
	if(argc>3) rend = atoi(argv[3]);
	int nloop = 10;
	if(argc>4) nloop = atoi(argv[4]);
	int nrand = randnum(1e8+1);
	if(argc>5)
		nrand = atoi(argv[5]);
	strcpy(dirname,"");
	if(argc>6)
		strcpy(dirname,argv[6]);
	else if(argv[1][0] != '/')
		strcpy(dirname,".");
	if(argc>8) {
		char *lline = new char[100000];
		Rprintf("reading file %s\n",argv[7]);
		FILE *fp1 = fopen(argv[7],"r");
		if(!fp1) { fprintf(stderr,"Error reading file %s\n",argv[7]); return 1; }
		origmat = matrix<>::fromcsv(fp1);
		fclose(fp1);
		myassert(origmat.min() >= 0.);
		origsum = 1./origmat.sum();
		comparetoorigmat = origmat;
		comparetoorigmat *= origsum;
		comparetoorigmat += 1.e-8;
		comparetoorigmat.applyloge();
		comparetoorigmat *= origmat;
		klorig = comparetoorigmat.sum() * origsum;

		comparetoorigmat = origmat;
		FILE *fp = fopen(argv[8],"r");
		if(!fp) {
			fprintf(stderr,"Error reading file %s\n",argv[8]);
			return 1;
		}
#define NMAXC 5010
		filtermat = matrix<>(mat.N(),NMAXC);
		for(int i=0; i<mat.N(); i++) {
			if(fgets(lline,99998,fp) <= 0) {
				fprintf(stderr,"Error reading file %s line %d\n",argv[8],i);
				return 1;
			}
			char *s = strtok(lline,",");
			int j=1;
			while(s && j<NMAXC) {
				filtermat(i,j) = atof(s);
				s = strtok(0,",");
				j++;
			}
			if(j==NMAXC) {
				fprintf(stderr,"Error reading file %s line %d (too many columns)\n",argv[8],i);
				return 1;
			}
			filtermat(i,0) = j;
		}
		filtermat = filtermat.submat(0,filtermat.N(),0,filtermat.max(-1,0)+1);
//		filtermat.print("FILTER MATRIX");
		delete [] lline;
		fclose(fp);
		if(argc>9) nmfold=1;
	}

	matrix<> *opt_w, *opt_h, *consensus, *opt_conn, **ww=0, **hh=0;
	if(nloop > 1 && nrand >= 0) {
		nmfconsensusCpp(mat, rini, rend, nloop, opt_w, opt_h, consensus, opt_conn, ww, hh, true, DIVERGENCE, "argmax", 0.1, 0.1, 0.1, 0, false);
//		for(int i=rini; i<=rend; i++) {
//		consensus[i].print("consensus");
//			opt_w[i].print("W");
//			opt_h[i].print("H");
//			(opt_w[i].multiply(opt_h[i])).print("W * H");
//		}
	} else if(nrand >= 0) {
		for(int i=rini; i<=rend; i++) {
			matrix<> w0, h0, w, h;
			Rprintf("run(%s,%s,%d,%d)\n", dirname, argv[1], i, nrand);
			nmfCpp(mat, i, true, DIVERGENCE, 0.1, 0.1, 0.1, 0, w0, h0, w, h);
			sprintf(sspace, "%s/%s-W-%02d-%08d", dirname, argv[1], i, nrand);
			Rprintf("writing %s\n", sspace);
			FILE *outf = fopen(sspace, "w");
			myassert(outf);
			w.tocsv(0, outf);
			fclose(outf);
			sprintf(sspace, "%s/%s-H-%02d-%08d", dirname, argv[1], i, nrand);
			Rprintf("writing %s\n", sspace);
			outf = fopen(sspace, "w");
			myassert(outf);
			h.tocsv(0, outf);
			fclose(outf);
			matrix<> atm((mat - w.multiply(h)));
			Rprintf("rmse = %lf\n", sqrt((atm * atm).sum())/sqrt(atm.size()));
		}
		return 0;
	} else if(nrand < -2) {
	  nmfconsensusCpp(mat, rini, rend, nloop, opt_w, opt_h, consensus, opt_conn,  ww, hh, true, DIVERGENCE, "argmax", 0.1, 0.1, 0.1, argv[1], false);

	  // validation of the NMF runs
	  double *err1 = new double[rend+1];
	  for(int i=rini; i<=rend; i++) err1[i] = 0.;

	  nrand = -nrand;
	  std::vector<int> looplist;
	  for(int i=0; i<nloop; i++) looplist.push_back(i+1);
	  std::random_shuffle(looplist.begin(),looplist.end());
	  for(int i=0; i< nrand; i++) {
	    int ini = (int) (nloop * (((double) i)/nrand)), fim = (int) (nloop * ((double) i+1.)/nrand);
	    matrix<> *opt_w1, *opt_h1, *consensus1, *opt_conn1;
	    nmfconsensusCppList(mat, rini, rend, fim%nloop, ini, looplist, opt_w1, opt_h1, consensus1, opt_conn1, ww, hh, true, DIVERGENCE, "argmax", 0.1, 0.1, 0.1,  argv[1], false, false);

	    for(int ii=rini; ii<=rend; ii++) {
	      double err=0.;
	      for(int k1=0; k1<consensus[ii].N(); k1++) {
		for(int k2=0; k2<consensus[ii].M(); k2++) {
		  err += (consensus[ii](k1,k2) - consensus1[ii](k1,k2))*(consensus[ii](k1,k2) - consensus1[ii](k1,k2));
		}
	      }
	      Rprintf("round=%d (from %d to %d), rank=%d, rmse=%.12lf\n", i, fim%nloop, ini, ii, sqrt((err/consensus[ii].N()) / consensus[ii].M()));
	      err1[ii] += sqrt((err/consensus[ii].N()) / consensus[ii].M());
	    }
	  }
	  for(int i=rini; i<=rend; i++) {
	    err1[i] /= nrand;
	    Rprintf("rank %d: average rmse = %.12lf\n", i, err1[i]);
	  }
	  delete [] err1;
	  return 0;
	} else {
	  if(nrand == -2) mat = mat.T();
	  nmfconsensusCpp(mat, rini, rend, nloop, opt_w, opt_h, consensus, opt_conn, ww, hh,  true, DIVERGENCE, "argmax", 0.1, 0.1, 0.1,argv[1], nrand==-2);
	}
	for(int i=rini; i<=rend; i++) {
		sprintf(sspace, "%s/%s-%02d-consensus", dirname, argv[1], i);
		Rprintf("writing %s\n", sspace);
		FILE *outf = fopen(sspace, "w");
		myassert(outf);
		consensus[i].tocsv(0, outf);
		fclose(outf);
		sprintf(sspace, "%s/%s-%02d-opt_w", dirname, argv[1], i);
		Rprintf("writing %s\n", sspace);
		outf = fopen(sspace, "w");
		myassert(outf);
		opt_w[i].tocsv(0, outf);
		fclose(outf);
		sprintf(sspace, "%s/%s-%02d-opt_h", dirname, argv[1], i);
		Rprintf("writing %s\n", sspace);
		outf = fopen(sspace, "w");
		myassert(outf);
		opt_h[i].tocsv(0, outf);
		fclose(outf);
		sprintf(sspace, "%s/%s-%02d-opt_conn", dirname, argv[1], i);
		Rprintf("writing %s\n", sspace);
		outf = fopen(sspace, "w");
		myassert(outf);
		opt_conn[i].tocsv(0, outf);
		fclose(outf);
	}
	return 0;
}
#endif
