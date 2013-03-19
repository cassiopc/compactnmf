/*
## nmf.cc (package compactnmf) 
## Copyright (C) 2012 Ivo Kwee and C. P. de Campos (cassiopc@acm.org)

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

#include "nmf.hh"
extern matrix<> origmat,filtermat,comparetoorigmat;
extern int nmfold;
extern double origsum,klorig;
static const double EPS = 1e-8;

void nmfCpp(matrix<> &v, int r, int verbose, int method, double lambda, double mu,
			double gamma, int reorder, matrix<> &w0, matrix<> &h0, matrix<> &w, matrix<> &h) {
	if (verbose) {
		Rprintf("method=%d\n",method);
	}

	//## test for negative values in v
	if ( v.min() < 0 ){
		error("*error* matrix entries can not be negative\n");
	}

	if ( (v.apply(LINEWISE,SUM)).min() == 0 ){
		error("*error* not all entries cain a row can be zero\n");
	}

	//## Modified algorithm needs to work on transpose!
	if (method==MODIFIED) {
		v = v.T();
		matrix<> w0_old(w0);
		if (!w0.empty()) w0 = h0.T();
		if (!h0.empty()) h0 = w0_old.T();
	}

	//## initialize program parameters
	const tint n = v.N();
	const tint m = v.M();
	const int stopconv = 10; //40
	const int niter = 3000;
  
	matrix<int> cons(m,m,0);
	matrix<int> consold = cons;
	int inc = 0;
  
	//## initialize random w and h
	if (w0.empty()) {
		w0.realoc(n, r);
        //runif(n*r);
		for(int i=0; i<w0.N(); i++)
			for(int j=0; j<w0.M(); j++)
				w0(i,j) = randnum(1e6+1)/1.e6; 
	}
	if (h0.empty()) {
		h0.realoc(r, m);
        //runif(r*m);
		for(int i=0; i<h0.N(); i++)
			for(int j=0; j<h0.M(); j++)
				h0(i,j) = randnum(1e6+1)/1.e6; 
	}

	w = w0;
	h = h0;
	w0.clear();
	h0.clear();

	//## start loop
	for(int it=1; it<=niter; it++) {
		switch(method) {
			case EUCLIDEAN:
				//## Euclidean distance reducing NMF iterations (Lee-Seung)
				//h  = h * ( t(w) %*% v )/ ( t(w) %*% w %*% h );
				h *= ( w.T().multiply(v) ) / ( (w.T().multiply(w)).multiply(h) ); 
				//w  = w * ( (v %*% t(h)) / (w %*% h %*% t(h)) );
				w *= ( v.multiply(h.T()) ) / ( (w.multiply(h)).multiply(h.T()) );
				break;
			case KMEANS:
				//## Frobeneus distance using H*H^t (Ding-He-Simon). Works as a relaxed k-means
                //H ← max (V H (H^t H )^−1,0)
				if(it==1) {
					w = v.multiply(w).multiply(w.T().multiply(w).inverse());
					for(int i=0; i<w.N(); i++) 
						for(int j=0; j<w.M();j++)
							if(w(i,j)<EPS) w(i,j) = EPS;
				}
				else {
					double beta=.5;
					h = v.multiply(w);
					matrix<> tmp = w.multiply(w.T()).multiply(w);
					for(int i=0; i<w.N(); i++) 
						for(int j=0; j<w.M();j++) {
							w(i,j) *= 1-beta+beta*(h(i,j)/tmp(i,j));
							if(w(i,j)<EPS) w(i,j) = EPS;
						}
				}
				break;
			case DIVERGENCE:
				//## Entropy-based divergence-reducing NMF iterations (Lee-Seung)
				//x1 = apply( w, 2, sum )
				//h  = h * ( t(w) %*% (v / (w %*% h))) / x1
				h *= ( (w.T()).multiply(v / w.multiply(h)) ) / w.apply(COLUMNWISE, SUM);

				//w  = w * ( (v / (w %*% h)) %*% t(h) ) 
				w *= (v / w.multiply(h)).multiply(h.T());
				
				//x2 = apply( h, 1, sum )
				//w  = t ( t(w) / x2 )
				w /= h.apply(LINEWISE,SUM);
				break;
			case MODIFIED:
				//## modified NNSC algorithm (see Badea and Tilivea,
                //## Pac. Symposium on Biocomputing, 2005) 
				//h  = h * ( t(w) %*% v ) / ( t(w) %*% w %*% h + lambda * h )
				h *= ( w.T().multiply(v) ) / ( (w.T().multiply(w)).multiply(h) + h * lambda);
				//w  = w + mu * (v - w %*% h) %*% t(h)
				w += (v - w.multiply(h) ).multiply(h.T()) * mu;

                //x2 = sqrt( apply(w*w,2,sum) )
                //w  = t( t(w) / x2 )
				w /= (w*w).apply(COLUMNWISE,SUM).T().sqrt();
				break;
			case LOCAL:
				//## Local NMF (Feng et al. 2002)
				//h  = sqrt( h * ( t(w) %*% (v / (w %*% h)))  )
				h = (h * w.T().multiply(v / w.multiply(h)) ).sqrt();
				//w  = w * ( (v / (w %*% h)) %*% t(h) ) 
				w *= (v / w.multiply(h)).multiply(h.T());
				//x1 = apply( w, 2, sum )
				//w  = t ( t(w) / x1 )
				w /= w.apply(COLUMNWISE,SUM).T();
				break;
			case NNSC:
                //## NNSC algorithm (Hoyer)
				//w  = w - mu * (w %*% h - v) %*% t(h)
				w -= (w.multiply(h) - v).multiply(h.T()) * mu;

				//w  = pmax(w,EPS)
				for(int i=0; i<w.N(); i++) 
					for(int j=0; j<w.M();j++)
						if(w(i,j)<EPS) w(i,j) = EPS;

				//x2 = sqrt( colSums(w*w) )
				//w  = t( t(w) / x2 )
				w /= (w*w).apply(COLUMNWISE,SUM).T().sqrt();
				//h  = h * ( t(w) %*% v ) / ( t(w) %*% w %*% h + lambda * h )
				h *= w.T().multiply(v) / ( (w.T().multiply(w)).multiply(h) + h * lambda);
				break;
			case NNSC2:
				//## NNSC algorithm (Shang-Zheng-Sun): double penalty on |w| and |h|.
				//w  = (1 - gamma) * w + mu * (v - w %*% h) %*% t(h)
				w = w * (1-gamma) + (v - w.multiply(h)).multiply(h.T()) * mu;
				//w  = pmax(w,EPS)
				for(int i=0; i<w.N(); i++) 
					for(int j=0; j<w.M();j++)
						if(w(i,j)<EPS) w(i,j) = EPS;

				//x2 = sqrt( colSums(w*w) )
				//w  = t( t(w) / x2 )
				w /= (w*w).apply(COLUMNWISE,SUM).T().sqrt();
				//h  = h * ( t(w) %*% v ) / ( t(w) %*% w %*% h + lambda * h )
				h *= w.T().multiply(v) / ( (w.T().multiply(w)).multiply(h) + h * lambda);
				break;
			default:
				error("nmf::ERROR:unknown algorithm\n");
		}
        
		//## test convergence every 10 iterations
		if ( (it%10) == 0 ) {
			//## adjust small values to avoid underflow
			//w = pmax(w,EPS)
			for(int i=0; i<w.N(); i++) 
				for(int j=0; j<w.M();j++)
					if(w(i,j)<EPS) w(i,j) = EPS;

			if(method == KMEANS) 
				h = w.T();
			else {
				//h = pmax(h,EPS)
				for(int i=0; i<h.N(); i++) 
					for(int j=0; j<h.M();j++)
						if(h(i,j)<EPS) h(i,j) = EPS;
			}

			//## construct connectivity matrix
			//index = apply( h, 2, which.max )
			matrix<> index = h.apply(COLUMNWISE,WHICHMAX);
			//mat1  = matrix( rep(index,m),m ,m, byrow=TRUE )			
			//mat2  = matrix( rep(index,m),m ,m, byrow=FALSE)
			//cons  = (mat1 == mat2)
			matrix<int> cons(m,m,0);
			for(int i=0; i<m; i++)
				for(int j=0; j<m; j++)
					if(index(i)==index(j)) cons(i,j)=1;

			if(cons == consold) inc++;
			else inc = 0;

			if (verbose) Rprintf("%d\t%d\t%d\n",it,inc,(int)(cons-consold).nnonzero());
			if (inc>stopconv) break;

			consold = cons;
		}
	}
	
	//## Modified algorithm needs to work on transpose!
	if (method==MODIFIED) {
		matrix<> tmp = h;
		h = w.T();
		w = tmp.T();
		v = v.T();
	}
	if(method == KMEANS)
		h = w.T();
}


matrix<> nmfconnectivityCpp(const matrix<> &h, std::string method="argmax") {
	int m = h.M();
	matrix<> conn(m,m,0);

	if (method == "argmax") {
        //    ## Compute m x m matrix which is 1 if samples are together, 0
		//    ## elsewhere. Determine sample assignment by its largest metagene
        //    ## expression value.
		//index = apply( h, 2, which.max )
		matrix<> index = h.apply(COLUMNWISE, WHICHMAX);
		//mat1  = matrix( rep(index,m),m ,m, byrow=TRUE )
		//mat2  = matrix( rep(index,m),m ,m, byrow=FALSE)
		//conn  = (mat1 == mat2)
		for(int i=0; i<m; i++)
			for(int j=0; j<m; j++)
				if(((int) (index(i)+0.5)) == ((int) (index(j)+0.5))) conn(i,j)=1;
	}
	else if (method == "average") {
		//hh = t( t(h) / apply(h,2,sum) )  ## normalize
		matrix<> hh = (h.T() / h.apply(COLUMNWISE, SUM)).T();
		matrix<> index = h.apply(COLUMNWISE, WHICHMAX);
		//mat1  = matrix( rep(index,m),m ,m, byrow=TRUE )
		//mat2  = matrix( rep(index,m),m ,m, byrow=FALSE)
		//conn  = (mat1 == mat2)
		for(int i=0; i<m; i++)
			for(int j=i+1; j<m; j++) {
				int ii = (int) (index(i)+.5);
				int jj = (int) (index(j)+.5);
				if(ii==jj) conn(i,j) = conn(j,i) = 1;
				else {
					double demi = hh(i,ii), demj = hh(j,jj);
					if(demi<EPS) demi=EPS;
					if(demj<EPS) demj=EPS;
					conn(i,j) = (hh(i,jj)/demi + hh(j,ii)/demj)/2.;
					conn(i,j) = conn(j,i) = conn(i,j)*conn(i,j);
				}
			}
	}
	else if (method == "covariance") {
        //    ## normalized covariance
		//hh = t( t(h) / apply(h,2,sum) )  ## normalize
		matrix<> hh = (h.T() / h.apply(COLUMNWISE, SUM)).T();
		//conn = t(hh) %*% (hh);
		conn = hh.T().multiply(hh);
	}
	else {
		error("nmfconnectivity::*error*:unknown method\n");
	}
	return conn;
}
  
/*
## nmfconsensu
## compute mean connectivity matrix
*/
void nmfconsensusCpp(matrix<> &a, int kstart, int kend, int nloop, 
					 matrix<>*&opt_w, matrix<>*&opt_h, matrix<>*&consensus, matrix<>*&opt_conn, matrix<>**&ww, matrix<>**&hh,
					 int verbose=false, int intmethod=DIVERGENCE, std::string conn_method="argmax",
					 double lambda=0.1, double mu=0.1, double gamma=0.1,  char *inputfile=0, bool transp=false) {
  std::vector<int> looplist;
  for(int i=0; i<nloop; i++) looplist.push_back(i+1);
  nmfconsensusCppList(a, kstart, kend, 0, nloop, looplist, opt_w, opt_h, consensus, opt_conn, ww, hh, verbose, intmethod, conn_method, lambda, mu, gamma, inputfile, true, transp);
}

void nmfconsensusCppList(matrix<> &a, int kstart, int kend, int iniloop, int fimloop, std::vector<int> looplist,
			 matrix<>*&opt_w, matrix<>*&opt_h, matrix<>*&consensus, matrix<>*&opt_conn, matrix<>**&ww, matrix<>**&hh,
			 int verbose=false, int intmethod=DIVERGENCE, std::string conn_method="argmax",
			 double lambda=0.1, double mu=0.1, double gamma=0.1,  char *inputfile=0, bool nopt=true, bool transp=false) {
    char sspace[200];

	//## test for negative values in v
	if ( a.min() < 0 ) {
		error("*error* matrix entries can not be negative\n");
		return;
	}

	//if ( min( apply(a,1,sum) ) == 0 ){
	if(a.apply(LINEWISE,SUM).min() == 0) {
		error("*error* not all entries of a row can be zero\n");
		return;
	}

	int m = a.M();

	if(inputfile && (ww==0 || hh==0)) {
	  nmfconsensusCppReadFiles(a, kstart, kend, iniloop, fimloop, looplist, ww, hh, true, inputfile, transp);
	}

    //  opt.w     = rep( list(NA), kend )
	opt_w = new matrix<>[kend+1];
    //  opt.h     = rep( list(NA), kend )
	opt_h = new matrix<>[kend+1];
    //  consensus = rep( list(matrix(0,m,m)), kend )
	consensus = new matrix<>[kend+1];
    //  opt.conn  = rep( list(matrix(0,m,m)), kend )
	opt_conn = new matrix<>[kend+1];
    //  conn      = matrix(0,m,m)
	matrix<> conn(m,m,0);
	matrix<> w0, h0, w, h;
	double asum = 1./a.sum();
	matrix<> masum = a*asum;
    double kl1 = (masum * masum.loge()).sum();
	matrix<> meps(a.N(),a.M(),1.e-8);
	int nloop = looplist.size();
	for (int j=kstart; j<=kend; j++) {
		if (verbose) Rprintf("\nrank %d\n",j);

		matrix<> connac(m,m,0), connacold(m,m,0);
		double min_err = 1e12;

		FILE *fp=0;
		if(inputfile) {
			sprintf(sspace, "%s-%s-%02d", inputfile, "E", j);
			fp = fopen(sspace, "wt");
			if(!fp) { Rprintf("error opening %s\n", sspace); exit(1); }
		}
		int ilooplist=iniloop;
		int qtd=0;
		do {
		  int iloop = looplist[ilooplist];
		  if (verbose) Rprintf(" iter%d ",iloop);
		  if(inputfile == 0) {
			w0.clear();
			h0.clear();
			w.clear();
			h.clear();
			nmfCpp(a, j, verbose, intmethod, lambda, mu, gamma, 0, w0, h0, w, h);
		  }
		  else {
			  w = ww[j][iloop];
			  h = hh[j][iloop];
		  }
		  double dv, dv2, kl, kl2;
		  if(nopt || qtd==0) {
			  matrix<> tmp1 = w.multiply(h);
			  double dvA = 0., dv2A = 0., klori2=0.;
			  double sum1 = 0.;
			  if(filtermat.size() > 0) {
				  comparetoorigmat.fill(-1.e8);
				  for(int l=0; l<tmp1.N(); l++) {
					  int nl = filtermat(l,0);
					  matrix<> curl = tmp1.line(l) / (double)(nmfold?1.:(nl-1.));
					  for(int ll=1; ll<nl; ll++) {
						  comparetoorigmat.updateline(filtermat(l,ll)-1,curl);
					  }
				  }
				  myassert(comparetoorigmat.min() >= 0.);
				  comparetoorigmat -= origmat;
				  dvA = comparetoorigmat.abs().sum()/(comparetoorigmat.N()*comparetoorigmat.M());
				  dv2A = comparetoorigmat.meansquare();
				  comparetoorigmat += origmat;

				  sum1 = 1./comparetoorigmat.sum();
				  comparetoorigmat *= sum1;
				  comparetoorigmat += 1.e-8;
				  comparetoorigmat.applyloge();
				  comparetoorigmat *= origmat;			  
				  klori2 = klorig - comparetoorigmat.sum() * origsum;
			  }
			  sum1 = 1./tmp1.sum();
			  matrix<> tmp = a - tmp1;
			  dv = tmp.abs().sum();
			  dv2 = tmp.meansquare();
			  kl = kl1 - (masum * (tmp1*sum1 + meps).loge()).sum();
			  kl2 = -(masum * (tmp1 + meps).loge()).sum();		  
			  if(inputfile)
				  fprintf(fp,"%d %d %lf %lf %lf %lf %lf %lf %lf %lf",qtd,iloop,dv,dv/(a.M()*a.N()),dv2,kl,kl2,dvA,dv2A,klori2);
		  }
			//      ##conn = r$c
			conn   = nmfconnectivityCpp(h, conn_method);
			connac += conn;

			//   ##err = sum(abs(a - w %*% h))
			if(nopt || qtd==0) {
			  double err = dv2;
			  if (err < min_err ) {
			    min_err = err;
			    opt_h[j] = h;
			    opt_w[j] = w;
			    opt_conn[j] = conn;
			  }
			  if(!nopt) {
			    min_err = 1e100;
			    opt_h[j].fill(0);
			    opt_w[j].fill(0);
			    opt_conn[j].fill(0);
			  }
			}
			ilooplist=(ilooplist+1)%nloop;
			qtd++;
			if(inputfile) {
				if(qtd > 1) {
					matrix<> tmp = connac / ((double) qtd) - connacold;
					dv = tmp.abs().sum();
					dv2 = tmp.meansquare();
					double dv3=dv/(tmp.M()*tmp.N());
					fprintf(fp," %lf %.8lf %.8lf %.8lf\n", dv, dv3, dv2, sqrt(dv2-dv3*dv3));
				} else
					fprintf(fp," 0 0 0 0\n");
				fflush(fp);
				connacold = connac / ((double) qtd);
			}
		} while (ilooplist!=(fimloop%nloop));
		if(inputfile)
			fclose(fp);

		consensus[j] = connac / ((double) qtd);
	}
}

void nmfconsensusCppReadFiles(matrix<> &a, int kstart, int kend, int iniloop, int fimloop, std::vector<int> looplist, matrix<>**&w, matrix<>**&h, int verbose, char *inputfile, bool transp) {
	//## test for negative values in v
	if ( a.min() < 0 ) {
		error("*error* matrix entries can not be negative\n");
		return;
	}

	//if ( min( apply(a,1,sum) ) == 0 ){
	if(a.apply(LINEWISE,SUM).min() == 0) {
		error("*error* not all entries of a row can be zero\n");
		return;
	}
	Rprintf("READING FILES WITH MATRICES W AND H\n");
	int m = a.M();
	int nn = a.N();
	
	if(transp) {
		m = a.N();
		nn= a.M();
	}
    //  opt.w     = rep( list(NA), kend )
	w = new matrix<>*[kend+1];
    //  opt.h     = rep( list(NA), kend )
	h = new matrix<>*[kend+1];
	int nloop = looplist.size();
	for (int j=kstart; j<=kend; j++) {
	  w[j] = new matrix<>[nloop+1];
	  h[j] = new matrix<>[nloop+1];
	  if (verbose) Rprintf("rank %d\n",j);
	  int ilooplist=iniloop;
	  do {
	    int iloop = looplist[ilooplist];
	    if (verbose) Rprintf(" iteration %d\n",iloop);
	    char sspace[200];
	    sprintf(sspace, "%s-%s-%02d-%08d", inputfile, "W", j, iloop);
	    Rprintf("reading %s\n",sspace);
	    FILE *fp = fopen(sspace, "r");
	    if(!fp) { Rprintf("error opening %s\n", sspace); exit(1); }
		w[j][iloop] = matrix<>::fromcsv(fp, nn, j, true,  false);
	    fclose(fp);
	    sprintf(sspace, "%s-%s-%02d-%08d", inputfile, "H", j, iloop);
	    Rprintf("reading %s\n",sspace);
	    fp = fopen(sspace, "r");
	    if(!fp) { Rprintf("error opening %s\n", sspace); exit(1); }
            h[j][iloop] = matrix<>::fromcsv(fp, j, m, true, false);
	    fclose(fp);

            if(transp) {
               matrix<> tmp = h[j][iloop].T();
               h[j][iloop] = w[j][iloop].T();
               w[j][iloop] = tmp;
            }

	    ilooplist=(ilooplist+1)%nloop;
	  } while (ilooplist!=(fimloop%nloop));
	}
}


#ifndef CPURE
extern "C" {
	
SEXP nmfCpp_wrapper(SEXP v, SEXP r, SEXP verbose, SEXP method, SEXP lambda, 
					SEXP mu, SEXP gamma, SEXP reorder, SEXP w0, SEXP h0)
{
	matrix<> a, b, c(v), w0_, h0_;
	if(w0 != NILSXP)
		w0_.set(w0);
	if(h0 != NILSXP)
		h0_.set(h0);
		
	nmfCpp(c,INTEGER(r)[0],INTEGER(verbose)[0],INTEGER(method)[0],
		   REAL(lambda)[0],REAL(mu)[0],REAL(gamma)[0],INTEGER(reorder)[0], w0_,h0_, a, b);
	
	SEXP w;
	PROTECT(w = allocMatrix(REALSXP, a.N(), a.M()));
	SEXP dimv = getAttrib(w, R_DimSymbol);
	for(int i=0; i<a.N(); i++)
		for(int j=0; j<a.M(); j++)
			getp(w,dimv,i,j) = a(i,j);
	
	SEXP h;
	PROTECT(h = allocMatrix(REALSXP, b.N(), b.M()));
	dimv = getAttrib(h, R_DimSymbol);
	for(int i=0; i<b.N(); i++)
		for(int j=0; j<b.M(); j++)
			getp(h,dimv,i,j) = b(i,j);
	
    SEXP rl = PROTECT(allocVector(VECSXP,2));
    SEXP nm = PROTECT(allocVector(STRSXP,2));
	SET_VECTOR_ELT(rl, 0, w);
	SET_STRING_ELT(nm, 0, mkChar("w"));
	SET_VECTOR_ELT(rl, 1, h);
	SET_STRING_ELT(nm, 1, mkChar("h"));
	
    setAttrib(rl, R_NamesSymbol, nm);
    UNPROTECT(4);
    return rl;
}
}
#endif
