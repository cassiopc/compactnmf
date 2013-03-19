## compactnmf     Copyright (C) 2012 C. P. de Campos (cassiopc@acm.org)

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

###### function compactnmf #######
## NON-NEGATIVE MATRIX FACTORIZATION
## PARAMETERS ARE AS FOLLOWS:
## v    : matrix to be factorized
## r    : rank of the factorization
## w    : N x r (w0 is an initial guess for w)
## h    : r x M (h0 is an initial guess for h)
## mu   : update step size
## lambda : regularization for h (if applicable)
## gamma  : regularization for w (if applicable)
compactnmf <- function(v, r, verbose=0, method="divergence", lambda=0.1,
                       mu=0.1, gamma=0.1, w0=NULL, h0=NULL )
{
  dyn.load('nmf.so');
  ## it's easier and faster to deal with numbers
  if (method=="euclidean") {
    intmethod = 1
  }
  else if (method=="divergence") {
    intmethod = 2
  }
  else if (method=="modified") {
    intmethod = 3
  }
  else if (method=="local") {
    intmethod = 4
  }
  else if (method=="nnsc") {
    intmethod = 5
  }
  else if (method=="nnsc2") {
    intmethod = 6
  }
  else {
    error("Invalid method")
  }
  if(!is.null(w0)) w0 = as.matrix(w0);
  if(!is.null(h0)) h0 = as.matrix(h0);
  
  v = .Call('nmfCpp_wrapper',as.matrix(v),as.integer(r),as.integer(verbose),as.integer(intmethod),as.real(lambda),
    as.real(mu),as.real(gamma),as.integer(0),w0,h0) 
  res = c()
  res$w = v$w
  res$h = v$h
  return(res)
}
