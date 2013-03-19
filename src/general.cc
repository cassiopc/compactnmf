/*
## general.cc (package compactnmf) 
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

#include "general.hh"

bool samesign(double a, double b) {
  return a*b>0;
}

int intcmp(const void *a, const void *b) {
  return *(int *)a - *(int *)b;
}
int doublecmp(const void *a, const void *b) {
  double d = *(double *)a - *(double *)b;
  return d<0? -1: (d>0? 1: 0);
}
