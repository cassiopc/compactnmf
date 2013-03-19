/*
## general.hh (package compactnmf) 
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

#ifndef GENERAL_HH
#define GENERAL_HH

#include <cstdlib>
#include <ctime>

#define DEBMACRO(x)
#define initrand() { struct timespec s; if( clock_gettime( CLOCK_REALTIME, &s) == -1 ) { perror( "clock gettime" ); exit(1); } srand(s.tv_nsec); }
#define randnum(x) ((int) (((double) x)*rand()/(RAND_MAX+1.0)))
#define eps 0.00000001
#define maxi(X,Y) ((X)>(Y)?(X):(Y))
#define mini(X,Y) ((X)<(Y)?(X):(Y))
#define atoi(x) strtol(x,(char **) 0, 10)
#define atof(x) strtod(x,(char **) 0)
bool samesign(double a, double b);
int intcmp(const void *a, const void *b);
int doublecmp(const void *a, const void *b);

#endif
