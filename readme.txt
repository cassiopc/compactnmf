## compactnmf package
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

TO COMPILE THE STANDALONE VERSION:

$ cd src
$ CPURE=1 make

TO COMPILE THE STANDALONE VERSION WITH CBLAS SUPPORT:
(CBLAS is a fast matrix multiplication library available under
 ubuntu linux by the package libblas-dev)

$ cd src
$ CPURE=1 CBLAS=1 make

TO COMPILE THE LIBRARY FOR R:

$ cd src
$ make

TO COMPILE THE LIBRARY FOR R WITH CBLAS:

$ cd src
$ CBLAS=1 make

The code is built into the directory bin/ For a quick help of the
stand-alone version, just type its command without arguments. For use
within R, it is necessary to place the file nmf.so in the appropriate
directory. The R code which calls the library is at the folder R/.
Please refer to the R source code for more details, it is documented.
Any questions, do not hesitate to contact me: cassiopc@acm.org
