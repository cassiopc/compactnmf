## compactlines     Copyright (C) 2012 C. P. de Campos (cassiopc@acm.org)

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

###### function compactlines #######
## PARAMETERS ARE AS FOLLOWS:

## datamat: each line contain one probe (feature), each column contain
## one sample (patient)

## verbose: show number of selected probes every 100 analyzed

## equalbynormalaberrated:
##    TRUE=consider values to be different (so counted in the hamming
## distance) if one value is zero (normal) and the other is not
## (aberrated)
##    FALSE=consider a value to be different by usual equality

## tolerance: how much two values have to numerically differ to be
## considered different

## maxhammingdiff: the maximum number of differences such that two
## features are still considered to be put together

## neighborhood: if defined as a number, then only features that are
## within this distance (in lines) from a given feature in the matrix
## will be analyzed to be combined with it.

compactlines <- function(datamat,maxhammingdiff=0,equalbynormalaberrated=FALSE,tolerance=0.01,verbose=TRUE,neighborhood=NULL)
{
  ## number of lines in the datamat
  n = length(datamat[,1])
  ## indicator of lines still to be considered
  totest = rep(TRUE,n)
  ## keeplines will contain the array of line numbers that are to be kept in the end
  keeplines = rep(0,n)

  ## if equalbynormalaberrated, then it only matters if it is zero or not
  if(equalbynormalaberrated) {
    datamat1 = abs(datamat) >= tolerance
  }
  ## i and t run over the line numbers, but t steps as lines are kept,
  ## so it is used as an index for the array keeplines, while i steps over
  ## lines that are analyzed, so runs over datamat
  t <- 0
  ## loop over the lines (probes)
  for(i in 1:n) {
    if(totest[i]) {
      totest[i]=FALSE
      if(verbose) {
        if((t %% 100)==0) {
          cat(i)
          cat(',')
        }
      }
      ## line i will be kept
      t <- t + 1
      keeplines[t] <- i
      ## lin.i is a copy of the i-th line, already adapted to the comparison according to normal vs. aberrated, or usual diff
      if(!equalbynormalaberrated) {
        lin.i <- datamat[i,]
      }
      ## for each of the lines still not considered or not included in a group
      ## check if they "match" with line i
      if(!is.null(neighborhood)) {
        ## if neighborhood is set, only look to nearby lines
        inij = max(1,i-neighborhood)
        endj = min(i+neighborhood,n)
        jset = (inij-1) + which(totest[inij:endj])
      } else {
        ## otherwise look into all not already used lines
        jset = which(totest)
      }
      ## loop over lines
      for(j in jset) {
        ## find the hamming distance
        if(equalbynormalaberrated) {
          ## based on zero/non zero
          numberdiffs <- sum(datamat1[i,] != datamat1[j,])
        } else {
          ## or based on actual hamming distance with a tolerance
          numberdiffs <- sum(abs(lin.i-datamat[j,])>=tolerance)
        }
        ## if the distance is small enough, merge the lines
        if (numberdiffs <= maxhammingdiff) {
          datamat[i,] <- datamat[i,] + datamat[j,]
          totest[j]=FALSE
        }
      }
    }
  }
  return(datamat[keeplines[1:t],])
}
