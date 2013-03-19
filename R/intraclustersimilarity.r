## intracluster.similarity     Copyright (C) 2012 C. P. de Campos (cassiopc@acm.org)

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

###### function intracluster.similarity #######
## PARAMETERS ARE AS FOLLOWS:

## consensusmat: similarity matrix between elements such that 0 means
##               no similarity and 1 total similarity this matrix is
##               supposed to be simmetric

## clusterid: an array with the cluster numbers that were assigned to
##            the elements

## RETURNS: an array with the intra cluster score of each cluster. The
##          numbers in the output come ordered by the cluster number
##          itself, so one has to do the mapping between score and the
##          cluster it relates to (1st number relates to the cluster
##          of lowest ID, and so on).

intracluster.similarity <- function(consensusmat,clusterid) {
  clusts <- sort(unique(clusterid))
  nclusts <- length(clusts)
  score <- rep(0,nclusts)
  ## for each cluster
  for(cl in 1:nclusts) {
      ## nr is the number of elements with the ID clusts[cl]
      nr <- sum(clusterid==clusts[cl]);
      ## compute the sum of the similarities of the elements that are with the given ID (subtracting nr to
      ## account for the pairs (i,i) in the table which are always equal to one, as the similarity to
      ## themselves is always perfect, and then this number is averaged
      score[cl] <- (sum(consensusmat[clusterid==clusts[cl],clusterid==clusts[cl]])-nr)/(nr*(nr-1))
  }
  ## later on, to achieve a single score, one might take the mean, median, or maximum value from the
  ## array score, for example.
  return(list(clusters=clusts,score=score))
}
