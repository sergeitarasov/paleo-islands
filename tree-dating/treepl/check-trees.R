library(ape)
library(MCMCtreeR)
library(phytools)
library(phangorn)
library(geiger)
library(dplyr)

library(dispRity)

setwd("treepl/")
source('R/make-treepl.r')
source('R/trans-tips.r')

tr_mc <-  read.tree("data/reduced_tree_run1.tre")[[2]]
write.tree(tr_mc, "data/mcmctree_53-tips.tre")


