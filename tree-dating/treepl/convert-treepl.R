library(ape)
library(MCMCtreeR)
library(phytools)
library(phangorn)
library(geiger)
library(dispRity)

setwd("treepl/")
source('R/make-treepl.r')

tree <- read.tree("data/mcmctree_53-tips.tre")
plot(tree)
# Get node information
node_info <- get_node_info(tree)

# Assuming you have node_info from get_node_info() function
template_file <- "data/treepl_temp.ini"
output_file <- "treepl_run/prime.ini"

write_treepl_input(node_info, template_file, output_file)


