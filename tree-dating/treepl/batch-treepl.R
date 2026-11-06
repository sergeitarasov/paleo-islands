library(ape)
library(MCMCtreeR)
library(phytools)
library(phangorn)
library(geiger)
library(dispRity)

setwd("treepl/")
source('R/make-treepl.r')


# Read tree and get node info
tree <- read.tree("data/mcmctree_53-tips.tre")
node_info <- get_node_info(tree)
node_info

# Read MCMC samples
samples <- read_mcmc_samples("../latest_MCMC/mcmc.txt")
dim(samples)

# Write batch of TreePL input files
write_treepl_batch(
  node_info = node_info,
  samples = samples,
  template_file = "data/treepl_temp_batch.ini",
  generations = 1000:2000,
  output_dir = "treepl_batch/input",
  treefile = "../data/sunken_run7_70p_raw-rooted.tre",
  outfile_pattern = "../timetrees/timetree_%d.tre"
)

#------------------
samples[1001,]

samples[["54"]]
samples[["54"]][1000:2000]
mean(samples[["54"]][1000:2000])


