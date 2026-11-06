
library(ape)
library(MCMCtreeR)
library(phytools)
library(phangorn)
library(geiger)

library(dispRity)

# Run from the tree-dating directory.
setwd("treepl/")
source('R/make-treepl.r')


# Read tree and get node info
tree <- read.tree("data/mcmctree_53-tips.tree")
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
  treefile = "../data/iqtree_316-tips-rooted.tree",
  outfile_pattern = "../timetrees/timetree_%d.tre"
)

#------------------
samples[1001,]

samples[["54"]]
samples[["54"]][1000:2000]
mean(samples[["54"]][1000:2000])

#-------
# samples[1000,]
# updated_info <- update_node_info_with_sample(node_info, samples, 1000)
# #write_treepl_input(node_info, template_file, output_file)
# updated_info

#------------------

