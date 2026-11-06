get_node_info <- function(tree) {
  # Check if tree is of class phylo
  if (!inherits(tree, "phylo")) {
    stop("Input must be a phylo object")
  }

  # Get number of nodes
  n_nodes <- tree$Nnode
  n_tips <- length(tree$tip.label)

  # Get node ages using tree.age
  ages <- tree.age(tree, order = "past")

  # Create empty dataframe to store results
  results <- data.frame(
    node = integer(),
    taxon1 = character(),
    taxon2 = character(),
    age = numeric(),
    stringsAsFactors = FALSE
  )

  # For each internal node
  for (node in (n_tips + 1):(n_tips + n_nodes)) {
    # Get all descendant tips
    desc <- try(phangorn::Descendants(tree, node, type="tips"), silent = TRUE)

    if (inherits(desc, "try-error") || length(desc[[1]]) == 0) {
      warning(sprintf("No descendants found for node %d", node))
      next
    }

    # Get tip labels
    tips <- tree$tip.label[desc[[1]]]

    # Find the most distant pair of tips
    if (length(tips) >= 2) {
      dist_matrix <- cophenetic(tree)[tips, tips]
      max_dist <- which(dist_matrix == max(dist_matrix), arr.ind = TRUE)[1,]
      taxon1 <- tips[max_dist[1]]
      taxon2 <- tips[max_dist[2]]
    } else {
      taxon1 <- tips[1]
      taxon2 <- tips[1]
    }

    # Get node age from tree.age results
    node_age <- ages$ages[node - n_tips + n_tips]

    results <- rbind(results, data.frame(
      node = node,
      taxon1 = taxon1,
      taxon2 = taxon2,
      age = node_age,
      stringsAsFactors = FALSE
    ))
  }

  return(results)
}

# Example usage:
# library(ape)
# tree <- read.tree("path/to/tree.tre")
# node_info <- get_node_info(tree)
# print(node_info)


write_treepl_input <- function(node_info, template_file, output_file) {
  # Read template
  template <- readLines(template_file)

  # Create calibration lines
  calibrations <- character()
  for(i in 1:nrow(node_info)) {
    # Create node name using node number instead of taxa names
    node_name <- sprintf("node_%d", node_info$node[i])

    # Add mrca, min and max constraints
    mrca_line <- sprintf("mrca = %s %s %s",
                        node_name,
                        node_info$taxon1[i],
                        node_info$taxon2[i])
    min_line <- sprintf("min = %s %.4f",
                       node_name,
                       node_info$age[i])
    max_line <- sprintf("max = %s %.4f",
                       node_name,
                       node_info$age[i])

    # Add empty line between node definitions
    calibrations <- c(calibrations, mrca_line, min_line, max_line, "")
  }

  # Combine all parts
  output <- c(
    template,
    "",
    calibrations
  )

  # Write output file
  writeLines(output, output_file)
}

read_mcmc_samples <- function(filename = "mcmc.txt") {
  # Function to read MCMC samples and clean column names to numbers only
  #
  # Args:
  #   filename: Path to the MCMC samples file (default: "mcmc.txt")
  #
  # Returns:
  #   A data frame with cleaned column names (numbers only)

  if (!file.exists(filename)) {
    stop(paste("File not found:", filename))
  }

  tryCatch({
    # Read the data with headers
    samples <- read.table(filename, header = TRUE)

    if (nrow(samples) == 0) {
      stop("The file is empty")
    }

    # Get columns that start with "t_n"
    t_n_cols <- grep("^t_n", names(samples), value = TRUE)

    if (length(t_n_cols) == 0) {
      stop("No columns starting with 't_n' found")
    }

    # Extract only the t_n columns
    samples_subset <- samples[, t_n_cols]

    # Clean column names by removing "t_n" prefix
    new_names <- gsub("t_n", "", names(samples_subset))
    names(samples_subset) <- new_names

    return(samples_subset)

  }, error = function(e) {
    stop(paste("Error reading MCMC samples:", e$message))
  })
}


update_node_info_with_sample <- function(node_info, samples, generation) {
  # Function to update node ages with values from a specific MCMC generation
  #
  # Args:
  #   node_info: Node information data frame from get_node_info()
  #   samples: MCMC samples data frame from read_mcmc_samples()
  #   generation: Which MCMC generation to use (row number)
  #
  # Returns:
  #   Updated node_info data frame with ages from the specified generation

  # Create a copy of original node_info
  new_node_info <- node_info

  # Get ages from specified generation
  sample_ages <- samples[generation,]

  # Update ages in node_info
  for(col in names(sample_ages)) {
    # Find matching node in node_info
    idx <- which(node_info$node == col)
    if(length(idx) > 0) {
      new_node_info$age[idx] <- sample_ages[[col]]
    }
  }

  return(new_node_info)
}

write_treepl_batch <- function(node_info, samples, template_file,
                               generations = 1000:2000,
                               output_dir = "treepl_batch/input",
                               treefile = "data/mcmctree_316-tips-rooted.tree",
                               outfile_pattern = "timetrees/timetree_%d.tre") {
  # Function to write multiple TreePL input files for selected MCMC generations
  #
  # Args:
  #   node_info: Node information data frame from get_node_info()
  #   samples: MCMC samples from read_mcmc_samples()
  #   template_file: Template file for TreePL input
  #   generations: Vector of generation numbers to process
  #   output_dir: Directory for output files
  #   treefile: Path to the input tree file
  #   outfile_pattern: Pattern for output files, use %d for generation number

  # Create output directory if it doesn't exist
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

  # Process each generation
  for(gen in generations) {
    # Update node ages with samples from this generation
    updated_info <- update_node_info_with_sample(node_info, samples, gen)

    # Create output filename for TreePL input
    output_file <- file.path(output_dir, paste0("input_", gen, ".ini"))

    # Read template file
    template <- readLines(template_file)

    # Update treefile and outfile paths
    template[grep("^treefile", template)] <- paste("treefile = ", treefile)
    template[grep("^outfile", template)] <- paste("outfile = ", sprintf(outfile_pattern, gen))

    # Write modified template with updated node ages
    writeLines(template, output_file)

    # Write TreePL input file with node info
    write_treepl_input(updated_info, output_file, output_file)

    # Print progress
    cat(sprintf("Written file for generation %d\n", gen))
  }

  cat(sprintf("\nCompleted writing %d files\n", length(generations)))
}
