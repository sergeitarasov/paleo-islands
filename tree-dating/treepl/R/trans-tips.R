translate_tips <- function(tree, dictionary_file) {
  # Read the dictionary
  dict <- read.csv(dictionary_file,
                   header = FALSE,
                   sep = ";",
                   stringsAsFactors = FALSE,
                   col.names = c("new_name", "old_name"))  # Note the reversed column names

  # Get current tip labels
  old_tips <- tree$tip.label

  # Create new tips vector
  new_tips <- old_tips

  # Print debugging info
  cat("Dictionary entries:", nrow(dict), "\n")
  cat("Tree tips:", length(old_tips), "\n\n")

  # Translate tips
  for (i in seq_along(old_tips)) {
    # Check if tip is in old names (second column)
    match_idx <- match(old_tips[i], dict$old_name)

    if (!is.na(match_idx)) {
      new_tips[i] <- dict$new_name[match_idx]
    }
  }

  # Find untranslated tips
  untranslated <- old_tips[!old_tips %in% dict$old_name]
  if (length(untranslated) > 0) {
    cat("Untranslated tips:\n")
    for (tip in untranslated) {
      cat(sprintf("- %s\n", tip))
    }
  }

  # Update tree tip labels
  tree$tip.label <- new_tips

  return(tree)
}
