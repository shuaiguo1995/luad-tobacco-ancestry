# Runner: executes an Rmd's code chunks without pandoc.
#   Rscript run.R 01_four_groups.Rmd
# Run from inside scripts/ so that the relative paths in each Rmd resolve.
args = commandArgs(trailingOnly = TRUE)
target = args[1]
cat("=== running", target, "===\n")
p = file.path(tempdir(), paste0(sub("\\.Rmd$", "", target), "_chunks.R"))
knitr::purl(target, output = p, quiet = TRUE)
source(p, echo = FALSE)
cat("=== finished", target, "===\n")
