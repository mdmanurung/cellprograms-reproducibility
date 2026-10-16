#!/usr/bin/env Rscript
# Install the cellprograms analysis environment (R >= 4.3).
# Pinned exactly as in DESCRIPTION Remotes; see environment/sessionInfo.txt
# for the full as-run versions.

repos <- "https://cloud.r-project.org"

# 1. Core dependencies (Imports)
install.packages(c("cli", "digest", "ebnm", "remotes"), repos = repos)

# 2. Pinned backends (same pins as the package's CI)
remotes::install_github("willwerscheid/flashier@fd44811ab22d61ee905e1964175571c452fe7683")
remotes::install_github("mdmanurung/flashier.utils")

# 3. The package itself — pinned to the exact commit the analyses were run with
#    (repo: https://github.com/mdmanurung/cellprograms)
remotes::install_github("mdmanurung/cellprograms@413e290")

# 4. Analysis-side Suggests used by the vignettes/benchmarks
install.packages(c("ggplot2", "data.table", "igraph", "mashr", "testthat"),
                 repos = repos)

# 5. Sanity check
stopifnot(packageVersion("cellprograms") >= "0.2.0")
cat("cellprograms environment ready:\n")
print(packageVersion("cellprograms"))
