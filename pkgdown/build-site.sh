#!/usr/bin/env bash
set -euo pipefail

Rscript --vanilla pkgdown/mmes-options/build-mmes-options.R
Rscript --vanilla pkgdown/mmes-options/check-mmes-options.R

Rscript -e 'pkgdown::init_site(); pkgdown::build_home(preview = FALSE); pkgdown::build_reference(examples = FALSE, preview = FALSE)'

Rscript -e 'pkg <- pkgdown:::as_pkgdown("."); writeLines(as.character(pkg$vignettes$name))' |
  while IFS= read -r article; do
    Rscript -e 'pkgdown::build_article(commandArgs(trailingOnly = TRUE)[1])' "$article"
  done

Rscript -e 'pkgdown::build_site(preview = FALSE, install = FALSE, examples = FALSE, lazy = TRUE, new_process = FALSE)'