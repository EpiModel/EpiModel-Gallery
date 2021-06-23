
# 1. renv package management
if (dir.exists("renv/")) {
  source("renv/activate.R")
  if (interactive()) {
    renv::status()
  }
} else {
  cat("* Run renv::init() to install the R packages for this project\n")
}

# 2. directory structure
if (!dir.exists("data/input/")) {
  dir.create("data/input/", recursive = TRUE)
}
if (!dir.exists("data/output/")) {
  dir.create("data/output/")
}
if (!dir.exists("out/")) {
  dir.create("out/")
}

# 3. Helpful aliases
rs <- function() .rs.restartR()

# 4. Standard options
options(deparse.max.lines = 5)
