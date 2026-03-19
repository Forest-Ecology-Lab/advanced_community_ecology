# --------------------------------- Set up Script --------------------------- *
# Project setup for Advanced Community Ecology.
#
#IMPORTANT!!!!!!
# Run this script before running everything else in this project.
# This installs all required libraries and packages using renv.
# --------------------------------------------------------------------------- *

{
  # install renv if needed
  if (!requireNamespace("renv", quietly = TRUE)) {
    install.packages("renv")
  }

  # restore project library
  cat("Restoring project environment using renv...\n")
  renv::restore(prompt = FALSE)

  cat("Environment ready.\n")
}
