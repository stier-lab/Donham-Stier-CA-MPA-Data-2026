# =============================================================================
# 00_libraries.R  (Data Processing Repo)
# =============================================================================
#
# PURPOSE:
#   Load R packages required for data processing (scripts 03-07).
#   This is a subset of the analysis repo's libraries. Only packages
#   needed for data import, wrangling, and response ratio calculation.
#
# AUTHORS: Emily Donham & Adrian Stier
# PROJECT: CA MPA Kelp Forest. Data Processing
# =============================================================================

# Load plyr BEFORE dplyr/tidyverse (masking order matters)
library(plyr)

# Core data manipulation
library(tidyverse)
library(here)
library(dplyr)  # Re-load so dplyr masks plyr

# Summary statistics (Rmisc::summarySE used in 07)
library(Rmisc)

# Spatial (sf used in data import for MPA boundaries)
library(sf)

# Confirmation
cat("Libraries loaded (data-processing subset).\n")
