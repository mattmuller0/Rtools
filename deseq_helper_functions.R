###########################################################################
#
#                                 HEADER
#
###########################################################################
# Author: Matthew Muller
# 
# Date: 2023-02-18
# 
# Script Name: DESeq2 Helper Functions
# 
# Notes:
# This script contains helper functions that can be used to assist in differential
# expression analysis using DESeq2. Most inputs are expected as DDS or SE objects.

message("This script is being deprecated. Please use the rna_functions and enrichment_functions scripts instead.")

###########################################################################
#
#                                 LIBRARIES
#
###########################################################################
packages <- c(
  "tidyverse",
  "ggplot2",
  "devtools",
  "BiocManager",
  "SummarizedExperiment",
  "DESeq2",
  "edgeR",
  "AnnotationDbi",
  "apeglm",
  "ggpubr",
  "EnhancedVolcano",
  "cowplot",
  "ggtree",
  "ggrepel",
  "ggbiplot"
)

for (pkg in packages) {suppressPackageStartupMessages(library(pkg, character.only = T, quietly = T))} # nolint

# LOAD FUNCTIONS
# space reserved for sourcing in functions
options(stringsAsFactors = FALSE)

source_url('https://raw.githubusercontent.com/mattmuller0/Rtools/main/converting_functions.R')
source_url('https://raw.githubusercontent.com/mattmuller0/Rtools/main/plotting_functions.R')
source_url('https://raw.githubusercontent.com/mattmuller0/Rtools/main/stats_functions.R')

###########################################################################
#
#                                 CODE
#
###########################################################################

# the code has been moved to the filtering_functions, rna_functions, and enrichment_functions scripts
source_url('https://raw.githubusercontent.com/mattmuller0/Rtools/main/filtering_functions.R')
source_url('https://raw.githubusercontent.com/mattmuller0/Rtools/main/rna_functions.R')
source_url('https://raw.githubusercontent.com/mattmuller0/Rtools/main/enrichment_functions.R')