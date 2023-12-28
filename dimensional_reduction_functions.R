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
# Add code here
#

# Function to run pca on a dds object
# Arguments:
#   - dds: dds object
#   - group: character, column name of the group to use
