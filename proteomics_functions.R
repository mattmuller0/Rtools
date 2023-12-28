###########################################################################
#
#                            proteomics_functions
#
###########################################################################
# Author: Matthew Muller
# Date: 2023-12-28
# Script Name: proteomics_functions

#======================== LIBRARIES ========================#
library(lintr) #nolint
library(httpgd) #nolint
library(languageserver) #nolint
library(devtools)
library(tidyverse)
library(glue)
library(OlinkAnalyze)
library(clusterProfiler)
library(org.Hs.eg.db)
library(AnnotationDbi)

# LOAD FUNCTIONS
# space reserved for sourcing in functions
source_url('https://raw.githubusercontent.com/mattmuller0/Rtools/main/general_functions.R')
source_url('https://raw.githubusercontent.com/mattmuller0/Rtools/main/plotting_functions.R')

#======================== CODE ========================#






















#======================== END ========================#
writeLines(capture.output(sessionInfo()), file.path(outdir, "session.log"))
save.image(file.path(outdir, "image.RData"))
