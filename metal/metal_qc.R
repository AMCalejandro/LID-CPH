#!/bin/Rscript


#--- Setting up the libpath, loading necesary libraries, and get the META file to QC
args<-commandArgs(trailingOnly=TRUE)
library.path <- .libPaths()
library("tidyr", lib.loc = library.path)
library("dplyr", lib.loc = library.path)
library("stringr", lib.loc = library.path)
library("data.table", lib.loc = library.path)

# Load utils
source("../utils/utils.R")

# Load data
meta_df = fread(args[1])

# Run the metal harmonise function
mata_harmonised = qc_metal(meta_df, N = 5)

# Show the header of the QCed file on stdout
head(meta_harmonised)


print("QC DONE")