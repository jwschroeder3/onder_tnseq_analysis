#!/usr/bin/env Rscript
# this script will just write the coefficient names that we have available
library(optparse)
library(tidyverse)
library(limma)
library(DESeq2)
library(configr)

# parse config file
args = commandArgs(trailingOnly=TRUE)
configs <- read.config( file=args[1] )
contrast.type.1 <- args[2]
contrast.type.2 <- args[3]
outname <- args[4]

# load needed helper functions
source(file.path( configs$srcdir, "helpers.R"))

# set up input/output locations
data_direc = configs$input$data_dir
out_direc = configs$outdir
deseq_fname = file.path(out_direc, paste( configs$outprefix, "deseq_fit.RData", sep="_"))
 
load(deseq_fname)
print(resultsNames(dds))
