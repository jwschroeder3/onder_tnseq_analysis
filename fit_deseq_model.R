#!/usr/bin/env Rscript
library(optparse)
library(tidyverse)
library(DESeq2)
library(GenomicAlignments)
library(Rsamtools)
library(BiocParallel)
library(configr)

# parse config file
args = commandArgs(trailingOnly=TRUE)
configs <- read.config( file=args[1] )

# load needed helper functions
source(file.path( configs$srcdir, "helpers.R"))

# set up parallelization
cores = configs$compute$cores
print(paste0("Using ", cores, " cores."))
register(MulticoreParam(workers=cores)) # bpprogressbar=true
 
# set up input/output locations
out_direc = configs$outdir
se_fname = file.path(
    out_direc,
    paste( configs$outprefix, "summarized_experiment.RData", sep="_")
)
dds_fname = file.path(
    out_direc,
    paste(configs$outprefix, "deseq_fit.RData", sep="_")
)

# read in conditions and summarized expreiments object (se)
se_var_name = load(se_fname)
#print(se_var_name)
#print(rowData(se))
#coldat=colData(se)
#print(coldat)
#count_mat = assays(se)[["counts"]]
#print(head(count_mat))
#all_zero = apply(count_mat==0, all, MAR=2)
#print(all_zero)
#print(coldat[all_zero,c("genotype","condition")])
#stop()
#assays(se)[["counts"]] = count_mat[,!all_zero]
#colData(se) = coldat[!all_zero,]
#se = se[,!all_zero]



#stop()
#my_data[["RNA"]]@counts <- as.matrix(my_data[["RNA"]]@counts)+1
dds = DESeqDataSet(
    se,
    design = as.formula(configs$deseq$design)
)

rownames(dds) = replace_na(mcols(dds)$locus_tag, "unknown")

dds = DESeq(
    dds,
    parallel=TRUE
)

save(dds, file=dds_fname)
