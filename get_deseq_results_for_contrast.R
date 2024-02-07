#!/usr/bin/env Rscript
# this script will output the results of a specified comparison from our deseq fit
library(optparse)
library(tidyverse)
library(limma)
library(DESeq2)
library(configr)
library(IHW)

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

print(
    paste0("Calculating log2(FC_", contrast.type.1, ") - log2(FC_", contrast.type.2, ")")
)

if (contrast.type.2 == "") {
    possible_contrasts = resultsNames(dds)
    #print(possible_contrasts)
    cont = as.integer(possible_contrasts == contrast.type.1)
    #print(contrast.type.1)
    #print(cont)
    #stop()
    res.IHW <- results(
        dds,
        contrast=cont,
        filterFun=ihw
    )
    res.LFC <- lfcShrink(
        dds,
        contrast=cont,
        type="ashr"
    )
    resmat <- results(
        dds,
        contrast=cont
    )
} else {
    res.IHW <- results(
        dds,
        contrast=list(c(contrast.type.1), c(contrast.type.2)),
        filterFun=ihw
    )
    res.LFC <- lfcShrink(
        dds,
        contrast=list(c(contrast.type.1), c(contrast.type.2)),
        type="ashr"
    )
    resmat <- results(
        dds,
        contrast=list(c(contrast.type.1), c(contrast.type.2))
    )
}
resmat$padj.IHW = res.IHW$padj
resmat$log2fc.shrunk = res.LFC$log2FoldChange
out_fname = file.path(out_direc, outname)
print(paste0("Writing results to ", out_fname))
write.csv( resmat[order(resmat$padj.IHW),], out_fname, quote=FALSE)

