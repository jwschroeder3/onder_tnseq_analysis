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
source(file.path( configs$srcdir, "/corexfs/schroedj/src/rnap_chip_analysis/src/helpers.R"))

# set up parallelization
cores = configs$compute$cores
print(paste0("Using ", cores, " cores."))
register(MulticoreParam(workers=cores))

# set up input/output locations
out_direc = configs$outdir
se_fname = file.path(
    out_direc, paste( configs$outprefix, "summarized_experiment.RData", sep="_")
)
if (!dir.exists(out_direc)) {
    dir.create(out_direc)
}
 
# set up our features for analysis with DeSeq2
conditions = read_csv(configs$deseq$sample_table, col_names = TRUE, show_col_types=FALSE) %>% 
    mutate(sampletype = fct_relevel(sampletype, configs$deseq$reference_sampletype))

# set desired sample type as the baseline
bamfiles = BamFileList(conditions$bamfile, asMates=TRUE)

if (configs$input$annotation_type == "gff") {
    features = gffRead(configs$input$target_locs)
}
if (configs$input$annotation_type == "narrowpeak") {
    features = narrowpeakRead(configs$input$target_locs)
}
if (configs$input$annotation_type == "bed") {
    features = bedRead(configs$input$target_locs, cols=configs$input$cols)
}

if (configs$input$feature_filter == "") {
    CDSs = features
} else {
    CDSs = features %>% filter(feature == configs$input$feature_filter)
}

#CDSs$gene = replace_na(
#    getAttributeField(CDSs$attributes, configs$input$locus_name_field),
#    "unknown"
#)
CDSs$locus_tag = CDSs$name

CDS_gr = as(CDSs, "GRanges")

print(paste0("Summarizing overlaps between reads and annotations..."))
se = summarizeOverlaps(
    CDS_gr,
    bamfiles,
    singleEnd = TRUE,
    inter.feature = FALSE,
    ignore.strand = TRUE,
    mode = "IntersectionStrict",
)

colData(se) = DataFrame(conditions)

save(conditions, se, file=se_fname)
