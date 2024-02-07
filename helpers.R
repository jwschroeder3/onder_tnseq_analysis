# this file contains numerous helper functions useful for analyzing RNA pol ChIP data
# return TRUE iff a given file exists
check_file_exists = function(in_fname) {
    return(file.exists(in_fname))
}

read_bedgraph = function(in_fname, score_colname="score") {
    print(in_fname)
    require(tidyverse)
    data = read_table(
        in_fname,
        col_names=c("seqname", "start", "end", score_colname),
        col_types = "ciin"
    )
    
    return(data)
}

# read coverage information from a file, by default in bedgraph format
read_coverage = function(in_fname, .col_names=c("chr", "start", "end", "score")) {
    if (!file.exists(in_fname)) {
        inp_fname = str_replace(in_fname, "input", "inp")
        if (!file.exists(inp_fname)) {
            print(paste0("Could not find either ", in_fname, " or ", inp_fname, "!"))
        }
        in_fname = inp_fname
    }
    print(paste0("Reading ", in_fname))
    coverage = read_table(in_fname, col_names=.col_names) %>%
        mutate(bin = factor(1:n()))
    return(coverage)
}

# return TRUE iff any value in $x is NA
any_is_na = function(x) {
    return(any(is.na(x)))
}

# return TRUE iff every value in $x is NA
all_is_na = function(x) {
    return(all(is.na(x)))
}

# return TRUE iff the fraction of zeros in $x is greater than $frac
frac_zero = function(x, frac=0.5) {
    total = length(x)
    frac_zero = length(which(x < 0)) / total
    return (frac_zero >= frac)
}

# read a gff file and return it as a data frame
gffRead = function(gffFile, nrows = -1) {
    cat("Reading ", gffFile, sep="")
    gff = read.table(
        gffFile,
        sep="\t",
        as.is=TRUE,
        quote="",
        header=FALSE,
        comment.char="#",
        nrows = nrows,
        colClasses=c("character", "character", "character", "integer","integer",
             "character", "character", "character", "character")
    )
    colnames(gff) = c("seqname", "source", "feature", "start", "end",
             "score", "strand", "frame", "attributes")
    cat("found", nrow(gff), "rows with classes:",
         paste(sapply(gff, class), collapse=", "), "\n")
    return(gff)
}

bedRead = function(bedFile, nrows=-1, cols=NA) {
    cat("reading ", bedFile, sep="")
    bed = read.table(
        bedFile,
        sep="\t",
        as.is=TRUE,
        quote="",
        header=FALSE,
        comment.char="#",
        nrows=nrows
    )
    if (is.na(cols[1])) {
        cols = c("seqname", "start", "end")
        if (dim(bed)[2] > 3) {
            cols = c(cols, "name", "score", "strand")
        }
        if (dim(bed)[2] > 6) {
            cols = c(cols, "thickStart", "thickEnd", "itemRgb", "blockCount", "blockSizes", "blockStarts")
        }
    }
    colnames(bed) = cols
    return(bed)
}

narrowpeakRead = function(narrowpeakFile, nrows=-1) {
    cat("Reading ", narrowpeakFile, "\n", sep="")
    np = read.table(
        narrowpeakFile,
        sep="\t",
        as.is=TRUE,
        quote="",
        header=FALSE,
        comment.char="#",
        nrows = nrows,
        colClasses=c("character", "integer", "integer", "character","numeric",
             "character", "numeric", "numeric", "numeric", "numeric")
    )
    colnames(np) = c("seqname", "start", "end", "name", "display",
             "strand", "score", "pval", "qval", "peak")
    if (all(np$name == ".")) {
        np$name = paste0("peak_", 1:nrow(np))
    }
    cat("found", nrow(np), "rows with classes:",
         paste(sapply(np, class), collapse=", "), "\n")
    return(np)
}

# pull an attribute field named $field from a data frame that was read using gffRead
# attribute fields occur in the ninth column of a gff file, and are separated
# using $attrsep
getAttributeField = function(x, field, attrsep = ";") {
    s = strsplit(x, split = attrsep, fixed = TRUE)
    sapply(s, function(atts) {
        a = strsplit(atts, split = "=", fixed = TRUE)
        m = match(field, sapply(a, "[", 1))
        if (!is.na(m)) {
            rv = a[[m]][2]
            return(rv)
        }

        a = strsplit(atts, split = " ", fixed = TRUE)
        m = match(field, sapply(a, "[", 1))
        if (!is.na(m)) {
            rv = a[[m]][2]
            return(rv)
        }

        else {
            rv = as.character("unknown")
        }
        return(as.character("unknown"))
    })
}

# check whether either $in_fname exists or if it exists with 'input' replaced by 'inp'
# return the correct filename if it exists, otherwise raise an error
set_fname = function(in_fname) {
    if (!file.exists(in_fname)) {
        inp_fname = str_replace(in_fname, "input", "inp")
        if (!file.exists(inp_fname)) {
            stop(paste0("Could not find either ", in_fname, " or ", inp_fname, "!"))
        }
        in_fname = inp_fname
    }
    return(in_fname)
}

# obtain the results from a deseq fitted object $deseq_data for the specified $contrast
get_contrast_results = function(contrast, deseq_data) {
    res = results(
        deseq_data,
        contrast = contrast
    )
    rowdat = rowData(deseq_data)
    res = as_tibble(res) %>%
        bind_cols(as_tibble(rowdat) %>% dplyr::select(peak_id))
    
    return(res)
}

# get a table of results for a given contrast and fitted deseq object
get_results = function(.contrast, .deseq_dat) {
    res = results(
        .deseq_dat,
        contrast = .contrast
    )
    res = as_tibble(res) %>%
        bind_cols(as_tibble(rowData(dds)) %>% dplyr::select(peak_id))
    return(as_tibble(res))
}

# plot a given locus' data
plotLocus = function(feats_df, signal_df, plotStart, plotEnd, chr_name,
                     lineAlpha=1, lineSize=1.25,
                     yvar="log2fc",
                     color_var="strain",
                     ylabel="",
                     feat_var = "locus_tag",
                     facet = NULL,
                     ylims="detect", plotFeatures=TRUE, arrowLength=1.5,
                     plotMotifLocs=FALSE, motifs_df=NULL,
                     strand_colors=c("#4667FC", "#FB4B13"),
                     ordinal_color=FALSE
) {

    #print("here")
    feat_var = sym(feat_var)
    yvar = sym(yvar)
    color_var = sym(color_var)
  
    plot_feats = feats_df %>% 
        dplyr::filter(seqname==chr_name, end > plotStart, start < plotEnd) %>%
        mutate(midpoint = as.numeric((start + end)/2))
    #print("signal_df")
    #print(signal_df %>% head)
    plot_sig = signal_df %>%
        dplyr::filter(seqname==chr_name, start<plotEnd, end>plotStart) %>%
        mutate(position = (start + end) / 2)
    #print("plot_sig")
    #print(plot_sit %>% head)
    #stop()

    #print("there")
    if (ordinal_color) {
        colors = sort(unique(unlist(plot_sig[,color_var], use.names=FALSE)))
        plot_sig = plot_sig %>%
            mutate(!!color_var := factor(!!color_var, levels=sort(unique(!!color_var))))
    }
    if (plotMotifLocs) {
        plotMotifs = motifs_df %>% 
            dplyr::filter(seqname==chr_name, end > plotStart, start < plotEnd)
    }
  
    #print("there there")
    #print(plot_sig[,yvar])
    min_chip = min(plot_sig[,yvar])
    max_chip = max(plot_sig[,yvar])
    chip_range = max_chip - min_chip
    geneTop = min_chip - 0.05*chip_range
    geneBottom = geneTop - 0.075*chip_range
    geneText = geneBottom - 0.05 * chip_range
    
    #print("everywhere")
    plot = ggplot()

    if (plotMotifLocs) {
        #print(plotMotifs %>% head)
        plot = plot +
            geom_rect(
                data=plotMotifs,
                aes(
                    xmin=start/1e6,
                    xmax=end/1e6,
                    ymin=geneTop,
                    ymax=Inf
                ),
                alpha=0.4
            )
    }
   
    #print("all at once")
    if (plotFeatures) {
        #print(plot_feats %>% head)
        plot = plot + 
            geom_rect(
                data=plot_feats,
                aes(
                    xmin=start/1e6,
                    xmax=end/1e6,
                    ymin=geneBottom,
                    ymax=geneTop,
                    fill=strand
                ), 
                color="black"
            ) +
            scale_fill_manual(values=strand_colors) +
            geom_text(
                data=plot_feats,
                aes(
                    x=midpoint/1e6,
                    y=geneText,
                    label=!!feat_var
                ),
                parse=T
            )
    }

    #print("down here now")
    #print(plot_sig %>% head)
    plot = plot + 
        geom_line(
            data=plot_sig,
            aes(
                x=position/1e6,
                y=!!yvar,
                color=!!color_var
            ),
            size=lineSize,
            alpha=lineAlpha
        ) +
        theme_classic() +
        theme(
            text = element_text(size=18),
            axis.line = element_line(size=2),
            axis.text = element_text(size=14, color="black"),
            axis.ticks = element_line(color="black")
        ) +
        labs(y=ylabel, x="Genome position (Mb)") +
        coord_cartesian(xlim = c(plotStart/1e6, plotEnd/1e6)) +
        guides(colour = guide_legend(override.aes = list(alpha = 1)))

    #print("keep going")
    if (ordinal_color) {
        cols = grDevices::colorRampPalette(c("#2b2c2e","#ed8002","#ed1400"))(length(colors))
        plot = plot + scale_color_manual(values = cols)
    }
  
    if (!(ylims == "detect")) {
        plot = plot + coord_cartesian(
            xlim = c(plotStart/1e6, plotEnd/1e6),
            ylim=ylims
        )
    }

    if (!is.null(facet)) {
        plot = plot + facet_grid(facet)
    }

    return(plot)
}

