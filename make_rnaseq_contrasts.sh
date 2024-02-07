#!/usr/bin/bash

set -e

for geno in wt vpsR vpsT flrA; do
    for time in t0 t1 t2 t3; do

        active_qrgB="genotype${geno}.time${time}.sampletypep2"

        Rscript /corexfs/schroedj/src/rnap_chip_analysis/src/get_deseq_results_for_contrast.R \
            rnaseq_analysis.cfg \
            ${active_qrgB} \
            "" \
            ${geno}_${time}_p2_vs_p1_results.csv

    done
done

