Below is the code run to produce the log2 fold-changes in the first submission of Onder. \<data\_location\> was replaced with the actual directory containing the alignment files.

```bash
cd <data_location>

samtools view -SH dksa1_tn_alignment_synlethal.sam > dksa1_insertions.sam
samtools view -S -F 4 -f 64 dksa1_tn_alignment_synlethal.sam | awk -v OFS="\t" -v FS="\t" '{print $1,$2,$3,$4,$4+1,$6,$7,$8,$9,$10,$11}' >> dksa1_insertions.sam
samtools view -bS dksa1_insertions.sam > dska1_insertions.bam

samtools view -SH dksa2_tn_alignment_synlethal.sam > dksa2_insertions.sam
samtools view -S -F 4 -f 64 dksa2_tn_alignment_synlethal.sam | awk -v OFS="\t" -v FS="\t" '{print $1,$2,$3,$4,$4+1,$6,$7,$8,$9,$10,$11}' >> dksa2_insertions.sam
samtools view -bS dksa2_insertions.sam > dska2_insertions.bam

samtools view -SH wt1_tn_alignment_synlethal.sam > wt1_insertions.sam
samtools view -S -F 4 -f 64 wt1_tn_alignment_synlethal.sam | awk -v OFS="\t" -v FS="\t" '{print $1,$2,$3,$4,$4+1,$6,$7,$8,$9,$10,$11}' >> wt1_insertions.sam
samtools view -bS wt1_insertions.sam > wt1_insertions.bam

samtools view -SH wt2_tn_alignment_synlethal.sam > wt2_insertions.sam
samtools view -S -F 4 -f 64 wt2_tn_alignment_synlethal.sam | awk -v OFS="\t" -v FS="\t" '{print $1,$2,$3,$4,$4+1,$6,$7,$8,$9,$10,$11}' >> wt2_insertions.sam
samtools view -bS wt2_insertions.sam > wt2_insertions.bam

cd <analysis_direc>

Rscript summarize_overlaps.R rnaseq_analysis.cfg
Rscript fit_deseq_model.R rnaseq_analysis.cfg

Rscript get_deseq_possible_contrasts.R rnaseq_analysis.cfg

Rscript get_deseq_results_for_contrast.R rnaseq_analysis.cfg sampletypewt sampletypedksA wt_vs_dksA_results.csv
```
