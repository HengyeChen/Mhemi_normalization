# Normalize Mhemi-seq data
We present a method named Mhemi-seq (MspJI-assisted hemi-methylation sequencing) to detect CpG methylation and hemi-methylation in genome. The script here is used to correct the cutting bias of MspJI in Mhemi-seq. 

# Input file
## Raw Mhemi-seq results of GM12878
GM12878_Mhemi_rep123_CpG_dyad.5000.bed
## The coordination of YNCGNR motifs in human hg38 reference genome
hg38_YNCGNR.fa.1000.bed


# Output file
## Normalized Mhemi-seq results
GM12878_Mhemi_rep123_CpG_dyad.YNCGNR.v4.3.bed
