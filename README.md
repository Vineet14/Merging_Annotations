# Merging_Annotations

### An effort to merge Ensembl and UCSC "hg19" annotations for use in splicing analysis and filter Miso v2.0 annotations.

Extract the files in the "input_files.tar.xz". On extraction the "input_files" directory should contain:
  all_miso.gtf,
  overlapping_genes.gtf,
  removed_overlapped_ensembl.gtf,
  ucsc_refseq.gtf

hg19 versions hg37.87 annotations "removed_overlapped_ensembl.gtf" from Ensembl were used with cotemporary UCSC refseq annotations "ucsc_refseq.gtf" to be merged.

Ensembl annotations were used as a base for the merging of annotations.

All ucsc transcripts with overlap with Ensembl annotations were removed while preparing the merge.

Genes which overlapped with each other "overlapping_genes.gtf" were removed before the merge 
    and subsequently added to the final merge. 

Miso annotations "all_miso.gtf" from MISO version v2.0 were used and were mapped back to the merged gtf files. 
    This was done to only identify miso annotations for splicing analysis that were represented in modern upto-date hg19 annotations.

Run command

### python3 merging_UCSC_Ensembl_gtfs.py <path_to-ucsc_refseq.gtf> <path_to-all_miso.gtf> <path_to-removed_overlapped_ensembl.gtf> <output_merged_file>

### cat <output_merged_file> <path_to-overlapping_genes.gtf> > merged.gtf

# For generation of the miso file with filtered miso annotations present in the merged file

This script is to identify miso annotations present in the merge of the upto0date standard annotation from ensembl and ucsc repositories

### python3 Modifying_Miso_annotation.py <output_merged_file> <path_to-all_miso.gtf> <final miso filtered file>
