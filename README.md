# mRNA and DNA reads

Simulate both __DNA__ and __RNA-seq__ data with specified transcript abundances in order to benchmark allele-specific expression algorithms.

#### Inputs:

1. Reference file (fasta)
2. GFF (Gene and transcript coordinates)
3. dbSNP database.

#### Outputs:

1. Paired-end DNA fastq files
2. Paired-end mRNA fastq files
3. TSV file with SNPs


### Description

From the GFF file we take exon coordinates and intersect them with the dbSNP to select SNPs that belong to exons. We also add a number of random SNVs per exon and we randomly assign a phase to each mutation. Then, we create two contigs (haplotypes) from a reference contig by introducing mutations based on the phase value.
