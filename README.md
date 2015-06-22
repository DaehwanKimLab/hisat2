# hisat2 
HISAT 2 is a fast and sensitive alignment program for mapping next-generation sequencing reads (whole-genome, transcriptome and exome sequencing data) against general population (as well as a single referece genome). Based on generalized BWT [1], we designed and implemented generalized FM index (GFM), which is the first implementation of GFM to the best of our knowledge, and in addition to one global GFM index that represents general population, HISAT 2 uses a large set of small GFM indexes that collectively cover the whole genome (each index represents a genomic region of ~57,000 bp and ~54,000 indexes are needed to cover human population). These small indexes (called local indexes) combined with several alignment strategies enable effective alignment of sequencing reads. This new indexing scheme is called Hierarchical Generalized FM index (HGFM).  We have developed HISAT 2 based on the HISAT and Bowtie 2 implementations.

A few notes: 

1) HISAT 2's index (HGFM) size for the human reference genome and ~12.3 million common SNPs is ~6.3GB.  The SNPs consist of ~11 million single nucelotide polymorphisms, ~728,000 deletions and ~555,000 insertions.  Insertions and deletions used in this index are small (usaully <20bp).  We plan to incorporate structural variations (SV) into this index.

2) The memory footprint of HISAT 2 is relatively low, ~6.9GB.

3) The runtime of HISAT 2 is estimated to be slightly slower than HISAT (30 ~ 100% slower for some data sets).

4) HISAT 2 provides greater accuracy for alignment of reads containing SNPs.


We plan to release an alpha version of HISAT 2 in August this year.


[1] Sirén J, Välimäki N, Mäkinen V (2014) Indexing graphs for path queries with applications in genome research. IEEE/ACM Transactions on Computational Biology and Bioinformatics 11: 375–388. doi: 10.1109/tcbb.2013.2297101 
