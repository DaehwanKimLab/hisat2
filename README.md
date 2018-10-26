# Graph-Based Genome Alignment and Genotyping with HISAT2 and HISAT-genotype

For more information, see the following websites:
* [HISAT2 website](http://ccb.jhu.edu/software/hisat2)
* [HISAT-genotype website](http://ccb.jhu.edu/software/hisat-genotype)

## Contact

[Daehwan Kim](https://kim-lab.org) (infphilo@gmail.com) and [Chanhee Park](https://www.linkedin.com/in/chanhee-park-97677297/) (parkchanhee@gmail.com)

## Abstract

Rapid advances in next-generation sequencing technologies have dramatically changed our ability to perform genome-scale analyses. The human reference genome used for most genomic analyses represents only a small number of individuals, limiting its usefulness for genotyping. We designed a novel method, HISAT2, for representing and searching an expanded model of the human reference genome, in which a large catalogue of known genomic variants and haplotypes is incorporated into the data structure used for searching and alignment. This strategy for representing a population of genomes, along with a fast and memory-efficient search algorithm, enables more detailed and accurate variant analyses than previous methods. We demonstrate two initial applications of HISAT2: HLA typing, a critical need in human organ transplantation, and DNA fingerprinting, widely used in forensics. These applications are part of HISAT-genotype, with performance not only surpassing earlier computational methods, but matching or exceeding the accuracy of laboratory-based assays.

![](HISAT2-genotype.png)

## HISAT2
HISAT2 is a fast and sensitive alignment program for mapping next-generation sequencing reads (whole-genome, transcriptome, and exome sequencing data) to a population of human genomes (as well as to a single reference genome). Based on an extension of BWT for a graph [1], we designed and implemented a graph FM index (GFM), an original approach and its first implementation to the best of our knowledge. In addition to using one global GFM index that represents general population, HISAT2 uses a large set of small GFM indexes that collectively cover the whole genome (each index representing a genomic region of 56 Kbp, with 55,000 indexes needed to cover human population). These small indexes (called local indexes) combined with several alignment strategies enable effective alignment of sequencing reads. This new indexing scheme is called Hierarchical Graph FM index (HGFM). We have developed HISAT2 based on the HISAT [2] and Bowtie 2 [3] implementations.  See the [HISAT2 website](http://ccb.jhu.edu/software/hisat2/index.shtml) for
more information.

A few notes:

1) HISAT2's index (HGFM) size for the human reference genome and 12.3 million common SNPs is 6.2GB. The SNPs consist of 11 million single nucleotide polymorphisms, 728,000 deletions, and 555,000 insertions. Insertions and deletions used in this index are small (usually <20bp). We plan to incorporate structural variations (SV) into this index.

2) The memory footprint of HISAT2 is relatively low, 6.7GB.

3) The runtime of HISAT2 is estimated to be slightly slower than HISAT (30–100% slower for some data sets).

4) HISAT2 provides greater accuracy for alignment of reads containing SNPs.

5) We released a first (beta) version of HISAT2 in September 8, 2015.

## License

[GPL-3.0](LICENSE)

# For reviwers:

## Code

The Code directory (`/code`) contains a specific version of [HISAT2 and HISAT-genotype](http://github.com/infphilo/hisat2) at GitHub (a branch named hisat2_v2.2.0_beta).

## Data

The Data directory (`/data`) contains all input files for reproducing some of our results such as from the 
evaluation of HISAT2 and other programs using both simulated and real reads, from typing and assembling
HLA genes of Illumina Platinum Genomes using HISAT-genotype, and from building a HISAT2 graph index.

* **Simulated read pairs**

| Type | Number of pairs | Path |
| - | - | - |
| SNPs and sequencing errors included | 10,000,000 | /data/reads/simulation/10M_DNA_mismatch_snp_reads_genome |
| SNPs included | 10,000,000 | /data/reads/simulation/10M_DNA_snp_reads_genome |
| Sequencing errors included | 10,000,000 | /data/reads/simulation/10M_DNA_mismatch_reads_genome |
| No SNPs nor sequencing errors included | 10,000,000 | /data/reads/simulation/10M_DNA_reads_genome |

Each directory comes with a true alignment file in SAM format so that users know where the reads were generated in the human reference genome.

* **Real read pairs**

| Number of read pairs | Path |
| - | - |
| 10,000,000 | /data/reads/real/10M |

* **Human reference genome, SNPs, haplotypes, and HISAT2's indexes**

| Type | Path |
| - | - |
| GRCh38 reference | /data/data/genome.fa |
| SNPs | /data/data/genome.snp |
| Haplotypes | /data/data/genome.haplotype |
| HISAT2's prebuilt graph index for comparison with other aligners | /data/indexes/HISAT2/genome.[1-8].ht2 |
| HISAT2's prebuilt graph index for genotyping | /data/indexes/HISAT2/genotype_genome.[1-8].ht2 |

## Running the Code Ocean pipeline

The [run.sh](run.sh) script includes the steps required to reproduce some of the results. 
It creates symbolic links to all data files in appropriate directories and executes several command lines 
depending on user-inputted options as detailed in the table below.

| Test case | A comma-separated list of aligners | Reference genome | Comment |
| - | - | - | - |
| 1 | default: hisat2 or use hisat2,bwa,bowtie2,vg | default: human genome or use chr22 | Evaluation using simulated read pairs with SNPs and sequencing errors |
| 2 | default: hisat2 or use hisat2,bwa,bowtie2,vg | default: human genome or use chr22 | Evaluation using simulated read pairs with SNPs and without sequencing errors |
| 3 | default: hisat2 or use hisat2,bwa,bowtie2,vg | default: human genome or use chr22 | Evaluation using simulated read pairs without no SNPs and with sequencing errors |
| 4 | default: hisat2 or use hisat2,bwa,bowtie2,vg | default: human genome or use chr22 | Evaluation using simulated read pairs without SNPs and without sequencing errors |
| 5 | default: hisat2 or use hisat2,bwa,bowtie2,vg | N/A | Evaluation using 10 million real read pairs |
| 6 | N/A | N/A | Typing and assembling HLA genes of Illumina Platinum Genomes |
| 7 | N/A | default: human genome or use chr22 | Building a HISAT2 graph index (a hierarchical graph FM index) |

Test cases 1-5 are expected to take several hours to several days depending on whether you use the whole genome or a single chromosome and which aligners you use. Test case 6 can be completed within one hour. Test case 7 takes about 3 hours and requires a compute node with at least 256 GB of RAM if you want to build a graph index for the human reference genome with 14 million common SNPs.

## Results

The Results directory (`/results`) is organized as follows:

## References

[1] Sirén J, Välimäki N, Mäkinen V (2014) Indexing graphs for path queries with applications in genome research. IEEE/ACM Transactions on Computational Biology and Bioinformatics 11: 375–388. doi: 10.1109/tcbb.2013.2297101

[2] Kim D, Langmead B, and Salzberg SL  HISAT: a fast spliced aligner with low memory requirements, Nature methods, 2015

[3] Langmead B, Salzberg SL: Fast gapped-read alignment with Bowtie 2. Nat Methods 2012, 9:357-359
