---
layout: page
title: Main
permalink: /main/
order: 1
share: false
---

**HISAT2** is a fast and sensitive alignment program for mapping next-generation sequencing reads (both DNA and RNA) to a population of human genomes as well as to a single reference genome. Based on an extension of BWT for graphs ([Sir&eacute;n et al. 2014](http://dl.acm.org/citation.cfm?id=2674828)), we designed and implemented a graph FM index (GFM), an original approach and its first implementation. In addition to using one global GFM index that represents a population of human genomes, **HISAT2** uses a large set of small GFM indexes that collectively cover the whole genome. These small indexes (called local indexes), combined with several alignment strategies, enable rapid and accurate alignment of sequencing reads. This new indexing scheme is called a Hierarchical Graph FM index (HGFM).

### [The HISAT-3N paper](https://genome.cshlp.org/content/31/7/1290.abstract) published at *Genome Research*. 7/1/2021

### HISAT-3N beta release 12/14/2020

HISAT-3N is a software system for analyzing nucleotide conversion sequencing reads. See the [HISAT-3N] for more details.

[HISAT-3N]:	{{ site.baseurl }}{% link _pages/hisat-3n.md %}

### Index files are moved to the AWS Public Dataset Program. 9/3/2020

We have moved HISAT2 index files to the AWS Public Dataset Program. See the [link](https://registry.opendata.aws/jhu-indexes/) for more details.


### HISAT 2.2.1 release 7/24/2020

This patch version includes the following changes.
* Python3 support
* Remove the HISAT-genotype related scripts. HISAT-genotype moved to [http://daehwankimlab.github.io/hisat-genotype/](http://daehwankimlab.github.io/hisat-genotype/)
* Fixed bugs related to `--read-lengths` option


### HISAT 2.2.0 release 2/6/2020

This major version update includes a new feature to handle “repeat” reads. Based on sets of 100-bp simulated and 101-bp real reads that we tested, we found that 2.6-3.4% and 1.4-1.8% of the reads were mapped to >5 locations and >100 locations, respectively. Attempting to report all alignments would likely consume a prohibitive amount of disk space. In order to address this issue, our repeat indexing and alignment approach directly aligns reads to repeat sequences, resulting in one repeat alignment per read. HISAT2 provides application programming interfaces (API) for C++, Python, and JAVA that rapidly retrieve genomic locations from repeat alignments for use in downstream analyses.  
Other minor bug fixes are also included as follows:  

* Fixed occasional sign (+ or -) issues of template lengths in SAM file
* Fixed duplicate read alignments in SAM file
* Skip a splice site if exon's last base or first base is ambiguous (N) 


### Index files are moved to a different location. 8/30/2019

Due to a high volume of index downloads, we have moved HISAT2 index files to a different location in order to provide faster download speed. If you use wget or curl to download index files, then you may need to use the following commands to get the correct file name.
* `wget --content-disposition` *download_link*
* `curl -OJ` *download_link*


### [The HISAT2 paper](https://www.nature.com/articles/s41587-019-0201-4) is out in *Nature Biotechnology*. 8/2/2019


### HISAT 2.1.0 release 6/8/2017

* This major version includes the first release of HISAT-genotype, which currently performs HLA typing,
  DNA fingerprinting analysis, and CYP typing on whole genome sequencing (WGS) reads. 
  We plan to extend the system so that it can analyze not just a few genes, but a whole human genome. 
  Please refer to [the HISAT-genotype website](https://daehwankimlab.github.io/hisat-genotype) for more details.
* HISAT2 can be directly compiled and executed on Windows system using Visual Studio, thanks to [Nigel Dyer](http://www2.warwick.ac.uk/fac/sci/systemsbiology/staff/dyer/).
* Implemented `--new-summary` option to output a new style of alignment summary, which is easier to parse for programming purposes.
* Implemented `--summary-file` option to output alignment summary to a file in addition to the terminal (e.g. stderr).
* Fixed discrepancy in HISAT2’s alignment summary.
* Implemented `--no-templatelen-adjustment` option to disable automatic template length adjustment for RNA-seq reads.


### HISAT2 2.0.5 release 11/4/2016
Version 2.0.5 is a minor release with the following changes.
* Due to a policy change (HTTP to HTTPS) in using SRA data (`--sra-option`), users are strongly encouraged to use this version. As of 11/9/2016, NCBI will begin a permanent redirect to HTTPS, which means the previous versions of HISAT2 no longer works with `--sra-acc` option soon.
* Implemented `-I` and `-X` options for specifying minimum and maximum fragment lengths.  The options are valid only when used with `--no-spliced-alignment`, which is used for the alignment of DNA-seq reads.
* Fixed some cases where reads with SNPs on their 5' ends were not properly aligned.
* Implemented `--no-softclip` option to disable soft-clipping.
* Implemented `--max-seeds` to specify the maximum number of seeds that HISAT2 will try to extend to full-length alignments (see [the manual] for details).


### [HISAT, StringTie and Ballgown protocol](http://www.nature.com/nprot/journal/v11/n9/full/nprot.2016.095.html) published at Nature Protocols 8/11/2016

  
### HISAT2 2.0.4 Windows binary available [here](http://www.di.fc.ul.pt/~afalcao/hisat2_windows.html), thanks to [Andre Osorio Falcao](http://www.di.fc.ul.pt/~afalcao/) 5/24/2016

	  
### HISAT2 2.0.4 release 5/18/2016
Version 2.0.4 is a minor release with the following changes.
* Improved template length estimation (the 9th column of the SAM format) of RNA-seq reads by taking introns into account.
* Introduced two options, `--remove-chrname` and `--add-chrname`, to remove "chr" from reference names or add "chr" to reference names in the alignment output, respectively (the 3rd column of the SAM format).
* Changed the maximum of mapping quality (the 5th column of the SAM format) from 255 to 60. Note that 255 is an undefined value according to the SAM manual and some programs would not work with this value (255) properly.
* Fixed NH (number of hits) in the alignment output.
* HISAT2 allows indels of any length pertaining to minimum alignment score (previously, the maximum length of indels was 3 bp).
* Fixed several cases that alignment goes beyond reference sequences.
* Fixed reporting duplicate alignments.


### HISAT2 2.0.3-beta release 3/28/2016
Version 2.0.3-beta is a minor release with the following changes.
* Fixed graph index building when using both SNPs and transcripts. As a result, genome_snp_tran indexes here on the HISAT2 website have been rebuilt.
* Included some missing files needed to follow the small test example (see [the manual] for details).


### HISAT2 2.0.2-beta release 3/17/2016
**Note (3/19/2016):** this version is slightly updated to handle reporting splice sites with the correct chromosome names.
Version 2.0.2-beta is a major release with the following changes.
* Memory mappaped IO (`--mm` option) works now.
* Building linear index can be now done using multi-threads.
* Changed the minimum score for alignment in keeping with read lengths, so it's now `--score-min L,0.0,-0.2`, meaning a minimum score of -20 for 100-bp reads and -30 for 150-bp reads.
* Fixed a bug that the same read was written into a file multiple times when `--un-conc` was used.
* Fixed another bug that caused reads to map beyond reference sequences.
* Introduced `--haplotype` option in the hisat2-build (index building), which is used with `--snp` option together to incorporate those SNP combinations present in the human population.  This option also prevents graph construction from exploding due to exponential combinations of SNPs in small genomic regions.
* Provided a new python script to extract SNPs and haplotypes from VCF files, <i>hisat2_extract_snps_haplotypes_VCF.py</i>
* Changed several python script names as follows<
  * *extract_splice_sites.py* to *hisat2_extract_splice_sites.py*
  * *extract_exons.py* to *hisat2_extract_exons.py*
  * *extract_snps.py* to *hisat2_extract_snps_haplotypes_UCSC.py*


### HISAT2 2.0.1-beta release 11/19/2015
Version 2.0.1-beta is a maintenance release with the following changes.
* Fixed a bug that caused reads to map beyond reference sequences.
* Fixed a deadlock issue that happened very rarely.
* Fixed a bug that led to illegal memory access when reading SNP information.
* Fixed a system-specific bug related to popcount instruction.


### HISAT2 2.0.0-beta release 9/8/2015 - first release
We extended the BWT/FM index to incorporate genomic differences among individuals into the reference genome, while keeping memory requirements low enough to fit the entire index onto a desktop computer. Using this novel Hierarchical Graph FM index (HGFM) approach, we built a new alignment system, HISAT2, with an index that incorporates ~12.3M common SNPs from the dbSNP database. HISAT2 provides greater alignment accuracy for reads containing SNPs.
* HISAT2's index size for the human reference genome and 12.3 million common SNPs is 6.2GB (the memory footprint of HISAT2 is 6.7GB). The SNPs consist of 11 million single nucleotide polymorphisms, 728,000 deletions, and 555,000 insertions. The insertions and deletions used in this index are small (usually <20bp).
* HISAT2 comes with several index types:
  * Hierarchical FM index (HFM) for a reference genome (index base: <i>genome</i>)
  * Hierarchical Graph FM index (HGFM) for a reference genome plus SNPs (index base: <i>genome_snp</i>)
  * Hierarchical Graph FM index (HGFM) for a reference genome plus transcripts (index base: <i>genome_tran</i>)
  * Hierarchical Graph FM index (HGFM) for a reference genome plus SNPs and transcripts (index base: <i>genome_snp_tran</i>)
* HISAT2 is a successor to both [HISAT](http://ccb.jhu.edu/software/hisat) and [TopHat2](http://ccb.jhu.edu/software/tophat). We recommend that HISAT and TopHat2 users switch to HISAT2.
  * HISAT2 can be considered an enhanced version of HISAT with many improvements and bug fixes. The alignment speed and memory requirements of HISAT2 are virtually the same as those of HISAT when using the HFM index (<i>genome</i>).
  * When using graph-based indexes (HGFM), the runtime of HISAT2 is slightly slower than HISAT (30~80% additional CPU time).
  * HISAT2 allows for mapping reads directly against transcripts, similar to that of TopHat2 (use <i>genome_tran</i> or <i>genome_snp_tran</i>).
* When reads contain SNPs, the SNP information is provided as an optional field in the SAM output of HISAT2 (e.g., **<code>Zs:Z:1|S|rs3747203,97|S|rs16990981</code>** - see [the manual] for details).  This feature enables fast and sensitive genotyping in downstream analyses. Note that there is no alignment penalty for mismatches, insertions, and deletions if they correspond to known SNPs.
* HISAT2 provides options for transcript assemblers (e.g., StringTie and Cufflinks) to work better with the alignment from HISAT2 (see options such as `--dta` and `--dta-cufflinks`).
* Some slides about HISAT2 are found [here]({{ '/assets/data/HISAT2-first_release-Sept_8_2015.pdf' | prepend: site.baseurl }}) and we are preparing detailed documention.
* We plan to incorporate a larger set of SNPs and structural variations (SV) into this index (e.g., long insertions/deletions, inversions, and translocations).

[the manual]: {{ site.baseurl }}{% link _pages/manual.md %}

### The HISAT2 source code is available in a [public GitHub repository](https://github.com/DaehwanKimLab/hisat2) (5/30/2015).


