---
layout: page
title: HISAT-3N 
permalink: /hisat-3n/
order: 4
share: false
---

HISAT-3N
============

Overview
-----------------
**HISAT-3N** (hierarchical indexing for spliced alignment of transcripts - 3 nucleotides)
is designed for nucleotide conversion sequencing technologies and implemented based on [HISAT2]. 
There are two strategies for HISAT-3N to align nucleotide conversion sequencing reads: *standard mode* and *repeat mode*. 
The standard mode aligns reads with a standard 3N index only, so it is fast and requires only a small amount of memory (~9GB for human genome alignment).
The repeat mode aligns reads with both a standard 3N index and a repeat 3N index, then outputs up to 1,000 alignment results (the number of outputted alignments can be adjusted by `--repeat-limit`).
The repeat mode can also align nucleotide conversion reads more accurately, 
and it is only 10% slower than the standard mode with slightly more memory requirements (the repeat mode uses about ~10.5GB).

HISAT-3N can be used for any nucleotide-converted sequencing reads, including [BS-seq], [SLAM-seq], [TAB-seq], [oxBS-seq], [TAPS], [scBS-seq], and [scSLAM-seq].

[HISAT2]:https://github.com/DaehwanKimLab/hisat2
[BS-seq]: https://en.wikipedia.org/wiki/Bisulfite_sequencing
[SLAM-seq]: https://www.nature.com/articles/nmeth.4435
[scBS-seq]: https://www.nature.com/articles/nmeth.3035
[scSLAM-seq]: https://www.nature.com/articles/s41586-019-1369-y
[TAPS]: https://www.nature.com/articles/s41587-019-0041-2
[TAB-seq]: https://doi.org/10.1016/j.cell.2012.04.027
[oxBS-seq]: https://science.sciencemag.org/content/336/6083/934


Getting started
============
The HISAT-3N alignment process requires a 64-bit computer running either Linux or Mac OS and at least 16GB of RAM. 

A few notes:  

1. Building the standard 3N index requires 16GB of RAM or less.
2. Building the repeat 3N index requires 256GB of RAM.
3. The alignment process using either the standard or repeat index requires less than 16GB of RAM.
4. [SAMtools] is required to sort SAM files in order to generate a HISAT-3N table.

[SAMtools]:        http://samtools.sourceforge.net

Install
------------
   
    git clone https://github.com/DaehwanKimLab/hisat2.git
    cd hisat2
    git checkout -b hisat-3n origin/hisat-3n
    make


Make sure that you select the `hisat-3n` branch


Build a 3N index with `hisat-3n-build`
-----------
`hisat-3n-build` builds a 3N-index, which internally contains two [HISAT2] indexes for a set of DNA sequences. For the standard 3N-index, 
each index contains 16 files with suffix `.3n.*.*.ht2`.
For the repeat 3N-index, there are 16 more files in addition to the standard 3N-index, and these files have the suffix 
`.3n.*.rep.*.ht2`. 
These files constitute the entirety of the HISAT-3N index.

* An example for building a standard HISAT-3N index:  
`hisat-3n-build genome.fa genome`  

* An example for building a repeat HISAT-3N index, which requires 256GB memory:  
`hisat-3n-build --repeat-index genome.fa genome` 

It is optional to make a graph index and add SNP or splice site information to the index, which can increase the alignment accuracy.
For more details, please refer to the [HISAT2 manual].

[HISAT2 manual]: {{ site.baseurl }}{% link _pages/manual.md %}

    # Standard HISAT-3N index with SNPs included
    hisat-3n-build --exons genome.exon genome.fa genome 
    
    # Standard HISAT-3N index with splice sites included
    hisat-3n-build --ss genome.ss genome.fa genome 
    
    # Repeat HISAT-3N index with SNPs included
    hisat-3n-build --repeat-index --exons genome.exon genome.fa genome 
    
    # Repeat HISAT-3N index with splice sites included
    hisat-3n-build --repeat-index --ss genome.ss genome.fa genome 


Alignment with `hisat-3n`
------------
After building the HISAT-3N index, you are ready to use `hisat-3n` for alignment. 
HISAT-3N has the same set of parameters as in HISAT2 with some additional arguments. Please refer to the [HISAT2 manual] for more details.

For the human reference genome, HISAT-3N requires about 9GB for alignment with the standard 3N-index and 10.5GB for the repeat 3N-index.

* `--base-change <nt1,nt2>`  
    Specify the nucleotide conversion type (e.g., C to T in bisulfite-sequencing reads). The parameter option is two characters separated by ','.  Type the original nucleotide for the first character (nt1) and type the converted nucleotide as the second character (nt2). For example, if performing [SLAM-seq] where some 'T's are converted to 'C's, input `--base-change T,C`. 
As another example, if performing bisulfite-seq, where some 'C's are converted to 'T's, please input `--base-change C,T`.
    If you want to align non-converted reads to the regular HISAT2 index, then omit this command.
       
* `--index/-x <hisat-3n-idx>`  
    Specify the index file basename for HISAT-3N.  The basename is the name of the index files up to but not including the suffix `.3n.*.*.ht2` / etc. 
    For example, if you build your index with basename 'genome' using a HISAT-3N-build, please input `--index genome`.
      
* `--repeat-limit <int>` 
    You can set up the number of alignments to be checked for each repeat alignment. You may increase the number to direct hisat-3n 
    to output more, if a read has multiple mapping locations. We suggest that you limit the repeat number for paired-end read alignment to no more 
    than 1,000,000. default: 1000.

* `--unique-only` 
    Only output uniquely aligned reads.

    
#### Examples:
* Single-end [SLAM-seq] read (T to C conversion) alignment with standard 3N-index:  
`hisat-3n --index genome -f -U read.fa -S output.sam --base-change T,C`

* Paired-end bisulfite-seq read (C to T conversion) alignment with repeat 3N-index:   
`hisat-3n --index genome -f -1 read_1.fa -2 read_2.fa -S output.sam --base-change C,T`

* Single-end TAPS reads (C to T conversion) alignment with repeat 3N-index and only output unique aligned results:   
`hisat-3n --index genome -q -U read.fq -S output.sam --base-change C,T --unique`



#### Extra SAM tags generated by HISAT-3N:

* `Yf:i:<N>`: Number of conversions detected in the read.

* `YZ:A:<A>`: The value `+` or `–` indicates the read is mapped to REF-3N (`+`) or REF-RC-3N (`-`), respectively.

Generate a 3N-conversion-table with `hisat-3n-table`
------------
### Preparation

To generate a 3N-conversion-table, users need to sort the `hisat-3n` generated SAM alignment file. 

[SAMtools] is required for this sorting process.

Use `samtools sort` to convert the SAM file into a sorted SAM file.

    samtools sort output.sam -o output_sorted.sam -O sam
    
Generate a 3N-conversion-table with `hisat-3n-table`:

### Usage
    hisat-3n-table [options]* --sam <samFile> --ref <refFile> --table-name <tableFile> --base-change <char1,char2>

#### Main arguments
* `--sam <samFile>`   
  Specify a sorted SAM filename

* `--ref <refFile>`  
  Specify the reference genome file (FASTA format), which was used to generate the HISAT-3N index. 
  
* `--table-name <tableFile>`  
  Specify the filename in which to write the 3N-conversion-table (tsv format).
  
* `--base-change <char1,char2>`  
  Specify the base-change rule. You should enter the exact same `--base-change` arguments that you used  in hisat-3n.
  For example, please input `--base-change C,T` for bisulfite-seq reads.
  
#### Input options
* `-u/--unique-only`  
  This directs the program to only count the unique aligned reads when constructing the 3N-conversion-table.
  
* `-m/--multiple-only`  
  This directs the program to only count the multiple aligned reads when constructing the 3N-conversion-table.
  
* `-c/--CG-only`  
  This directs the program to count the CpG islands in the reference genome. This option is especially designed for bisulfite-seq reads.
  
* `-p/--threads <int>`  
  This directs the program to launch these many `<int>` parallel threads for building the table (default: 1). 

* `-h/--help`  
   This will display the program’s usage information and then quit the program.


#### Examples:
* Generate a 3N conversion table for bisulfite sequencing data:
`hisat-3n-table -p 16 --sam output_sorted.sam --ref genome.fa --table-name output.tsv --base-change C,T`

* Generate a 3N-conversion-table for TAPS data and only count bases in CpG islands covered by uniquely aligned reads:  
`hisat-3n-table -p 16 --sam output_sorted.sam --ref genome.fa --table-name output.tsv --base-change C,T --CG-only --unique-only`
  

#### Note:
There are 7 columns in the 3N-conversion-table:

1. `ref`: the chromosome name.
2. `pos`: 1-based position in `ref`.
3. `strand`: '+' for forward strand. '-' for reverse strand.
4. `convertedBaseQualities`: the qualities of the converted bases in read-level measurement. The length of this string is equal to the number of converted bases.
5. `convertedBaseCount`: the number of distinct read positions where converted bases in read-level measurements were found.
this number is equal to the length of convertedBaseQualities.
6. `unconvertedBaseQualities`: the qualities of the unconverted bases in read-level measurement. The length of this string is equal to the number of unconverted bases in read-level measurement.
7. `unconvertedBaseCount`: the number of distinct read positions where unconverted bases in read-level measurements were found.
this number is equal to the length of unconvertedBaseQualities.

##### Sample 3N-conversion-table:
    ref    pos    strand    convertedBaseQualities    convertedBaseCount    unconvertedBaseQualities    unconvertedBaseCount
    1      11874  +         FFFFFB<BF<F               11                                                0
    1      11877  -         FFFFFF<                   7                                                 0
    1      11878  +         FFFBB//F/BB               11                                                0
    1      11879  +                                   0                     FFFBB//FB/                  10
    1      11880  -         F                         1                     FFFF/                       5



Download
========
You can download pre-built HISAT-3N indexes.  

- GRCh38

{:.table-noborder}
| :--- | :--- |
| genome_3n | <https://zenodo.org/record/4711458/files/grch38_3n.tar.gz?download=1> |
| genome_3n_rep | <https://zenodo.org/record/4711458/files/grch38_3n_rep.tar.gz?download=1> |
  
Publication
============

* HISAT-3N paper  
Yun Zhang, Chanhee Park, Christopher Bennett, Micah Thornton, and Daehwan Kim <br/>
[HISAT-3N: a rapid and accurate three-nucleotide sequence aligner](https://doi.org/10.1101/2020.12.15.422906). _bioRxiv_ (2020) 

* HISAT2 paper  
Daehwan Kim, Joseph Paggi, Chanhee Park, Christopher Bennett, and Steven Salzberg <br/>
[Graph-based genome alignment and genotyping with HISAT2 and HISAT-genotype](https://doi.org/10.1038/s41587-019-0201-4). _Nat Biotechnol_ **37**, 907–915 (2019) 
