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
is designed for nucleotide conversion sequencing technologies and implemented based on HISAT2. 
There are two strategies for HISAT-3N to align nuleotide conversion sequencing reads: *standard mode* and *repeat mode*. 
The standard mode align reads with standard-3N index only, so it is fast and require small memory (~9GB for human genome alignment).
The repeat mode align reads with both standard-3N index and repeat-3N index, then output 1,000 alignment result (the output number can be adjust by `--repeat-limit`).
The repeat mode can align nucleotide conversion reads more accurately, 
and it is only 10% slower than the standard mode with tiny more memory (repeat mode use about ~10.5GB) usage than standard mode.

HISAT-3N is developed based on [HISAT2], which is particularly optimized for RNA sequencing technology. 
HISAT-3N can be used for any base-converted sequencing reads include [BS-seq], [SLAM-seq], [TAB-seq], [oxBS-seq], [TAPS], [scBS-seq], and [scSLAM-seq],.

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
HISAT-3N requires a 64-bit computer running either Linux or Mac OS X and at least 16 GB of RAM. 

A few notes:  

1. The repeat 3N index building process requires 256 GB of RAM.
2. The standard 3N index building requires no more than 16 GB of RAM.
3. The alignment process with either standard or repeat index requires no more than 16 GB of RAM.
4. [SAMtools] is required to sort SAM file for hisat-3n-table.

Install
------------
   
    git clone https://github.com/DaehwanKimLab/hisat2.git
    cd hisat2
    git checkout -b hisat-3n origin/hisat-3n
    make


Make sure that you are in the `hisat-3n` branch


Build a 3N index with `hisat-3n-build`
-----------
`hisat-3n-build` builds a 3N-index, which contains two hisat2 indexes, from a set of DNA sequences. For standard 3N-index,
each index contains 16 files with suffix `.3n.*.*.ht2`.
For repeat 3N-index, there are 16 more files in addition to the standard 3N-index, and they have the suffix 
`.3n.*.rep.*.ht2`. 
These files constitute the hisat-3n index and no other file is needed to alignment reads to the reference.

* Example for standard HISAT-3N index building:  
`hisat-3n-build genome.fa genome`  

* Example for repeat HISAT-3N index building (require 256 GB memory):  
`hisat-3n-build --repeat-index genome.fa genome` 

It is optional to make the graph index and add SNP or spicing site information to the index, to increase the alignment accuracy.
for more detail, please check the [HISAT2 manual].

[HISAT2 manual]:https://daehwankimlab.github.io/hisat2/manual/

    # Standard HISAT-3N integrated index with SNP information
    hisat-3n-build --exons genome.exon genome.fa genome 
    
    # Standard HISAT-3N integrated index with splicing site information
    hisat-3n-build --ss genome.ss genome.fa genome 
    
    # Repeat HISAT-3N integrated index with SNP information
    hisat-3n-build --repeat-index --exons genome.exon genome.fa genome 
    
    # Repeat HISAT-3N integrated index with splicing site information
    hisat-3n-build --repeat-index --ss genome.ss genome.fa genome 

Alignment with `hisat-3n`
------------
After we build the HISAT-3N index, you are ready to use `hisat-3n` for alignment. 
HISAT-3N uses the HISAT2 argument but has some extra arguments. Please check [HISAT2 manual] for more detail.

For human genome reference, HISAT-3N requires about 9GB for alignment with standard 3N-index and 10.5 GB for repeat 3N-index.

* `--base-change <chr1,chr2>`  
    Provide which base is converted in the sequencing process to another base. Please enter
    2 letters separated by ',' for this argument. The first letter(chr1) should be the converted base, the second letter(chr2) should be
    the converted to base. For example, during slam-seq, some 'T' is converted to 'C',
    please enter `--base-change T,C`. During bisulfite-seq, some 'C' is converted to 'T', please enter `--base-change C,T`.
    If you want to align non-converted reads to the regular HISAT2 index, do not use this option.
       
* `--index/-x <hisat-3n-idx>`  
    The index for HISAT-3N.  The basename is the name of the index files up to but not including the suffix `.3n.*.*.ht2` / etc. 
    For example, you build your index with basename 'genome' by HISAT-3N-build, please enter `--index genome`.
      
* `--repeat-limit <int>` 
    You can set up the number of alignment will be check for each repeat alignment. You may increase the number to let hisat-3n 
    output more, if a read has multiple mapping. We suggest the repeat limit number for paired-end reads alignment is no more 
    than 1,000,000. default: 1000.

* `--unique-only` 
    Only output uniquely aligned reads.
    
#### Examples:
* Single-end slam-seq reads (T to C conversion) alignment with standard 3N-index:  
`hisat-3n --index genome -f -U read.fa -S output.sam --base-change T,C`

* Paired-end bisulfite-seq reads (C to T conversion) alignment with repeat 3N-index:   
`hisat-3n --index genome -f -1 read_1.fa -2 read_2.fa -S output.sam --base-change C,T`

* Single-end TAPS reads (have C to T conversion) alignment with repeat 3N-index and only output unique aligned result:   
`hisat-3n --index genome -q -U read.fq -S output.sam --base-change C,T --unique`



#### Extra SAM tags generated by HISAT-3N:

* `Yf:i:<N>`: Number of conversions are detected in the read.

* `YZ:A:<A>`: The value `+` or `–` indicate the read is mapped to REF-3N (`+`) or REF-RC-3N (`-`).

Generate a 3N-conversion-table with `hisat-3n-table`
------------
### Preparation

To generate 3N-conversion-table, users need to sort the SAM file which generated by `hisat-3n`. 
[SAMtools] is required for this sorting process.

Use `samtools sort` to convert the SAM file to a sorted SAM file.

    samtools sort output.sam -o output_sorted.sam -O sam
    
Generate 3N-conversion-table with `hisat-3n-table`:

### Usage
    hisat-3n-table [options]* --sam <samFile> --ref <refFile> --table-name <tableFile> --base-change <char1,char2>

#### Main arguments
* `--sam <samFile>`   
  The sorted SAM file processed by samtools.

* `--ref <refFile>`  
  The reference genome file (FASTA format) for generating HISAT-3N index. 
  
* `--table-name <tableFile>`  
  Filename to write 3N-conversion-table (tsv format) to.
  
* `--base-change <char1,char2>`  
  The base-change rule. User should enter the exact same `--base-change` arguments in hisat-3n.
  For example, please enter `--base-change C,T` for bisulfite sequencing reads.
  
#### Input options
* `-u/--unique-only`  
  Only count the unique aligned reads into 3N-conversion-table.
  
* `-m/--multiple-only`  
  Only count the multiple aligned reads into 3N-conversion-table.
  
* `-c/--CG-only`  
  Only count the CpG island in reference genome. This option is designed for bisulfite sequencing reads.
  
* `-p/--threads <int>`  
  Launch `int` parallel threads (default: 1) for table building. 

* `-h/--help`  
  Print usage information and quit.
  

#### Examples:
* Generate 3N conversion table for bisulfite sequencing data:
`hisat-3n-table -p 16 --sam output_sorted.sam --ref genome.fa --table-name output.tsv --base-change C,T`

* Generate 3N-conversion-table for TAPS data and only count base in CpG island and uniquely aligned:  
`hisat-3n-table -p 16 --sam output_sorted.sam --ref genome.fa --table-name output.tsv --base-change C,T --CG-only --unique-only`
  

#### Note:
There are 7 columns in the 3N-conversion-table:

1. `ref`: the chromosome name.
2. `pos`: 1-based position in ref.
3. `strand`: '+' for forward strand. '-' for reverse strand.
4. `convertedBaseQualities`: the qualities for converted base in read-level measurement. Length of this string is equal to
the number of converted Base in read-level measurement.
5. `convertedBaseCount`: number of distinct read positions where converted base in read-level measurements were found.
this number should equal to the length of convertedBaseQualities.
6. `unconvertedBaseQualities`: the qualities for unconverted base in read-level measurement. Length of this string is equal to
the number of unconverted Base in read-level measurement.
7. `unconvertedBaseCount`: number of distinct read positions where unconverted base in read-level measurements were found.
this number should equal to the length of unconvertedBaseQualities.

##### Sample 3N-conversion-table:
    ref    pos    strand    convertedBaseQualities    convertedBaseCount    unconvertedBaseQualities    unconvertedBaseCount
    1      11874  +         FFFFFB<BF<F               11                                                0
    1      11877  -         FFFFFF<                   7                                                 0
    1      11878  +         FFFBB//F/BB               11                                                0
    1      11879  +                                   0                     FFFBB//FB/                  10
    1      11880  -         F                         1                     FFFF/                       5
[SAMtools]:        http://samtools.sourceforge.net

Publication
============

* HISAT-3N paper

* HIAST2 paper  
Kim, D., Paggi, J.M., Park, C. _et al._ [Graph-based genome alignment and genotyping with HISAT2 and HISAT-genotype](https://doi.org/10.1038/s41587-019-0201-4). _Nat Biotechnol_ **37**, 907–915 (2019)  
