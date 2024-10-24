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
There are two strategies for HISAT-3N to align nucleotide conversion sequencing reads: *standard mode* and *repeat mode*.
The standard mode align reads with standard-3N index only, so it is fast and require small memory (~9GB for human genome alignment).
The repeat mode align reads with both standard-3N index and repeat-3N index, then output 1,000 alignment result (the output number can be adjusted by `--repeat-limit`).
The repeat mode can align nucleotide conversion reads more accurately,
and it is only 10% slower than the standard mode with tiny more memory (repeat mode use about ~10.5GB) usage than the standard mode.

HISAT-3N is developed based on [HISAT2], which is particularly optimized for RNA sequencing technology.
HISAT-3N supports both strand-specific and non-strand reads.
HISAT-3N can be used for any base-converted sequencing reads include [BS-seq], [SLAM-seq], [TAB-seq], [oxBS-seq], [TAPS], [scBS-seq], and [scSLAM-seq].

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
4. [SAMtools] is required to sort SAM file for `hisat-3n-table`.

Install
------------

    git clone https://github.com/DaehwanKimLab/hisat2.git hisat-3n
    cd hisat-3n
    git checkout -b hisat-3n origin/hisat-3n
    make


Make sure that you are in the `hisat-3n` branch


Build a HISAT-3N index with `hisat-3n-build`
-----------
`hisat-3n-build` builds a 3N-index, which contains two hisat2 indexes, from a set of DNA sequences. For standard 3N-index,
each index contains 16 files with suffix `.3n.*.*.ht2`.
For repeat 3N-index, there are 16 more files in addition to the standard 3N-index, and they have the suffix
`.3n.*.rep.*.ht2`.
These files constitute the hisat-3n index and no other file is needed to alignment reads to the reference.

* `--base-change <chr1,chr2>` argument is required for `hisat-3n-build` and `hisat-3n`.   
  Provide which base is converted in the sequencing process to another base. Please enter
  2 letters separated by ',' for this argument. The first letter(chr1) should be the converted base, the second letter(chr2) should be
  the converted to base. For example, during slam-seq, some 'T' is converted to 'C',
  please enter `--base-change T,C`. During bisulfite-seq, some 'C' is converted to 'T', please enter `--base-change C,T`. 
* Different conversion types may build the same hisat-3n index. Please check the table below for more detail. 
  Once you build the hisat-3n index with C to T conversion (for example BS-seq). 
  You can align the T to C conversion reads (for example SLAM-seq reads) with the same index.


  | Conversion Types                   | HISAT-3N index suffix         |
  |:----------------------------------:|:-----------------------------:|
  |C -> T<br>T -> C<br>A -> G<br>G -> A|.3n.CT.\*.ht2 <br>.3n.GA.\*.ht2|
  |A -> C<br>C -> A<br>G -> T<br>T -> G|.3n.AC.\*.ht2 <br>.3n.TG.\*.ht2|
  |A -> T<br>T -> A                    |.3n.AT.\*.ht2 <br>.3n.TA.\*.ht2|
  |C -> G<br>G -> C                    |.3n.CG.\*.ht2 <br>.3n.GC.\*.ht2|

#### Examples:
    # Build the standard HISAT-3N index (with C to T conversion):  
    hisat-3n-build --base-change C,T genome.fa genome
  
    # Build the repeat HISAT-3N index (with T to C conversion, require 256 GB memory for human genome index):  
    hisat-3n-build --base-change T,C --repeat-index genome.fa genome

It is optional to make the graph index and add SNP or spice site information to the index, to increase the alignment accuracy. 
The graph index building may require more memory than the linear index building.
For more detail, please check the [HISAT2 manual].

[HISAT2 manual]:https://daehwankimlab.github.io/hisat2/manual/

#### Examples:
    # Build the standard HISAT-3N index integrated index with SNP information
    hisat-3n-build --base-change C,T --snp genome.snp genome.fa genome 
    
    # Build the standard HISAT-3N integrated index with splice site information
    hisat-3n-build --base-change C,T --ss genome.ss --exon genome.exon genome.fa genome 
    
    # Build the repeat HISAT-3N index integrated index with SNP information
    hisat-3n-build --base-change C,T --repeat-index --snp genome.snp genome.fa genome 
    
    # Build the repeat HISAT-3N integrated index with splice site information
    hisat-3n-build --base-change C,T --repeat-index --ss genome.ss --exon genome.exon genome.fa genome 

Alignment with `hisat-3n`
------------
After we build the HISAT-3N index, you are ready to use `hisat-3n` for alignment.
HISAT-3N uses the HISAT2 argument but has some extra arguments. Please check [HISAT2 manual] for more detail.

* `--base-change <chr1,chr2>`  
  Provide which base is converted in the sequencing process to another base. Please enter
  2 letters separated by ',' for this argument. The first letter(chr1) should be the converted base, the second letter(chr2) should be
  the converted to base. For example, during slam-seq, some 'T' is converted to 'C',
  please enter `--base-change T,C`. During bisulfite-seq, some 'C' is converted to 'T', please enter `--base-change C,T`.

* `-x <hisat-3n-idx>`  
  The index for HISAT-3N.  The basename is the name of the index files up to but not including the suffix `.3n.*.*.ht2` / etc.
  For example, you build your index with basename 'genome' by HISAT-3N-build, please enter `-x genome`.

* `--directional-mapping`  
  Make directional mapping. Please use this option only if your sequencing reads are generated from a strand-specific library. 
  The directional mapping mode is about 2x faster than the default (non-directional) mapping mode.

* `--repeat-limit <int>`
  You can set up the number of alignment will be checked for each repeat alignment. You may increase the number to let hisat-3n
  output more, if a read has multiple mapping. We suggest the repeat limit number for paired-end reads alignment is no more
  than 1,000,000. default: 1000.

* `--unique-only`
  Only output uniquely aligned reads. This option may conflict with the hisat2 option `--un`, `--al`,
  `--un-conc`, and `--al-conc`. We will fix this problem in future version. Please do not use `--unique-only` if you want to use the hisat2 options above.

#### Examples:
    # Single-end slam-seq reads (T to C conversion, RNA) alignment with the standard 3N-index:  
      hisat-3n -x genome -f -U read.fa -S alignment_result.sam --base-change T,C --no-repeat-index
    
    # Paired-end strand-specific bisulfite-seq read (C to T conversion) alignment with repeat 3N-index:   
      hisat-3n -x genome -f -1 read_1.fa -2 read_2.fa -S alignment_result.sam --base-change C,T --repeat --no-spliced-alignment --directional-mapping
    
    # Single-end TAPS reads (have C to T conversion， RNA) alignment with the repeat 3N-index:   
      hisat-3n -x genome -q -U read.fq -S alignment_result.sam --base-change C,T --repeat


#### Extra SAM tags generated by HISAT-3N:

* `Yf:i:<N>`: Number of conversions are detected in the read.
* `Zf:i:<N>`: Number of un-converted bases are detected in the read. `Yf` + `Zf` = total number of bases which can be converted in the read sequence.
* `YZ:A:<A>`: The value `+` or `–` indicate the read is mapped to REF-3N (`+`) or REF-RC-3N (`-`).

Generate a 3N-conversion-table with `hisat-3n-table`
------------
### Preparation

To generate 3N-conversion-table, users need to sort the SAM file which generated by `hisat-3n`.
[SAMtools] is required for this sorting process.

Use `samtools sort` to convert the SAM file to a sorted SAM file.

    samtools sort alignment_result.sam -o sorted_alignment_result.sam -O sam

Generate 3N-conversion-table with `hisat-3n-table`:

### Usage
    hisat-3n-table [options]* --alignments <alignmentFile> --ref <refFile> --base-change <char1,char2>

#### Main arguments
* `--alignments <alignmentFile>`   
  SORTED SAM file. Please enter `-` for standard input.

* `--ref <refFile>`  
  The reference genome file (FASTA format) for generating HISAT-3N index.
  
* `--output-name <outputFile>`  
  Filename to write 3N-conversion-table (tsv format) to.  By default, table is written to the “standard out” or “stdout” filehandle (i.e. the console).

* `--base-change <char1,char2>`  
  The base-change rule. User should enter the exact same `--base-change` arguments in hisat-3n.
  For example, please enter `--base-change C,T` for bisulfite sequencing reads.

#### Input options
* `-u/--unique-only`  
  Only count the unique aligned reads into 3N-conversion-table.

* `-m/--multiple-only`  
  Only count the multiple aligned reads into 3N-conversion-table.

* `-c/--CG-only`  
  Only count the CpG sites in reference genome. This option is designed for bisulfite sequencing reads.

* `--added-chrname`  
  Please add this option if you use `--add-chrname` during `hisat-3n` alignment. 
  During `hisat-3n` alignment, the prefix "chr" is added in front of chromosome name and shows on SAM output, when user choose `--add-chrname`.
  `hisat-3n-table` cannot find the chromosome name on reference because it has an additional "chr" prefix. This option is to help `hisat-3n-table`
  find the matching chromosome name on reference file. The 3n-table provides the same chromosome name as SAM file.

* `--removed-chrname`  
  Please add this option if you use `--remove-chrname` during `hisat-3n` alignment.
  During `hisat-3n` alignment, the prefix "chr" is removed in front of chromosome name and shows on SAM output, when user choose `--remove-chrname`.
  `hisat-3n-table` cannot find the chromosome name on reference because it has no "chr" prefix. This option is to help `hisat-3n-table`
  find the matching chromosome name on reference file. The 3n-table provides the same chromosome name as SAM file.

#### Other options:
* `-p/--threads <int>`  
  Launch `int` parallel threads (default: 1) for table building.
  
* `-h/--help`  
  Print usage information and quit.

#### Examples:
    # Generate the 3N-conversion-table for bisulfite sequencing data:  
      hisat-3n-table -p 16 --alignments sorted_alignment_result.sam --ref genome.fa --output-name output.tsv --base-change C,T
    
    # Generate the 3N-conversion-table for TAPS data and only count base in CpG site and uniquely aligned:  
      hisat-3n-table -p 16 --alignments sorted_alignment_result.sam --ref genome.fa --output-name output.tsv --base-change C,T --CG-only --unique-only
    
    # Generate the 3N-conversion-table for bisulfite sequencing data from sorted BAM file:  
      samtools view -h sorted_alignment_result.bam | hisat-3n-table --ref genome.fa --alignments - --output-name output.tsv --base-change C,T
    
    # Generate the 3N-conversion-table for bisulfite sequencing data from unsorted BAM file:  
      samtools sort alignment_result.bam -O sam | hisat-3n-table --ref genome.fa --alignments - --output-name output.tsv --base-change C,T


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

* HISAT-3N   
  Zhang, Y., Park, C., Bennett, C., Thornton, M. and Kim, D. [Rapid and accurate alignment of nucleotide conversion sequencing reads with HISAT-3N](https://doi.org/10.1101/gr.275193.120). _Genome Research_ **31(7)**: 1290-1295 (2021)


* HIAST2   
  Kim, D., Paggi, J.M., Park, C. _et al._ [Graph-based genome alignment and genotyping with HISAT2 and HISAT-genotype](https://doi.org/10.1038/s41587-019-0201-4). _Nat Biotechnol_ **37**, 907–915 (2019)  
