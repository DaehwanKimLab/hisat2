---
layout: page
title: HowTo
permalink: /howto/
order: 6
share: false
---

## HOWTO
{: .no_toc}

- TOC
{:toc}

### Building indexes
Depend on your purpose, you have to download reference sequence, gene annotation and SNP files.  
We also provides scripts to build indexes. [Download]({{ site.baseurl }}{% link _pages/download.md %})

#### Prepare data
1. Download reference
```
$ wget ftp://ftp.ensembl.org/pub/release-84/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz
$ gzip -d Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz
$ mv Homo_sapiens.GRCh38.dna.primary_assembly.fa genome.fa
```
   
1. Download GTF and make exon, splicesite file.  
   If you want to build HFM index, you can skip this step.
```
$ wget ftp://ftp.ensembl.org/pub/release-84/gtf/homo_sapiens/Homo_sapiens.GRCh38.84.gtf.gz  
$ gzip -d Homo_sapiens.GRCh38.84.gtf.gz
$ mv Homo_sapiens.GRCh38.84.gtf genome.gtf
$ hisat2_extract_splice_sites.py genome.gtf > genome.ss
$ hisat2_extract_exons.py genome.gtf > genome.exon
```

1. Download SNP  
   If you want to build HFM index, you can skip this step.  
```
$ wget http://hgdownload.cse.ucsc.edu/goldenPath/hg38/database/snp144Common.txt.gz
$ gzip -d snp144Common.txt.gz
```

   Convert chromosome names of UCSC Database to Ensembl Annotation
```
$ awk 'BEGIN{OFS="\t"} {if($2 ~ /^chr/) {$2 = substr($2, 4)}; if($2 == "M") {$2 = "MT"} print}' snp144Common.txt > snp144Common.txt.ensembl
```

   make SNPs and haplotype file
```
$ hisat2_extract_snps_haplotypes_UCSC.py genome.fa snp144Common.txt.ensembl genome
```

#### Build HFM index
It takes about 20 minutes(depend on HW spec) to build index, and requires at least 6GB memory.
```
$ hisat2-build -p 16 genome.fa genome
```

#### Build HGFM index with SNPs
```
$ hisat2-build -p 16 --snp genome.snp --haplotype genome.haplotype genome.fa genome_snp
```

#### Build HGFM index with transcripts
It takes about 1 hour(depend on HW spec) to build index, and requires at least 160GB memory.
```
$ hisat2-build -p 16 --exon genome.exon --ss genome.ss genome.fa genome_tran
```

#### Build HGFM index with SNPs and transcripts

```
$ hisat2-build -p 16 --snp genome.snp --haplotype genome.haplotype --exon genome.exon --ss genome.ss genome.fa genome_snp_tran
```



