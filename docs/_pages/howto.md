---
layout: page
title: HOWTO
permalink: /howto/
order: 5
share: false
---

## HOWTO
{: .no_toc}

- TOC
{:toc}

### Building indexes
Depend on your purpose, you have to download reference sequence, gene annotation and SNP files

#### Prepare data
1. Download reference
```
$ wget ftp://ftp.ensembl.org/pub/release-95/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz
$ gzip -d Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz
$ mv Homo_sapiens.GRCh38.dna.primary_assembly.fa genome.fa
```
   <p/>
   
1. Download GTF and make exon, splicesite file.  
   If you want to build HFM index, you can skip this step.
```
$ wget ftp://ftp.ensembl.org/pub/release-95/gtf/homo_sapiens/Homo_sapiens.GRCh38.95.gtf.gz  
$ gzip -d Homo_sapiens.GRCh38.95.gtf.gz
$ mv Homo_sapiens.GRCh38.95.gtf genome.gtf
$ hisat2_extract_splice_sites.py genome.gtf > genome.ss
$ hisat2_extract_exons.py genome.gtf  > genome.exon
```

1. Download SNP  
   If you want to build HFM index, you can skip this step.  
   ```
   Hello World3
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
$ hisat2_extract_splice_sites.py genome.gtf > genome.ss
$ hisat2_extract_exons.py genome.gtf > genome.exon
$ hisat2-build -p 16 --exon genome.exon --ss genome.ss genome.fa genome_tran
```

#### Build HGFM index with SNPs and transcripts

```
$ hisat2-build -p 16 --snp genome.snp --haplotype genome.haplotype --exon genome.exon --ss genome.ss genome.fa genome_snp_tran
```



