#!/bin/sh

#
# Downloads sequence for the T2T-CHM13v2.0 of H. sapiens (human) from NCBI. 
# ANNOTATION REPORT: 
# https://www.ncbi.nlm.nih.gov/genome/annotation_euk/Homo_sapiens/GCF_009914755.1-RS_2023_03
# 

GENOME_RELEASE=GCF_009914755.1_T2T-CHM13v2.0_genomic
GENOME_BASE=https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/009/914/755/GCF_009914755.1_T2T-CHM13v2.0/
GTF_FILE=$GENOME_RELEASE.gtf

get() {
	file=$1
	if ! wget --version >/dev/null 2>/dev/null ; then
		if ! curl --version >/dev/null 2>/dev/null ; then
			echo "Please install wget or curl somewhere in your PATH"
			exit 1
		fi
		curl -o `basename $1` $1
		return $?
	else
		wget $1
		return $?
	fi
}

HISAT2_BUILD_EXE=./hisat2-build
if [ ! -x "$HISAT2_BUILD_EXE" ] ; then
	if ! which hisat2-build ; then
		echo "Could not find hisat2-build in current directory or in PATH"
		exit 1
	else
		HISAT2_BUILD_EXE=`which hisat2-build`
	fi
fi

HISAT2_SS_SCRIPT=./hisat2_extract_splice_sites.py
if [ ! -x "$HISAT2_SS_SCRIPT" ] ; then
	if ! which hisat2_extract_splice_sites.py ; then
		echo "Could not find hisat2_extract_splice_sites.py in current directory or in PATH"
		exit 1
	else
		HISAT2_SS_SCRIPT=`which hisat2_extract_splice_sites.py`
	fi
fi

HISAT2_EXON_SCRIPT=./hisat2_extract_exons.py
if [ ! -x "$HISAT2_EXON_SCRIPT" ] ; then
	if ! which hisat2_extract_exons.py ; then
		echo "Could not find hisat2_extract_exons.py in current directory or in PATH"
		exit 1
	else
		HISAT2_EXON_SCRIPT=`which hisat2_extract_exons.py`
	fi
fi

rm -f genome.fa
F=$GENOME_RELEASE.fna
if [ ! -f $F ] ; then
	get ${GENOME_BASE}/$F.gz || (echo "Error getting $F" && exit 1)
	gunzip $F.gz || (echo "Error unzipping $F" && exit 1)
else
	cp $F genome.fa
fi


if [ ! -f $GTF_FILE ] ; then
	get ${GENOME_BASE}/${GTF_FILE}.gz || (echo "Error getting ${GTF_FILE}" && exit 1)
	gunzip ${GTF_FILE}.gz || (echo "Error unzipping ${GTF_FILE}" && exit 1)
else
       ${HISAT2_SS_SCRIPT} ${GTF_FILE} > genome.ss
       ${HISAT2_EXON_SCRIPT} ${GTF_FILE} > genome.exon
fi

${HISAT2_BUILD_EXE} -p 4 genome.fa --ss genome.ss --exon genome.exon genome_tran

