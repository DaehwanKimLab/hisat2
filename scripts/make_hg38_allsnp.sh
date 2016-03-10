#!/bin/sh

#
# Downloads sequence for the HG38 version of H. spiens (human) from
# UCSC.
#
# The base files, named ??.fa.gz
#
# By default, this script builds and index for just the base files,
# since alignments to those sequences are the most useful.  To change
# which categories are built by this script, edit the CHRS_TO_INDEX
# variable below.
#

UCSC_HG38_BASE=http://hgdownload.cse.ucsc.edu/goldenPath/hg38/bigZips
F=hg38.chromFa.tar.gz

DBSNP_RELEASE=144
SNP_FILE=snp${DBSNP_RELEASE}.txt
UCSC_SNP=http://hgdownload.cse.ucsc.edu/goldenPath/hg38/database/${SNP_FILE}

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

HISAT2_SNP_SCRIPT=./hisat2_extract_snps_haplotypes_UCSC.py
if [ ! -x "$HISAT2_SNP_SCRIPT" ] ; then
	if ! which hisat2_extract_snps.py ; then
		echo "Could not find hisat2_extract_snps_haplotypes_UCSC.py in current directory or in PATH"
		exit 1
	else
		HISAT2_SNP_SCRIPT=`which hisat2_extract_snps_haplotypes_UCSC.py`
	fi
fi

rm -f genome.fa
get ${UCSC_HG38_BASE}/$F || (echo "Error getting $F" && exit 1)
tar xvzf $F || (echo "Error unzipping $F" && exit 1)
for i in {1..22}; do cat chroms/chr$i.fa >> genome.fa; done
cat chroms/chr[XYM].fa >> genome.fa 
rm $F

if [ ! -f $SNP_FILE ] ; then
       get ${UCSC_SNP}.gz || (echo "Error getting ${UCSC_COMMON_SNP}" && exit 1)
       ${HISAT2_SNP_SCRIPT} genome.fa ${SNP_FILE}.gz genome
fi

CMD="${HISAT2_BUILD_EXE} -p 4 --large-index genome.fa --snp genome.snp --haplotype genome.haplotype genome"
echo Running $CMD
if $CMD ; then
	echo "genome index built; you may remove fasta files"
else
	echo "Index building failed; see error message"
fi
