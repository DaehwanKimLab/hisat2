#!/bin/sh

#
# Downloads sequence for Zea mays AGPv3.31 from
# Ensembl.
#
# The base files, named ??.fa.gz
#
# By default, this script builds and index for just the base files,
# since alignments to those sequences are the most useful.  To change
# which categories are built by this script, edit the CHRS_TO_INDEX
# variable below.
#

ENSEMBL_RELEASE=31
ENSEMBL_DNA_BASE=ftp://ftp.ensemblgenomes.org/pub/plants/release-${ENSEMBL_RELEASE}/fasta/zea_mays/dna/
ENSEMBL_GTF_BASE=ftp://ftp.ensemblgenomes.org/pub/plants/release-${ENSEMBL_RELEASE}/gtf/zea_mays/
GTF_FILE=Zea_mays.AGPv3.31.gtf

ENSEMBL_VCF_BASE=ftp://ftp.ensemblgenomes.org/pub/plants/release-${ENSEMBL_RELEASE}/vcf/zea_mays/
VCF_FILE=zea_mays.vcf.gz
ERCC_FILE=ERCC92
ERCC_FTP=https://tools.thermofisher.com/content/sfs/manuals/${ERCC_FILE}.zip

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

HISAT2_VCF_SCRIPT=./hisat2_extract_snps_haplotypes_VCF.py
if [ ! -x "$HISAT2_VCF_SCRIPT" ] ; then
	if ! which hisat2_extract_snps_haplotypes_VCF.py ; then
		echo "Could not find hisat2_extract_snps_haplotypes_VCF.py in current directory or in PATH"
		exit 1
	else
		HISAT2_SNP_SCRIPT=`which hisat2_extract_snps_haplotypes_VCF.py`
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
F=Zea_mays.AGPv3.31.dna.genome.fa
if [ ! -f $F ] ; then
	get ${ENSEMBL_DNA_BASE}/$F.gz || (echo "Error getting $F" && exit 1)
	gunzip $F.gz || (echo "Error unzipping $F" && exit 1)
	mv $F genome.fa
	get ${ERCC_FTP} || (echo "Error getting ${ERCC_FILE}.zip" && exit 1)
	unzip ${ERCC_FILE}.zip
	cat ${ERCC_FILE}.fa >> genome.fa
fi

if [ ! -f $GTF_FILE ] ; then
       get ${ENSEMBL_GTF_BASE}/${GTF_FILE}.gz || (echo "Error getting ${GTF_FILE}" && exit 1)
       gunzip ${GTF_FILE}.gz || (echo "Error unzipping ${GTF_FILE}" && exit 1)
       cat ${ERCC_FILE}.gtf >> ${GTF_FILE}
       ${HISAT2_SS_SCRIPT} ${GTF_FILE} > genome.ss
       ${HISAT2_EXON_SCRIPT} ${GTF_FILE} > genome.exon
fi

if [ ! -f $VCF_FILE ] ; then
       get ${ENSEMBL_VCF_BASE}/${VCF_FILE} || (echo "Error getting ${ENSEMBL_VCF_BASE}/${VCF_FILE}" && exit 1)
       ${HISAT2_VCF_SCRIPT} --non-rs genome.fa ${VCF_FILE} genome
fi

CMD="${HISAT2_BUILD_EXE} -p 4 genome.fa --snp genome.snp --haplotype genome.haplotype --ss genome.ss --exon genome.exon genome_snp_tran_ercc"
echo Running $CMD
if $CMD ; then
	echo "genome index built; you may remove fasta files"
else
	echo "Index building failed; see error message"
fi
