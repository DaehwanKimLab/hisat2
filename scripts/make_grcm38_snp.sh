#!/bin/sh

#
# Downloads sequence for the GRCm37 release 81 version of M. Musculus (mouse) from
# Ensembl.
#
# By default, this script builds and index for just the base files,
# since alignments to those sequences are the most useful.  To change
# which categories are built by this script, edit the CHRS_TO_INDEX
# variable below.
#

ENSEMBL_RELEASE=81
ENSEMBL_GRCm38_BASE=ftp://ftp.ensembl.org/pub/release-${ENSEMBL_RELEASE}/fasta/mus_musculus/dna

DBSNP_RELEASE=138
SNP_FILE=snp${DBSNP_RELEASE}Common.txt
UCSC_COMMON_SNP=http://hgdownload.cse.ucsc.edu/goldenPath/mm10/database/${SNP_FILE}

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

HISAT2_SNP_SCRIPT=./extract_snps.py
if [ ! -x "$HISAT2_SNP_SCRIPT" ] ; then
	if ! which extract_snps.py ; then
		echo "Could not find extract_snps.py in current directory or in PATH"
		exit 1
	else
		HISAT2_SNP_SCRIPT=`which extract_snps.py`
	fi
fi

rm -f genome.fa
F=Mus_musculus.GRCm38.dna.primary_assembly.fa
if [ ! -f $F ] ; then
	get ${ENSEMBL_GRCm38_BASE}/$F.gz || (echo "Error getting $F" && exit 1)
	gunzip $F.gz || (echo "Error unzipping $F" && exit 1)
	awk '{if($1 ~ /^>/) {print $1} else {print}}' $F > genome.fa
	rm $F
fi


if [ ! -f $SNP_FILE ] ; then
       get ${UCSC_COMMON_SNP}.gz || (echo "Error getting ${UCSC_COMMON_SNP}" && exit 1)
       gunzip ${SNP_FILE}.gz || (echo "Error unzipping ${SNP_FILE}" && exit 1)
       awk 'BEGIN{OFS="\t"} {if($2 ~ /^chr/) {$2 = substr($2, 4)}; if($2 == "M") {$2 = "MT"} print}' ${SNP_FILE} > ${SNP_FILE}.tmp
       mv ${SNP_FILE}.tmp ${SNP_FILE}
       ${HISAT2_SNP_SCRIPT} genome.fa ${SNP_FILE} > genome.snp
fi

CMD="${HISAT2_BUILD_EXE} -p 4 genome.fa --snp genome.snp genome_snp"
echo Running $CMD
if $CMD ; then
	echo "genome index built; you may remove fasta files"
else
	echo "Index building failed; see error message"
fi
