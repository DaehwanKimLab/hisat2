#!/bin/sh

#
# Downloads sequence for the GRCh37 release 75 version of H. spiens (human) from
# Ensembl.
#
# Note that Ensembl's GRCh37 build has three categories of compressed fasta
# files:
#
# The base files, named ??.fa.gz
#
# By default, this script builds and index for just the base files,
# since alignments to those sequences are the most useful.  To change
# which categories are built by this script, edit the CHRS_TO_INDEX
# variable below.
#

CHRS_TO_INDEX="$BASE_CHRS"

ENSEMBL_RELEASE=75
ENSEMBL_GRCh37_GTF_BASE=ftp://ftp.ensembl.org/pub/release-${ENSEMBL_RELEASE}/gtf/homo_sapiens
UCSC_HG19=http://hgdownload.cse.ucsc.edu/goldenPath/hg19/bigZips/chromFa.tar.gz

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

rm -f genome.fa
get ${ENSEMBL_GRCh37_BASE}/$F.gz || (echo "Error getting $F" && exit 1)
gunzip $F.gz || (echo "Error unzipping $F" && exit 1)
awk '{if($1 ~ /^>/) {print $1} else {print}}' $F >> genome.fa
rm $F

CMD="${HISAT2_BUILD_EXE} genome.fa genome"
echo Running $CMD
if $CMD ; then
	echo "genome index built; you may remove fasta files"
else
	echo "Index building failed; see error message"
fi
