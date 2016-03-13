#!/bin/sh

#
# Downloads sequence for the rn6 version of R. norvegicus (rat) from
# UCSC.
#
# Note that UCSC's rn6 build has two categories of compressed fasta
# files:
#
# 1. The base files, named chr??.fa.gz
# 2. The unplaced-sequence files, named chr??_random.fa.gz
#
# By default, this script indexes all these files.  To change which
# categories are built by this script, edit the CHRS_TO_INDEX
# variable below.
#

RN6_BASE=ftp://hgdownload.cse.ucsc.edu/goldenPath/rn6/bigZips
F=rn6.fa.gz

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
get ${RN6_BASE}/$F || (echo "Error getting $F" && exit 1)
gzip -cd $F > genome.fa || (echo "Error unzipping $F" && exit 1)
rm $F

CMD="${HISAT2_BUILD_EXE} -p 4 genome.fa genome"
echo Running $CMD
if $CMD ; then
	echo "genome index built; you may remove fasta files"
else
	echo "Index building failed; see error message"
fi
