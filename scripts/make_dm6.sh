#!/bin/sh

#
# Downloads sequence for a D. melanogaster from flybase.  Currently set
# to download 5.22, but F, REL, and IDX_NAME can be edited to reflect a
# different version number.  (But note that you will usually also have
# to change the date in REL.)
#

DM6_BASE=ftp://hgdownload.cse.ucsc.edu/goldenPath/dm6/bigZips
F=dm6.fa.gz

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
		wget -O `basename $1` $1
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
get ${DM6_BASE}/$F || (echo "Error getting $F" && exit 1)
gzip -cd $F > genome.fa || (echo "Error unzipping $F" && exit 1)
rm $F

CMD="${HISAT2_BUILD_EXE} genome.fa genome"
echo "Running $CMD"
if $CMD ; then
	echo "genome index built; you may remove fasta files"
else
	echo "Index building failed; see error message"
fi

