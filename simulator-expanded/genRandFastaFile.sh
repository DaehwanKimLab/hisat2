#!/bin/bash
numReads=$1;
readLenConst=$2;
readLen=$3;
trueRandom=$4;

if [[ ! $# =~ ^4$ ]]; then 
	echo "Incorrect # of Arguments: $# Should be 3";
	echo "Usage: ./genRandFastaFile.sh [Number of Reads] [Read Len Constant (0)-no (1)-yes] [Read Length] [Use True Random Data /dev/urandom (0)-no (1)-yes]"
	exit 1; 
fi 
if [[ ! $numReads =~ ^[0-9]+$ ]]; then
	echo "Error in Argument: $1";
	echo "Argument one should be the number of reads."; 
	exit 1;
fi
if [[ ! $readLenConst =~ ^[0-1]$ ]]; then
	echo "Error in Argument: $2";
	echo "Argument two should be 0 or 1 indicating (no) or (yes) to constant read length";
	exit 1;
fi 
if [[ ! $trueRandom =~ ^[0-1]$ ]]; then
	echo "Error in Argument: $2";
	echo "Argument four should be 0 or 1 indicating (no) or (yes) to use of true randoms";
	exit 1;
fi 
if [[ ! $readLen =~ ^[0-9]+$ ]]; then 
	echo "Error in Argument: $3"; 
        echo "Argument three should be a value indicating either read length, or average read length"; 
	exit 1; 
fi	
for ((i=1;i < $numReads;i++));  
do
	printf ">> Read $i\n";
	for ((j=1;j < $readLen;j++));
	do
		alp="ACGT";
		if [[ $trueRandom =~ ^1$ ]]; then
		rnd=$(od -A n -t d -N 1 /dev/urandom);
	else
		rnd=$RANDOM;
	fi
		loc=$(($rnd%4+1));
		pck=${alp:$(($loc-1)):1}
		printf "$pck";
	done
	printf  "\n";
done
