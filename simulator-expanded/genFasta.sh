#!/bin/bash
function usage {
	echo "./$(basename $0) -h --> Shows Usage"
}
verbose=false;
truerandom=false;
while getopts "g:n:s:o:vth" flag
do
	case "${flag}" in
		g) genfile=${OPTARG};;
		n) numread=${OPTARG};;
		s) readlen=${OPTARG};;
		o) outfile=${OPTARG};;
		v) verbose=true;;
		t) truerandom=true;;
		h) usage
		   exit;;
	esac
done

if $verbose; then
echo "Genome File: $genfile";
echo "Number of Reads: $numread"; 
echo "Read Length: $readlen"; 
echo "Output File: $outfile";  
fi;
if [[ ! -s $genfile ]]; then 
       echo "Error: File $genfile does not exist exiting"; 
       exit;
fi       

singleEntry=$(cat $genfile | grep ">" | wc -l);

if [[ ! $singleEntry =~ ^1$ ]]; then 
	echo "Error:  There are a total of $singleEntry > characters in $genfile";
	echo "        There should only be one header, and one assembled sequence";
        exit
fi

genomeFileLen=$(cat $genfile | tail -n +2 | tr -d '\n' |wc -c)
locsForReading=$(($genomeFileLen - $readlen))

if $verbose; then
echo "Genome File Length: $genomeFileLen Effective Number of Locations: $locsForReading";
fi;

