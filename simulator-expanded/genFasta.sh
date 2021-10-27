#!/bin/bash
function usage {
	echo "Quickly Generate a set of Simulated Reads from a FASTA Assembly file." 
	echo "Usage: ./genFasta.sh -g [Genome (.fasta) File Name] -n [Number of Reads to Generate]" 
        echo "                     -s [Length of Reads to Generate] -o [Output File Name]         "
	echo "                     -v <Run in verbose mode> -t <Use true random values WIP>       "
        echo "                     -f <Filter out lines with only N at top of genome file>        "
        echo "                     -h <Show this Help Message>	                                  "
	echo "------------------------------------------------------------------------------------"
	echo "Version 1.0                                                                         "
}
verbose=false;
truerandom=false;
filter=false;
while getopts "g:n:s:o:vthf" flag
do
	case "${flag}" in
		g) genfile=${OPTARG};;
		n) numread=${OPTARG};;
		s) readlen=${OPTARG};;
		o) outfile=${OPTARG};;
		v) verbose=true;;
		t) truerandom=true;;
		f) filter=true;;
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

if [[ -e $outfile ]]; then 
	echo "Warning:  The file $outfile already exists"
	echo "remove it? (y) (n)";
	read -p '[y/n]: ' choice;
	while [[ ! $choice =~ ^[yn]$ ]]; 
	do
		read -p '[y/n]: ' choice;
	done
if [[ $choice =~ ^y$ ]]; then 
	echo "Removing $outfile ...";	
	rm -f $outfile;
else 
	echo "Please specify alternative file."; 
	exit
fi


fi

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

beggenfile=1;
charperline=$(($(head -n 2 $genfile | grep NNNNNNN$ | wc -c)));
linesN=$(($(cat $genfile | grep NNNNNNN$ | wc -l)));

if [[ $filter ]]; then 
	beggenfile=$(( $charperline * $linesN ));
fi


shuf -i $beggenfile-$(($genomeFileLen - $readlen)) -n $((numread)) > temp_locs.txt
STARTS=$(cat temp_locs.txt)

scaf=$(head $genfile |grep ">" | tr -d ">")
for START in $STARTS
do
	echo "$scaf	$START	$(($START+$readlen))" >> temp_locs_ranges.bed;
done;


if $verbose; then 
	echo "List of Locations Generated: locations_temp.txt";
fi


 bedtools getfasta -fi $genfile -bed temp_locs_ranges.bed > $outfile
 rm temp_locs.txt
 rm temp_locs_ranges.bed
