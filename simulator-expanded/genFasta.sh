#!/bin/bash
function usage {
	echo "------------------------------------------------------------------------------------"
	echo "Quickly Generate a set of Simulated Reads from a FASTA Assembly file.               " 
	echo "------------------------------------------------------------------------------------"
	echo "------------------------------------------------------------------------------------"
	echo "Usage: ./genFasta.sh -g [Genome (.fasta) File Name] -n [Number of Reads to Generate]" 
        echo "                     -s [Length of Reads to Generate] -o [Output File Name]         "
	echo "                     -v <Run in verbose mode> -t <Use true random values WIP>       "
        echo "                     -f <Filter out lines with only N at top of genome file>        "
        echo "                     -h <Show this Help Message> -b <Generate Bisulfite Reads>	  "     
	echo "------------------------------------------------------------------------------------"
	echo "Version 1.0                                                                         "
	echo "------------------------------------------------------------------------------------"
}
verbose=false;
truerandom=false;
filter=false;
bisulfitereads=false;
while getopts "g:n:s:o:vthfb" flag
do
	case "${flag}" in
		g) genfile=${OPTARG};;
		n) numread=${OPTARG};;
		s) readlen=${OPTARG};;
		o) outfile=${OPTARG};;
		v) verbose=true;;
		t) truerandom=true;;
		f) filter=true;;
		b) bisulfitereads=true;;
		h) usage
		   exit;;
	esac
done

if $verbose; then
	echo "------------------------------------------------------------------------------------"
echo "Genome File: $genfile";
echo "Number of Reads: $numread"; 
echo "Read Length: $readlen"; 
echo "Output File: $outfile";  
echo "Generate Bisulfite Reads?: $bisulfitereads";
	echo "------------------------------------------------------------------------------------"
fi;

if [[ -e $outfile ]]; then 
	echo "------------------------------------------------------------------------------------"
	echo "Warning:  The file $outfile already exists"
	echo "remove it? (y) (n)";
	read -p '[y/n]: ' choice;
	while [[ ! $choice =~ ^[yn]$ ]]; 
	do
		read -p '[y/n]: ' choice;
	done
if [[ $choice =~ ^y$ ]]; then 
	echo "Removing $outfile ...";	
	echo "------------------------------------------------------------------------------------"
	rm -f $outfile;
else 
	echo "Please specify alternative file."; 
	echo "------------------------------------------------------------------------------------"
	exit
fi


fi

if [[ ! -s $genfile ]]; then 
	echo "------------------------------------------------------------------------------------"
       echo "Error: File $genfile does not exist exiting"; 
	echo "------------------------------------------------------------------------------------"
       exit;
fi       

singleEntry=$(cat $genfile | grep ">" | wc -l);

if [[ ! $singleEntry =~ ^1$ ]]; then 
	echo "------------------------------------------------------------------------------------"
	echo "Error:  There are a total of $singleEntry > characters in $genfile";
	echo "        There should only be one header, and one assembled sequence";
	echo "------------------------------------------------------------------------------------"
        exit
fi

genomeFileLen=$(cat $genfile | tail -n +2 | tr -d '\n' |wc -c)
locsForReading=$(($genomeFileLen - $readlen))

if $verbose; then
	echo "------------------------------------------------------------------------------------"
echo "Genome File Length: $genomeFileLen Effective Number of Locations: $locsForReading";
	echo "------------------------------------------------------------------------------------"
fi;

beggenfile=1;
charperline=$(($(head -n 2 $genfile | grep NNNNNNN$ | wc -c)));

if $filter; then 
	nlines=$( cat $genfile | tail -n +1 | grep -m1 -n [ACGT] | grep -o "[0-9]\+")
	beggenfile=$(( $charperline * $nlines ));
	if $verbose; then 
		echo "------------------------------------------------------------------------------------"
		echo "Filtering At the Beginning of the File (To Remove only N rows)                      "
		echo "First non-N only line: $nlines							  "
		echo "Characters per line: $charperline							  "
		echo "Starting Postion: $beggenfile							  "
		echo "------------------------------------------------------------------------------------"
	fi
fi


shuf -i $beggenfile-$(($genomeFileLen - $readlen)) -n $((numread)) > temp_locs.txt
STARTS=$(cat temp_locs.txt)

if [[ -e temp_locs_ranges.bed ]]; then 
	rm temp_locs_ranges.bed
fi

scaf=$(head $genfile |grep ">" | tr -d ">")
for START in $STARTS
do
	echo "$scaf	$START	$(($START+$readlen))" >> temp_locs_ranges.bed;
done;


if $verbose; then 
	echo "------------------------------------------------------------------------------------"
	echo "List of Locations Generated: locations_temp.txt";
	echo "------------------------------------------------------------------------------------"
fi

if $bisulfitereads; then 
	echo "------------------------------------------------------------------------------------"
	echo "Use Previously Generated Methylation Profile? (y) - yes (n) - no			  "
	read -p '[y/n]: ' methfilechoice;
	while [[ ! $methfilechoice =~ ^[yn]$ ]]; 
	do
		read -p '[y/n]: ' methfilechoice;
	done
	methpatternfile="temp_cpg_locs.txt";
	echo "Generating Random Methylation Profile...                                            ";
	cat $genfile | tail -n +2 | tr -d '\n' | grep -o -b "CG" | grep -Eo '[[:digit:]]+' > temp_cpg_locs.txt
        echo "Number of CpGs found in file: $(cat temp_cpg_locs.txt |wc -l)                       "	
	echo "Generate Random Methylation Profile for all locations? (y) - yes (n) -no            "; 
	read -p '[y/n]: ' methchoice;
	while [[ ! $methchoice =~ ^[yn]$ ]]; 
	do
		read -p '[y/n]: ' methchoice;
	done
	if [[ $methchoice =~ ^[n]$ ]]; then
		echo "Generate for How Many Locations (Max = $(cat temp_cpg_locs.txt | wc -l))?"
		read -p 'Number of Locations: ' nummethloc;
		while [[ ! $nummethloc =~ ^[0-9+] ]];
		do
			read -p 'Number of Locations: ' nummethloc;
		done
		cat temp_cpg_locs.txt | sort -R | head -n $nummethloc > temp_cpg_locs_select.txt
		echo "Randomly Selected $nummethloc of $(cat temp_cpg_locs.txt | wc -l)"
	        rm temp_cpg_locs.txt
		cp temp_cpg_locs_select.txt temp_cpg_locs.txt
	else 
		echo "Generating for Maximum Number of Locations"; 
		nummethloc=$(cat temp_cpg_locs.txt | wc -l)
	fi
	methlocs=$(cat temp_cpg_locs.txt);
	if [[ -e temp_cpg_locs_probs_from_R.csv ]]; then 
		rm temp_cpg_locs_probs_from_R.csv;
	fi;
	Rscript getMethylationProfile.R
	if [[ ! -s temp_cpg_locs_probs_from_R.csv ]]; then 
		echo "Error generating Methylation Profile, is Rscript installed?"; 
	fi	
	echo "------------------------------------------------------------------------------------"
fi
	echo "------------------------------------------------------------------------------------"
	echo "Generating the Reads File Now...                                                    "
 	bedtools getfasta -fi $genfile -bed temp_locs_ranges.bed > $outfile
 	echo "Reads Generated in file $outfile						          "
	echo "------------------------------------------------------------------------------------"

 methpatternfile="temp_cpg_locs_probs_from_R.csv";
 if $bisulfitereads; then 
	echo "------------------------------------------------------------------------------------"
	echo "Methylating Generated Reads File $outfile,                                          "
	echo "            according to methylation profile  $methpatternfile                      " 
	cat $outfile | grep '>' > temp_split1.txt
	cat $outfile | grep -v '>' > temp_split2.txt 
	paste temp_split1.txt temp_split2.txt -d, > temp_combo_reads.csv
	sleep 5
	rm Methylated_reads_out_of_R.csv
	Rscript methylateReads.R
	cat Methylated_reads_out_of_R.csv | tail -n +2 > methreadout.csv
	paste temp_combo_reads.csv methreadout.csv -d',' | grep -v ,NA | cut -f1,5 -d, | tr -d '\"' > methreads
	paste temp_combo_reads.csv methreadout.csv -d',' | grep ,NA | cut -f1,2 -d, | tr -d '\"' > nonmethreads
	cat methreads nonmethreads | sed 's/,/\n/' > $outfile
	rm methreads nonmethreads
	echo "Reads in file $outfile are now methylated according to ratios in $methpatternfile   " 
	echo "Saving methylation pattern file as meth_pattern.csv                                 " 
	cat $methpatternfile > meth_pattern.csv
	echo "------------------------------------------------------------------------------------"

fi
 rm temp_*
