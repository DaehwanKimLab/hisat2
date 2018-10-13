#!/bin/bash

rm -r hisatgenotype_db/RBG
./hisat2-edit/hisatgenotype_scripts/hisatgenotype_extract_RBG.py

if [[ ! -d "RBG/" ]] ; then
    echo "RBG did not generate"
    exit
fi

mv RBG/ hisatgenotype_db/

./hisat2-edit/hisatgenotype_extract_vars.py --base rbg --verbose

./hisat2-edit/hisatgenotype_locus.py --base rbg -p 8 --debug "pair,full" 2> ABO_syn_results.txt

./hisatgenotype_results_parser.py > ABO_syn_results_counts3.txt
