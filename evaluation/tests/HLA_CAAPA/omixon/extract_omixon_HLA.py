#!/usr/bin/env python

import sys

if __name__ == '__main__':
    line = open("hla_omixon_hisat_analysis/data/raw/OMIXON_HLAoutput.txt").readline()
    hla_types = {}
    genome_id = ""
    for column in line.strip().split():
        if column.startswith("LP"):
            genome_id = column
            assert genome_id not in hla_types
            hla_types[genome_id] = []
        else:
            if column.startswith("A") or \
                    column.startswith("B") or \
                    column.startswith("C") or \
                    column.startswith("DR") or \
                    column.startswith("DQ"):
                hla_types[genome_id].append(column)

    for genome_id, types in hla_types.items():
        print genome_id
        print "\t", ' '.join(types)
            
