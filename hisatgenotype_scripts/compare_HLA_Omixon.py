#!/usr/bin/env python

import sys, os
from argparse import ArgumentParser, FileType
use_message = '''
'''

def compare(hisatgenotype_fname, omixon_fname):
    hla_list = ["A", "B", "C", "DQA1", "DQB1", "DRB1"]
    
    # Read HISAT-genotype predicted HLA alleles for the CAAPA genomes
    hisat_hla = {}
    for line in open(hisatgenotype_fname):
        line = line.strip()
        fields = line.split('\t')
        if len(fields) == 2:
            sample, allele = fields
            abundance, vars_covered = 0.0, ""
        elif len(fields) == 3:
            sample, allele, abundance = fields
            vars_covered = ""
        else:
            assert len(fields) == 4
            sample, allele, abundance, vars_covered = fields
        abundance = float(abundance)
        if sample not in hisat_hla:
            hisat_hla[sample] = {}
        gene, allele = allele.split('*')
        if gene not in hisat_hla[sample]:
            hisat_hla[sample][gene] = []
        hisat_hla[sample][gene].append([allele, abundance])

    # Read Omixon predicted HLA alleles for the CAAPA genomes
    omixon_hla = {}
    for line in open(omixon_fname):
        line = line.strip()
        sample, allele1, allele2 = line.split('\t')
        gene1, allele1 = allele1.split('*')
        gene2, allele2 = allele2.split('*')
        
        assert gene1 == gene2
        if sample not in omixon_hla:
            omixon_hla[sample] = {}
        if gene1 not in omixon_hla[sample]:
            omixon_hla[sample][gene1] = []

        if len(omixon_hla[sample][gene1]) >= 2:
            continue
            
        omixon_hla[sample][gene1].append(allele1)
        omixon_hla[sample][gene1].append(allele2)

    for gene in hla_list:
        count, count_10 = [0, 0, 0], [0, 0, 0]
        print >> sys.stderr, gene
        for sample in omixon_hla.keys():
            if sample not in hisat_hla:
                continue
            hisat_sample = hisat_hla[sample]
            omixon_sample = omixon_hla[sample]
            if gene not in omixon_sample or gene not in hisat_sample:
                continue
            hisat_gene = hisat_sample[gene]
            omixon_gene = omixon_sample[gene]
            num_match, num_match_10 = 0, 0
            for omixon_allele in omixon_gene:
                omixon_allele = omixon_allele.split(':')
                for hisat_allele_idx in range(len(hisat_gene)):
                    hisat_allele = hisat_gene[hisat_allele_idx]
                    hisat_allele = hisat_allele[0].split(':')
                    equal = True
                    for i in range(min(len(omixon_allele), len(hisat_allele), 2)):
                        omixon_num = omixon_allele[i]
                        hisat_num = hisat_allele[i]
                        """
                        if not omixon_num[-1].isdigit():
                            omixon_num = omixon_num[:-1]
                        if not hisat_num[-1].isdigit():
                            hisat_num = hisat_num[:-1]
                        if int(hisat_num) != int(omixon_num):
                            equal = False
                            break
                        """
                        if hisat_num != omixon_num:
                            equal = False
                            break
                    if equal:
                        if hisat_allele_idx < 2:
                            num_match += 1
                        num_match_10 += 1
                        break
                    
            # DK - for debugging purposes
            """
            if gene in ["A", "B", "C", "DQA1", "DQB1", "DRB1"] and num_match < 2:
                print sample
                print "\t", omixon_gene
                print "\t", hisat_gene
                # sys.exit(1)
            """
                
            assert num_match < len(count)
            count[num_match] += 1
            count_10[num_match_10] += 1

        if sum(count) <= 0:
            continue
        
        print >> sys.stderr, "\tTop two\t0: %d, 1: %d, 2: %d (%.2f%%)" % (count[0], count[1], count[2], (count[1] + count[2] * 2) / float(sum(count) * 2) * 100.0)
        print >> sys.stderr, "\tTop ten\t0: %d, 1: %d, 2: %d (%.2f%%)" % (count_10[0], count_10[1], count_10[2], (count_10[1] + count_10[2] * 2) / float(sum(count_10) * 2) * 100.0)
        

if __name__ == "__main__":
    parser = ArgumentParser(
        description='Compare HISAT-genotype and Omixon HLA typing results')
    parser.add_argument('hisatgenotype_fname',
                        nargs='?',
                        type=str,
                        help='hisatgenotype file name (e.g. cp_hla.txt)')
    parser.add_argument('omixon_fname',
                        nargs='?',
                        type=str,
                        help='omixon file name (e.g. omixon_caapa_hla.txt)')

    args = parser.parse_args()

    compare(args.hisatgenotype_fname,
            args.omixon_fname)

