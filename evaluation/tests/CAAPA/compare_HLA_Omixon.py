#!/usr/bin/env python

import sys, os
use_message = '''
'''

def compare_HLA():
    # Read HISAT-genotype predicted HLA alleles for the CAAPA genomes
    hisat_hla = {}
    for line in open("HISAT.txt"):
        line = line.strip()
        sample, allele, abundance = line.split('\t')
        if sample not in hisat_hla:
            hisat_hla[sample] = {}
        gene, allele = allele.split('*')
        if gene not in hisat_hla[sample]:
            hisat_hla[sample][gene] = []
        hisat_hla[sample][gene].append([allele, abundance])

    # Read Omixon predicted HLA alleles for the CAAPA genomes
    omixon_hla = {}
    for line in open("Omixon.txt"):
        line = line.strip()
        sample, allele1, allele2 = line.split('\t')
        gene1, gene2 = allele1.split(':')[0], allele2.split(':')[0]
        allele1, allele2 = allele1.split(':')[1:], allele2.split(':')[1:]
        allele1, allele2 = ':'.join(allele1), ':'.join(allele2)
        assert gene1 == gene2
        if sample not in omixon_hla:
            omixon_hla[sample] = {}
        if gene1 not in omixon_hla[sample]:
            omixon_hla[sample][gene1] = []
        omixon_hla[sample][gene1].append(allele1)
        omixon_hla[sample][gene1].append(allele2)

    for gene in ["A", "B", "C", "DQA1", "DQB1", "DRB1"]:
        count = [0, 0, 0]
        print >> sys.stderr, gene
        for sample in omixon_hla.keys():
            assert sample in hisat_hla
            hisat_sample = hisat_hla[sample]
            omixon_sample = omixon_hla[sample]
            if gene not in omixon_sample:
                continue
            assert gene in hisat_sample
            hisat_gene = hisat_sample[gene]
            omixon_gene = omixon_sample[gene]
            num_match = 0
            for omixon_allele in omixon_gene:
                omixon_allele = omixon_allele.split(':')
                for hisat_allele_idx in range(len(hisat_gene)):
                    # daehwan - for testing purposes
                    if hisat_allele_idx >= 2:
                        break
                    hisat_allele = hisat_gene[hisat_allele_idx]
                    hisat_allele = hisat_allele[0].split(':')
                    equal = True
                    for i in range(min(len(omixon_allele), len(hisat_allele))):
                        omixon_num = omixon_allele[i]
                        hisat_num = hisat_allele[i]
                        if not hisat_num[-1].isdigit():
                            hisat_num = hisat_num[:-1]
                        if int(hisat_num) != int(omixon_num):
                            equal = False
                            break
                    if equal:
                        num_match += 1
                        break
            assert num_match < len(count)
            count[num_match] += 1
        print >> sys.stderr, "\t0: %d, 1: %d, 2: %d" % (count[0], count[1], count[2])
        
    
if __name__ == "__main__":
    compare_HLA()
