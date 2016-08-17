#!/usr/bin/env python

import sys, os
use_message = '''
'''

if __name__ == "__main__":
    hla_list = ["A", "B", "C", "DQA1", "DQB1", "DRB1"]
    
    # Read HISAT-genotype predicted HLA alleles for the CAAPA genomes
    hisat_hla = {}
    for line in open("hisat_caapa_hla.txt"):
        line = line.strip()
        sample, allele, abundance = line.split('\t')
        abundance = float(abundance)
        # DK - for debugging purposes
        if abundance < 5.0:
            continue
        if sample not in hisat_hla:
            hisat_hla[sample] = {}
        gene, allele = allele.split('*')
        if gene not in hisat_hla[sample]:
            hisat_hla[sample][gene] = []
        hisat_hla[sample][gene].append([allele, abundance])

    # Read Omixon predicted HLA alleles for the CAAPA genomes
    omixon_hla = {}
    for line in open("omixon_caapa_hla.txt"):
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
            if gene not in omixon_sample:
                continue
            assert gene in hisat_sample
            hisat_gene = hisat_sample[gene]
            omixon_gene = omixon_sample[gene]
            num_match, num_match_10 = 0, 0
            for omixon_allele in omixon_gene:
                omixon_allele = omixon_allele.split(':')
                for hisat_allele_idx in range(len(hisat_gene)):
                    hisat_allele = hisat_gene[hisat_allele_idx]
                    hisat_allele = hisat_allele[0].split(':')
                    equal = True
                    for i in range(min(len(omixon_allele), len(hisat_allele))):
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
            # """
            if num_match == 1 and num_match_10 == 2:
                print sample
                print "\t", omixon_gene
                print "\t", hisat_gene
                sys.exit(1)
            # """
                
            assert num_match < len(count)
            count[num_match] += 1
            count_10[num_match_10] += 1
            
        print >> sys.stderr, "\tTop two\t0: %d, 1: %d, 2: %d (%.2f%%)" % (count[0], count[1], count[2], (count[1] + count[2] * 2) / float(sum(count) * 2) * 100.0)
        print >> sys.stderr, "\tTop ten\t0: %d, 1: %d, 2: %d (%.2f%%)" % (count_10[0], count_10[1], count_10[2], (count_10[1] + count_10[2] * 2) / float(sum(count_10) * 2) * 100.0)
        
