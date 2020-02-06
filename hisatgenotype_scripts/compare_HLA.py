#!/usr/bin/env python

import sys, os
from argparse import ArgumentParser, FileType
use_message = '''
'''

def compare(hisatgenotype_fname,
            utsw_fname):
    hla_list = ["A", "B", "C", "DQA1", "DQB1", "DRB1"]
    for level in [1,2]:
        print >> sys.stderr, "Level: %d" % level
        def read_hla_types(fname):
            hla, hla_orig = {}, {}
            for line in open(fname):
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
                # sample = sample.split('_')[0]
                abundance = float(abundance)
                if sample not in hla:
                    hla[sample] = {}
                    hla_orig[sample] = {}
                gene, allele = allele.split('*')
                if gene not in hla[sample]:
                    hla[sample][gene] = []
                    hla_orig[sample][gene] = []
                hla_orig[sample][gene].append([allele, abundance])

                if level == 1:
                    allele = allele.split(':')[0]
                else:
                    assert level == 2
                    allele = ':'.join(allele.split(':')[:2])

                found = False
                for i in range(len(hla[sample][gene])):
                    cmp_allele, cmp_abundance = hla[sample][gene][i]
                    if level == 1 or allele.find(':') == -1:
                        one = two = allele
                        cmp_one = cmp_two = cmp_allele
                    else:
                        one, two = allele.split(':')
                        cmp_one, cmp_two = cmp_allele.split(':')
                    if one == cmp_one and two == cmp_two:
                        found = True
                        hla[sample][gene][i][1] = cmp_abundance + abundance
                        break

                if not found:
                    hla[sample][gene].append([allele, abundance])

            for sample_hla in hla.values():
                for gene, allele_list in sample_hla.items():
                    sample_hla[gene] = sorted(allele_list, key=lambda a: a[1], reverse=True)
                
            return hla, hla_orig
                    
        hla1, hla1_orig = read_hla_types(hisatgenotype_fname)
        hla2, hla2_orig = read_hla_types(utsw_fname)

        for gene in hla_list:
            count, count_10 = [0, 0, 0], [0, 0, 0]
            print >> sys.stderr, "\t%s" % gene
            for sample in hla2.keys():
                if sample not in hla1:
                    continue
                hla1_sample = hla1[sample]
                hla2_sample = hla2[sample]
                if gene not in hla1_sample or gene not in hla2_sample:
                    continue
                hla1_gene = hla1_sample[gene]
                hla2_gene = hla2_sample[gene]
                num_match, num_match_10 = 0, 0
                for hla2_allele, _ in hla2_gene:
                    hla2_allele = hla2_allele.split(':')
                    for allele_idx in range(len(hla1_gene)):
                        hla1_allele = hla1_gene[allele_idx][0]
                        hla1_allele = hla1_allele.split(':')
                        equal = True
                        for i in range(min(len(hla1_allele), len(hla2_allele), level)):
                            hla1_num = hla1_allele[i]
                            hla2_num = hla2_allele[i]
                            if hla1_num != hla2_num:
                                equal = False
                                break
                            
                        if equal:
                            if allele_idx < 2:
                                num_match += 1
                                if len(hla2_gene) == 1:
                                    num_match += 1
                            num_match_10 += 1
                            if len(hla2_gene) == 1:
                                num_match_10 += 1
                            break

                # DK - for debugging purposes
                # """
                # if gene in ["A", "B", "C", "DQA1", "DQB1", "DRB1"] and num_match < 2:
                if level == 3 and gene in ["B"] and num_match < 2:
                    print sample
                    print "\t", hla1_gene, "orig:", hla1_orig[sample][gene]
                    print "\t", hla2_gene, "orig:", hla2_orig[sample][gene]
                    # sys.exit(1)
                # """

                # DK - debugging purposes
                if num_match >= len(count) or num_match_10 >= len(count_10):
                    print sample, num_match, num_match_10

                assert num_match < len(count) and num_match_10 < len(count_10)
                count[num_match] += 1
                count_10[num_match_10] += 1

            if sum(count) <= 0:
                continue

            print >> sys.stderr, "\t\tTop two\t0: %d, 1: %d, 2: %d (%.2f%%)" % (count[0], count[1], count[2], (count[1] + count[2] * 2) / float(sum(count) * 2) * 100.0)
            print >> sys.stderr, "\t\tTop ten\t0: %d, 1: %d, 2: %d (%.2f%%)" % (count_10[0], count_10[1], count_10[2], (count_10[1] + count_10[2] * 2) / float(sum(count_10) * 2) * 100.0)


if __name__ == "__main__":
    parser = ArgumentParser(
        description='Compare HISAT-genotype and Utsw HLA typing results')
    parser.add_argument('hisatgenotype_fname',
                        nargs='?',
                        type=str,
                        help='hisatgenotype file name (e.g. cp_hla.txt)')
    parser.add_argument('utsw_fname',
                        nargs='?',
                        type=str,
                        help='utsw file name (e.g. utsw_caapa_hla.txt)')

    args = parser.parse_args()

    compare(args.hisatgenotype_fname,
            args.utsw_fname)

