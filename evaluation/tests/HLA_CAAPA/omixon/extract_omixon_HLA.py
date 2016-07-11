#!/usr/bin/env python

import sys, os, subprocess
import glob

if __name__ == '__main__':
    HLA_types = {}
    stat_HLA_types = {}
    HLA_fnames = glob.glob("HLAresults/LP*gz")
    for HLA_fname in HLA_fnames:
        genome_id = HLA_fname.split('/')[1]
        genome_id = genome_id.split("_HLA")[0]
        assert genome_id not in HLA_types
        HLA_types[genome_id] = {}

        gzip_cmd = ["gzip", "-cd", HLA_fname]
        proc = subprocess.Popen(gzip_cmd, stdout=subprocess.PIPE, stderr=open("/dev/null", 'w'))
        expect_allele = False
        for line in proc.stdout:
            line = line.strip()
            if line == "":
                continue
            if line.startswith("Gene"):
                try:
                    gene = line.split()[1]
                except:
                    gene = ""
                continue
            if gene not in ["HLA-A", "HLA-B", "HLA-C", "HLA-DQA1", "HLA-DQB1", "HLA-DRB1"]:
                continue
            if line.startswith("Pair"):
                expect_allele = True
            if not line.startswith("HLA"):
                continue
            if not expect_allele:
                continue
            expect_allele = False
            if len(line.split()) >= 12:
                allele, allele2 = line.split()[0][:-1], line.split()[6][:-1]
            else:
                allele = line.split(',')[0]
                allele2 = allele

            if allele == "HLA-DRB1*08:01:03":
                allele = "HLA-DRB1*08:01:01"
            elif allele == "HLA-DRB1*11:11:02":
                allele = "HLA-DRB1*11:11:01"
            if allele2 == "HLA-DRB1*08:01:03":
                allele2 = "HLA-DRB1*08:01:01"
            elif allele2 == "HLA-DRB1*11:11:02":
                allele2 = "HLA-DRB1*11:11:01"

            if gene not in HLA_types[genome_id]:
                HLA_types[genome_id][gene] = []
            HLA_types[genome_id][gene].append([allele, allele2])

            if gene not in stat_HLA_types:
                stat_HLA_types[gene] = set()
            stat_HLA_types[gene].add(allele)
            stat_HLA_types[gene].add(allele2)

    for genome_id, types in HLA_types.items():
        # print genome_id
        for gene, gene_types in types.items():
            None
            # print "\t", gene, gene_types[0]

    #
    IMGT_gen = {}
    for fname in glob.glob("../IMGTHLA/fasta/*_gen.fasta"):
        for line in open(fname):
            line = line.strip()
            if not line.startswith('>'):
                continue
            allele = "HLA-" + line.split()[1]
            gene = allele.split('*')[0]
            if gene not in IMGT_gen:
                IMGT_gen[gene] = set()
            IMGT_gen[gene].add(allele)
    IMGT_partial = {}
    for fname in glob.glob("../IMGTHLA/fasta/*_nuc.fasta"):
        for line in open(fname):
            line = line.strip()
            if not line.startswith('>'):
                continue
            allele = "HLA-" + line.split()[1]
            gene = allele.split('*')[0]
            if gene not in IMGT_partial:
                IMGT_partial[gene] = set()
            IMGT_partial[gene].add(allele)

    #
    for gene, alleles in stat_HLA_types.items():
        num_gen, num_partial = 0, 0
        for allele in alleles:
            if allele in IMGT_gen[gene]:
                num_gen += 1
            else:
                allele1 = allele.split('*')[1].split(':')
                Subset = False
                for allele2 in IMGT_gen[gene]:
                    allele2 = allele2.split('*')[1].split(':')
                    Subset2 = True
                    for i in range(min(len(allele1), len(allele2))):
                        if allele1[i] != allele2[i]:
                            Subset2 = False
                            break

                    if Subset2:
                        Subset = True
                        break
                if Subset:
                    num_gen += 1
            if allele in IMGT_partial[gene]:
                num_partial += 1
            else:
                allele1 = allele.split('*')[1].split(':')
                Subset = False
                for allele2 in IMGT_partial[gene]:
                    allele2 = allele2.split('*')[1].split(':')
                    Subset2 = True
                    for i in range(min(len(allele1), len(allele2))):
                        if allele1[i] != allele2[i]:
                            Subset2 = False
                            break

                    if Subset2:
                        Subset = True
                        break
                if Subset:
                    num_partial += 1
                else:
                    print "%s is not in IMGT_partial[%s]" % (allele, gene)
                    for allele2 in IMGT_partial[gene]:
                        allele2 = allele2.split('*')[1].split(':')
                        Subset2 = True
                        for i in range(min(len(allele1), len(allele2))):
                            if allele1[i] != allele2[i]:
                                Subset2 = False
                                break

                            if i == 1:
                                print "\tCandidate: %s" % allele2
                                break
                

        print "%s: %d alleles %d/%d %d/%d" % (gene, len(alleles), num_gen, len(IMGT_gen[gene]), num_partial, len(IMGT_partial[gene]))
    
