#!/usr/bin/env python

#
# Copyright 2016, Daehwan Kim <infphilo@gmail.com>
#
# This file is part of HISAT 2.
#
# HISAT 2 is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# HISAT 2 is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with HISAT 2.  If not, see <http://www.gnu.org/licenses/>.
#


import sys, os, subprocess, glob

if __name__ == '__main__':
    hla_list = ["A", "B", "C", "DQA1", "DQB1", "DRB1"]
    gen_alleles = {}
    for hla in hla_list:
        for line in open("IMGTHLA/fasta/%s_gen.fasta" % hla):
            if line.startswith(">"):
                allele = line.split()[1]
                gene = allele.split('*')[0]
                if gene not in gen_alleles:
                    gen_alleles[gene] = set()
                gen_alleles[gene].add(allele)
                
    nuc_alleles = {}
    for hla in hla_list:
        for line in open("IMGTHLA/fasta/%s_nuc.fasta" % hla):
            if line.startswith(">"):
                allele = line.split()[1]
                gene = allele.split('*')[0]
                if gene not in nuc_alleles:
                    nuc_alleles[gene] = set()
                nuc_alleles[gene].add(allele)

    print >> sys.stderr, "IMGTHLA"
    for gene, alleles in nuc_alleles.items():
        print >> sys.stderr, "\t%s: %d alleles" % (gene, len(alleles))

    # Read HLA alleles from Omixon data
    omixon_alleles = {}
    omixon_fnames = glob.glob("HLAresults/*.gz")
    for fname in omixon_fnames:
        genome = fname.split("/")[1].split("_HLA")[0]
        view_cmd = ["gzip", "-cd", fname]
        proc = subprocess.Popen(view_cmd, stdout=subprocess.PIPE, stderr=open("/dev/null", 'w'))
        allele_count = {}
        prev_allele1, prev_allele2 = "", ""
        for line in proc.stdout:
            if not line.startswith("HLA"):
                continue

            fields = line.strip().split()
            if len(fields) > 6:
                allele1, allele2 = fields[0][4:-1], fields[6][4:-1]
            else:
                allele1 = allele2 = fields[0][4:-1]

            gene = allele1.split("*")[0]
            if gene not in hla_list:
                continue
            if gene not in omixon_alleles:
                omixon_alleles[gene] = set()
            if gene not in allele_count:
                allele_count[gene] = 0
            if allele_count[gene] >= 10:
                continue

            if allele2 == "":
                allele2 = prev_allele2
            assert allele1 != "" and allele2 != ""

            def update_allele(allele):
                if allele == "DRB1*08:01:03":
                    allele = "DRB1*08:01:01"
                elif allele == "DRB1*11:11:02":
                    allele = "DRB1*11:11:01"
                return allele

            allele1, allele2 = update_allele(allele1), update_allele(allele2)
            
            allele_count[gene] += 1
            omixon_alleles[gene].add(allele1)
            omixon_alleles[gene].add(allele2)
            prev_allele1, prev_allele2 = allele1, allele2

            print "%s\t%s\t%s" % (genome, allele1, allele2)

    print >> sys.stderr, "Omixon"
    for gene, alleles in omixon_alleles.items():
        print >> sys.stderr, "\t%s: %d alleles" % (gene, len(alleles))
        for allele in alleles:
            if allele in nuc_alleles[gene]:
                continue
            found = False
            for allele_cmp in nuc_alleles[gene]:
                if allele_cmp.find(allele) != -1:
                    found = True
                    break                    

            if not found:
                print >> sys.stderr, "\t\t%s is missing" % allele

            
