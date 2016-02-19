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


import sys, os, subprocess
import re
from argparse import ArgumentParser, FileType


def generate_haplotypes(snp_file,
                        haplotype_file,
                        vars,
                        genotypes_list,
                        inter_gap,
                        intra_gap):
    assert len(vars) == len(genotypes_list)
    assert len(vars) > 0
    num_chromosomes = len(genotypes_list[0])

    # Write SNPs
    for var in vars:
        snpID, chr, type, pos, data = var
        print >> snp_file, "%s\t%s\t%s\t%s\t%s" % \
            (snpID, chr, type, pos, data)

    # genotypes_list looks like
    #    Var0: 000001000
    #    Var1: 010000000
    #    Var2: 001100000
    # Get haplotypes from genotypes_list
    haplotypes = set()
    cnv_genotypes = ["" for i in range(num_chromosomes)]
    for genotypes in genotypes_list:
        for i in range(len(genotypes)):
            genotype = genotypes[i]
            cnv_genotypes[i] += genotype
    cnv_genotypes = set(cnv_genotypes)
    for raw_haplotype in cnv_genotypes:
        if '1' not in raw_haplotype:
            continue
        haplotype = ""
        for i in range(len(raw_haplotype)):
            if raw_haplotype[i] == '1':
                if haplotype == "":
                    haplotype = str(i)
                else:
                    haplotype += ("#%d" % i)                    
        assert haplotype != ""            
        haplotypes.add(haplotype)
    # haplotypes look like
    #    '8#10#12#23', '8#12#23', '5#8#12#23#30'

    # Split some haplotypes that include large gaps inside
    def split_haplotypes(haplotypes):
        split_haplotypes = set()
        for haplotype in haplotypes:
            haplotype = haplotype.split('#')
            assert len(haplotype) > 0
            if len(haplotype) == 1:
                split_haplotypes.add(haplotype[0])
                continue
            prev_s, s = 0, 1
            while s < len(haplotype):
                _, _, prev_locus, prev_type, prev_data = vars[int(haplotype[s-1])]
                _, _, locus, type, data = vars[int(haplotype[s])]
                prev_locus, locus = int(prev_locus), int(locus)
                if prev_type == 'D':
                    prev_locus += (int(prev_data) - 1)
                if prev_locus + intra_gap < locus:
                    split_haplotypes.add('#'.join(haplotype[prev_s:s]))
                    prev_s = s
                s += 1
                if s == len(haplotype):
                    split_haplotypes.add('#'.join(haplotype[prev_s:s]))
        return split_haplotypes

    haplotypes2 = split_haplotypes(haplotypes)

    def cmp_haplotype(a, b):
        a = a.split('#')
        _, _, a1_locus, _, _ = vars[int(a[0])]
        _, _, a2_locus, a2_type, a2_data = vars[int(a[-1])]
        a_begin, a_end = int(a1_locus), int(a2_locus)
        if a2_type == 'D':
            a_end += (int(a2_data) - 1)
        b = b.split('#')
        _, _, b1_locus, _, _ = vars[int(b[0])]
        _, _, b2_locus, b2_type, b2_data = vars[int(b[-1])]
        b_begin, b_end = int(b1_locus), int(b2_locus)
        if b2_type == 'D':
            b_end += (int(b2_data) - 1)
        if a_begin != b_begin:
            return a_begin - b_begin
        return a_end - b_end
    haplotypes = sorted(list(haplotypes2), cmp=cmp_haplotype)

    # daehwan - for debugging purposes
    """
    dis = prev_locus - locus
    print "\n[%d, %d]: %d haplotypes" % (i, j, len(haplotypes)), dis
    if len(cur_vars) in range(0, 1000):
        # print "vars:", sorted(list(cur_vars), cmp=cmp_varKey
        print "num:", len(haplotypes)
        for haplotype in haplotypes:
            print haplotype.split('#')
        print "\nnum:", len(haplotypes2)
        for haplotype in haplotypes2:
            print haplotype.split('#')
    """

    # Write haplotypes
    for h_i in range(len(haplotypes)):
        h = haplotypes[h_i].split('#')
        _, chr, h1_locus, _, _ = vars[int(h[0])]
        _, _, h2_locus, h2_type, h2_data = vars[int(h[-1])]
        h_begin, h_end = int(h1_locus), int(h2_locus)
        if h2_type == 'D':
            h_end += (int(h2_data) - 1)
        assert h_begin <= h_end
        varIDs = []
        # for var in h:
        #    varIDs.append(str(var2ID[var]))
            # daehwan - for debugging purposes
            # varIDs.append(var)
        h_new_begin = h_begin
        for h_j in reversed(range(0, h_i)):
            hc = haplotypes[h_j].split('#')
            _, _, hc_begin, hc_type, hc_data = vars[int(hc[-1])]
            hc_begin = int(hc_begin)
            hc_end = hc_begin
            if hc_type == 'D':
                hc_end += (int(hc_data) - 1)
            if hc_end + inter_gap < h_begin:
                break
            if h_new_begin > hc_end:
                h_new_begin = hc_end
        assert h_new_begin <= h_begin
        num_haplotypes = 1
        print >> haplotype_file, "ht%d\t%s\t%d\t%d\t%s" % \
            (num_haplotypes, chr, h_new_begin, h_end, ','.join(varIDs))

        

"""
"""
def main(base_fname,
         species,
         inter_gap,
         intra_gap,
         verbose = False):
    VCF_url_bases = {}
    VCF_url_bases["human"] =  "ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502"
    VCF_fnames = {}
    VCF_fnames["human"] = ["ALL.chr%d.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz" % (i+1) for i in range(16,17)]

    if species not in VCF_url_bases or \
            species not in VCF_fnames:
        print >> sys.stderr, "Error: %s is not supported." % species
        sys.exit(1)

    SNP_file = open("%s.snp" % base_fname, 'w')
    haplotype_file = open("%s.haplotype" % base_fname, 'w')
        
    url_base = VCF_url_bases[species]
    for fname in VCF_fnames[species]:
        if not os.path.exists(fname):
            os.system("wget %s/%s" % (url_base, fname))
        assert os.path.exists(fname)
        gzip_cmd = ["gzip", "-cd", fname]
        gzip_proc = subprocess.Popen(gzip_cmd,
                                     stdout=subprocess.PIPE,
                                     stderr=open("/dev/null", 'w'))
        genomeIDs = []
        vars, genotypes_list = [], []
        curr_right = -1
        for line in gzip_proc.stdout:
            if line.startswith("##"):
                continue

            fields = line.strip().split()
            chr, pos, snpID, ref_allele, alt_alleles, qual, filter, info, format = fields[:9]
            genotypes = fields[9:]

            if line.startswith("#"):
                genomeIDs = genotypes
                continue

            assert len(genotypes) == len(genomeIDs)

            if not snpID.startswith("rs"):
                continue

            pos = int(pos)
            if curr_right + inter_gap < pos and len(vars) > 0:
                assert len(vars) == len(genotypes_list)
                generate_haplotypes(SNP_file,
                                    haplotype_file,
                                    vars,
                                    genotypes_list,
                                    inter_gap,
                                    intra_gap)
                vars, genotypes_list = [], []

            assert ',' not in ref_allele
            alt_alleles = alt_alleles.split(',')
            for a in range(len(alt_alleles)):
                alt_allele = alt_alleles[a]
                ref_allele2, pos2 = ref_allele, pos
                min_len = min(len(ref_allele2), len(alt_allele))
                assert min_len >= 1
                if min_len > 1:
                    ref_allele2 = ref_allele2[min_len - 1:]
                    alt_allele = alt_allele[min_len - 1:]
                    pos2 += (min_len - 1)
                    
                type, data = '', ''
                if len(ref_allele2) == 1 and len(alt_allele) == 1:
                    type = 'S'
                    data = alt_allele
                elif len(ref_allele2) == 1:
                    assert len(alt_allele) > 1
                    type = 'I'
                    data = alt_allele[1:]
                elif len(alt_allele) == 1:
                    assert len(ref_allele2) > 1
                    type = 'D'
                    data = len(ref_allele2) - 1
                else:
                    assert False
            
                cnv_genotypes = []
                for genotype in genotypes:
                    P1, P2 = genotype.split('|')
                    P1, P2 = int(P1), int(P2)
                    if P1 == a + 1:
                        cnv_genotypes.append('1')
                    else:
                        cnv_genotypes.append('0')
                    if P2 == a + 1:
                        cnv_genotypes.append('1')
                    else:
                        cnv_genotypes.append('0')
                # Skip SNPs not present in a given population (e.g. 2,504 genomes in 1000 Genomes Project)
                if '1' not in cnv_genotypes:
                    continue

                vars.append([snpID, chr, pos2, type, data])
                genotypes_list.append(''.join(cnv_genotypes))
                right = pos2
                if type == 'D':
                    right += (int(data) - 1)
                if curr_right < right:
                    curr_right = right

        if len(vars) > 0:
            assert len(vars) == len(genotypes_list)
            generate_haplotypes(SNP_file,
                                haplotype_file,
                                vars,
                                genotypes_list,
                                inter_gap,
                                intra_gap)
            vars, genotypes_list = [], []

    SNP_file.close()
    haplotype_file.close()
    

if __name__ == '__main__':
    parser = ArgumentParser(
        description='Extract haplotypes from a VCF file')
    """
    parser.add_argument('genome_file',
                        nargs='?',
                        type=FileType('r'),
                        help='input genome file')
    parser.add_argument('snp_file',
                        nargs='?',
                        type=FileType('r'),
                        help='input snp file')
    """
    parser.add_argument("-b", "--base",
                        dest="base_fname",
                        type=str,
                        default="",
                        help="base filename for SNPs and haplotypes (default: same as --species)")
    parser.add_argument('--species',
                        dest='species',
                        type=str,
                        default="human",
                        help='species (default: human)')
    parser.add_argument("--inter-gap",
                        dest="inter_gap",
                        type=int,
                        default=30,
                        help="Maximum distance for variants to be in the same haplotype")
    parser.add_argument("--intra-gap",
                        dest="intra_gap",
                        type=int,
                        default=50,
                        help="Break a haplotype into several haplotypes")
    parser.add_argument('-v', '--verbose',
                        dest='verbose',
                        action='store_true',
                        help='also print some statistics to stderr')

    args = parser.parse_args()
    if args.base_fname == "":
        args.base_fname = args.species

    main(args.base_fname,
         args.species,
         args.inter_gap,
         args.intra_gap,
         args.verbose)
