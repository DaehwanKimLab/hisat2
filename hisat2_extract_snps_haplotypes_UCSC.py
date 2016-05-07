#!/usr/bin/env python

#
# Copyright 2015, Daehwan Kim <infphilo@gmail.com>
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


import sys, subprocess
import re
from argparse import ArgumentParser, FileType


"""
"""
def reverse_complement(seq):
    result = ""
    for nt in seq:
        base = nt
        if nt == 'A':
            base = 'T'
        elif nt == 'a':
            base = 't'
        elif nt == 'C':
            base = 'G'
        elif nt == 'c':
            base = 'g'
        elif nt == 'G':
            base = 'C'
        elif nt == 'g':
            base = 'c'
        elif nt == 'T':
            base = 'A'
        elif nt == 't':
            base = 'a'

        result = base + result

    return result


"""
"""
def read_genome(genome_file):
    chr_dic = {}
    chr_name, sequence = "", ""
    for line in genome_file:
        if line.startswith(">"):
            if chr_name and sequence:
                chr_dic[chr_name] = sequence
            chr_name = line.strip().split()[0][1:]
            sequence = ""
        else:
            sequence += line.strip()
    if chr_name and sequence:
        chr_dic[chr_name] = sequence
    return chr_dic


"""
Compare two variants [chr, pos, type, data, dic]
"""
def compare_vars(a, b):
    a_chr, a_pos, a_type, a_data = a[:4]
    b_chr, b_pos, b_type, b_data = b[:4]

    # daehwan - for debugging purposes
    if a_chr != b_chr:
        print a
        print b
    
    assert a_chr == b_chr
    if a_pos != b_pos:
        return a_pos - b_pos
    if a_type != b_type:
         if a_type == 'I':
             return -1
         elif b_type == 'I':
             return 1
         if a_type == 'S':
             return -1
         else:
             return 1
    if a_data < b_data:
        return -1
    elif a_data > b_data:
        return 1
    else:
        return 0


"""
"""
def compatible_vars(a, b):
    a_chr, a_pos, a_type, a_data = a[:4]
    b_chr, b_pos, b_type, b_data = b[:4]
    assert a_chr == b_chr
    assert a_pos <= b_pos
    if a_pos == b_pos:
        return False
    if a_type == 'D':
        if b_pos <= a_pos + a_data:
            return False
    return True


"""
"""
def generate_haplotypes(snp_file,
                        haplotype_file,
                        vars,
                        inter_gap,
                        intra_gap,
                        num_haplotypes):
    assert len(vars) > 0

    # Sort variants and remove redundant variants
    vars = sorted(vars, cmp=compare_vars)
    tmp_vars = []
    v = 0
    while v < len(vars):
        var = vars[v]
        for v2 in range(v + 1, len(vars)):
            var2 = vars[v2]
            if compare_vars(var, var2) == 0:
                v += 1
            else:
                assert compare_vars(var, var2) < 0
                break
        tmp_vars.append(var)
        v += 1
    vars = tmp_vars

    # Create new variant ID for variants with the same ID
    # e.g. same two variant ID, rs60160543, are split into rs60160543.0 and rs60160543.1
    vars_count = {}
    for var in vars:
        id = var[4]["id"]
        if id not in vars_count:
            vars_count[id] = 0
        vars_count[id] += 1
    vars_duplicate = set()
    for id, count in vars_count.items():
        if count <= 1:
            continue
        vars_duplicate.add(id)
    vars_count = {}
    for var in vars:
        id = var[4]["id"]
        if id not in vars_count:
            vars_count[id] = 0
        else:
            vars_count[id] += 1
        if id not in vars_duplicate:
            var[4]["id2"] = id
        else:
            var[4]["id2"] = "%s.%d" % (id, vars_count[id])

    # variant compatibility
    vars_cmpt = [-1 for i in range(len(vars))]
    for v in range(len(vars)):
        var_chr, var_pos, var_type, var_data = vars[v][:4]
        if var_type == 'D':
            var_pos += (var_data - 1)
        for v2 in range(v + 1, len(vars)):
            if vars_cmpt[v2] >= 0:
                continue
            var2_chr, var2_pos = vars[v2][:2]
            if var_chr != var2_chr:
                break
            if var_pos + inter_gap < var2_pos:
                break
            vars_cmpt[v2] = v

    # Assign genotypes for those missing genotypes
    genotypes_list = []
    for v in range(len(vars)):
        var = vars[v]
        var_dic = var[4]
        freq = var_dic["freq"]
        used = [False for i in range(100)]
        if vars_cmpt[v] >= 0:
            v2 = v - 1
            while v2 >= vars_cmpt[v]:
                var2 = vars[v2]
                if not compatible_vars(var2, var) or \
                        freq >= 0.1:
                    var2_dic = var2[4]
                    assert "genotype" in var2_dic
                    genotype_num = var2_dic["genotype"]
                    used[genotype_num] = True
                v2 -= 1

        assert False in used
        for i in range(len(used)):
            if not used[i]:                
                var_dic["genotype"] = i
                break
        genotypes_list.append(var_dic["genotype"])

    # Write SNPs into a file (.snp)
    for var in vars:
        chr, pos, type, data, var_dic = var
        varID = var_dic["id2"]
        if type == 'S':
            type = "single"
        elif type == 'D':
            type = "deletion"
        else:
            assert type == 'I'
            type = "insertion"
        print >> snp_file, "%s\t%s\t%s\t%s\t%s" % \
            (varID, type, chr, pos, data)

    # genotypes_list looks like
    #    Var0: 0
    #    Var1: 0
    #    Var2: 1
    #    Var3: 2
    # Get haplotypes from genotypes_list
        
    max_genotype_num = max(genotypes_list)
    haplotypes = ["" for i in range(max_genotype_num + 1)]
    for i in range(len(genotypes_list)):
        num = genotypes_list[i]
        if haplotypes[num] == "":
            haplotypes[num] = str(i)
        else:
            haplotypes[num] += ("#%d" % i)
    haplotypes = set(haplotypes)
    
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
                _, prev_locus, prev_type, prev_data, _ = vars[int(haplotype[s-1])]
                _, locus, type, data, _ = vars[int(haplotype[s])]
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
        _, a1_locus, _, _, _ = vars[int(a[0])]
        _, a2_locus, a2_type, a2_data, _ = vars[int(a[-1])]
        a_begin, a_end = int(a1_locus), int(a2_locus)
        if a2_type == 'D':
            a_end += (int(a2_data) - 1)
        b = b.split('#')
        _, b1_locus, _, _, _ = vars[int(b[0])]
        _, b2_locus, b2_type, b2_data, _ = vars[int(b[-1])]
        b_begin, b_end = int(b1_locus), int(b2_locus)
        if b2_type == 'D':
            b_end += (int(b2_data) - 1)
        if a_begin != b_begin:
            return a_begin - b_begin
        return a_end - b_end
    
    haplotypes = sorted(list(haplotypes2), cmp=cmp_haplotype)

    # Write haplotypes
    for h_i in range(len(haplotypes)):
        h = haplotypes[h_i].split('#')
        chr, h1_locus, _, _, _ = vars[int(h[0])]
        _, h2_locus, h2_type, h2_data, _ = vars[int(h[-1])]
        h_begin, h_end = int(h1_locus), int(h2_locus)
        if h2_type == 'D':
            h_end += (int(h2_data) - 1)
        assert h_begin <= h_end
        h_new_begin = h_begin
        for h_j in reversed(range(0, h_i)):
            hc = haplotypes[h_j].split('#')
            _, hc_begin, hc_type, hc_data, _ = vars[int(hc[-1])]
            hc_begin = int(hc_begin)
            hc_end = hc_begin
            if hc_type == 'D':
                hc_end += (int(hc_data) - 1)
            if hc_end + inter_gap < h_begin:
                break
            if h_new_begin > hc_end:
                h_new_begin = hc_end
        assert h_new_begin <= h_begin
        h_add = []
        for id in h:
            var_dic = vars[int(id)][4]
            h_add.append(var_dic["id2"])
        print >> haplotype_file, "ht%d\t%s\t%d\t%d\t%s" % \
            (num_haplotypes, chr, h_new_begin, h_end, ','.join(h_add))
        num_haplotypes += 1

    return num_haplotypes


"""
"""
def main(genome_file,
         snp_fname,
         base_fname,
         inter_gap,
         intra_gap,
         verbose,
         testset):
    # load genomic sequences
    chr_dic = read_genome(genome_file)

    if testset:
        ref_testset_file = open(base_fname + ".ref.testset.fa", "w")
        alt_testset_file = open(base_fname + ".alt.testset.fa", "w")

    snp_out_file = open(base_fname + ".snp", 'w')
    haplotype_out_file = open(base_fname + ".haplotype", 'w')

    # load SNPs
    snp_list = []
    prev_chr, curr_right = "", -1
    num_haplotypes = 0
    if snp_fname.endswith(".gz"):
        snp_cmd = ["gzip", "-cd", snp_fname]
    else:
        snp_cmd = ["cat", snp_fname]
    snp_proc = subprocess.Popen(snp_cmd,
                                stdout=subprocess.PIPE,
                                stderr=open("/dev/null", 'w'))
    ids_seen = set()
    for line in snp_proc.stdout:
        if not line or line.startswith('#'):
            continue

        line = line.strip()
        try:
            fields = line.split('\t')
            """
            id, chr, start, end, rs_id, score, strand, refNCBI, refUCSC, observed, molType, classType, valid, \
                avHet, avHetSE, func, locType, weight, exceptions, submitterCount, submitters, \
                alleleFreqCount, alleles, alleleNs, alleleFreqs, bitfields = fields
            """
            id, chr, start, end, rs_id, score, strand, refNCBI, refUCSC, observed, molType, classType = fields[:12]
            alleleFreqs = fields[-2].split(',')[:-1]
            if len(alleleFreqs) > 0:
                try:
                    float(alleleFreqs[0])
                except ValueError:
                    alleleFreqs = []
        except ValueError:
            continue

        start, end = int(start), int(end)
        score = int(score)

        if molType != "genomic":
            continue

        if classType not in ["single", "deletion", "insertion"]:
            continue

        if classType == "single":
            if start + 1 != end:
                continue
        elif classType == "deletion":
            assert start < end
        else:
            assert classType == "insertion"
            if start != end:
                continue

        if chr not in chr_dic:
            continue
        chr_seq = chr_dic[chr]
        chr_len = len(chr_seq)

        if start >= len(chr_seq):
            continue

        if rs_id in ids_seen:
            continue
        ids_seen.add(rs_id)

        if (prev_chr != chr or curr_right + inter_gap < start) and \
                len(snp_list) > 0:
            num_haplotypes = generate_haplotypes(snp_out_file,
                                                 haplotype_out_file,
                                                 snp_list,
                                                 inter_gap,
                                                 intra_gap,
                                                 num_haplotypes)
            snp_list = []

        observed = observed.upper()
        allele_list = observed.split("/")
        if len(alleleFreqs) == 0:
            alleleFreqs = [0.0 for i in range(len(allele_list))]
        
        # Reverse complement alleles if strand is negative
        if strand == "-":
            tmp_allele_list = []
            for allele in allele_list:
                tmp_allele_list.append(reverse_complement(allele))
            allele_list = tmp_allele_list
            
        if classType == "single":
            allele_count = min(len(allele_list), len(alleleFreqs))
            ref_base = chr_seq[start].upper()
            if ref_base not in allele_list:
                continue
            for a in range(allele_count):
                allele = allele_list[a]
                freq = float(alleleFreqs[a])
                if allele not in "ACGT" or len(allele) != 1:
                    continue
                if allele == ref_base:
                    continue
                snp_list.append([chr, start, 'S', allele, {"id":rs_id, "freq":freq}])

                if testset:
                    ref_seq = chr_seq[start-50:start+50]
                    alt_seq = chr_seq[start-50:start] + allele + chr_seq[start+1:start+50]
                    print >> ref_testset_file, ">%s_single_%d" % (rs_id, start - 50)
                    print >> ref_testset_file, ref_seq
                    print >> alt_testset_file, ">%s_single_%d_%s" % (rs_id, start - 50, ref_seq)
                    print >> alt_testset_file, alt_seq
                
        elif classType == "deletion":
            if start > 0:
                prev_base = chr_seq[start-1].upper()
                if prev_base not in "ACGT":
                    continue

            if len(allele_list) != 2 or \
                    len(allele_list) != len(alleleFreqs):
                continue
                
            freq = 0.0
            if allele_list[0] == "-":
                freq = float(alleleFreqs[1])
            else:
                assert allele_list[1] == "-"
                freq = float(alleleFreqs[0])
            
            delLen = end - start
            snp_list.append([chr, start, 'D', delLen, {"id":rs_id, "freq":freq}])
            if testset and delLen > 0 and delLen <= 10:
                ref_seq = chr_seq[start-50:start+50]
                alt_seq = chr_seq[start-50:start] + chr_seq[start+delLen:start+50+delLen]
                print >> ref_testset_file, ">%s_deletion_%d" % (rs_id, start - 50)
                print >> ref_testset_file, ref_seq
                print >> alt_testset_file, ">%s_deletion_%d_%s" % (rs_id, start - 50, ref_seq)
                print >> alt_testset_file, alt_seq
        else:
            assert classType == "insertion"
            if start > 0:
                prev_base = chr_seq[start-1].upper()
                if prev_base not in "ACGT":
                    continue
            allele_count = min(len(allele_list), len(alleleFreqs))
            for a in range(allele_count):
                allele = allele_list[a]
                freq = float(alleleFreqs[a])
                if allele == "-" or len(allele) <= 0:
                    continue
                if re.match('^[ACGT]+$', allele):
                    snp_list.append([chr, start, 'I', allele, {"id":rs_id, "freq":freq}])
                    insLen = len(allele)
                    if testset and insLen > 0 and insLen <= 10:
                        ref_seq = chr_seq[start-50:start+50]
                        alt_seq = chr_seq[start-50:start] + allele + chr_seq[start:start+50-insLen]
                        print >> ref_testset_file, ">%s_insertion_%d" % (rs_id, start - 50)
                        print >> ref_testset_file, ref_seq
                        print >> alt_testset_file, ">%s_insertion_%d_%s" % (rs_id, start - 50, ref_seq)
                        print >> alt_testset_file, alt_seq

        if curr_right < end:
            curr_right = end

        if prev_chr != chr:
            curr_right = end
        prev_chr = chr

    if testset:
        ref_testset_file.close()
        alt_testset_file.close()

    if len(snp_list) > 0:
        generate_haplotypes(snp_out_file,
                            haplotype_out_file,
                            snp_list,
                            inter_gap,
                            intra_gap,
                            num_haplotypes)
        snp_list = []
    
    snp_out_file.close()
    haplotype_out_file.close()
                       
                     

if __name__ == '__main__':
    parser = ArgumentParser(
        description='Extract SNPs and haplotypes from a SNP file downloaded from UCSC (e.g. http://hgdownload.soe.ucsc.edu/goldenPath/hg38/database/snp144.txt.gz)')
    parser.add_argument('genome_file',
                        nargs='?',
                        type=FileType('r'),
                        help='input genome file (e.g. genome.fa)')
    parser.add_argument('snp_fname',
                        nargs='?',
                        type=str,
                        help='input snp file downloaded from UCSC (plain text or gzipped file is accepted: snp144Common.txt or snp144Common.txt.gz)')
    parser.add_argument("base_fname",
                        nargs='?',
                        type=str,
                        help="base filename for SNPs and haplotypes")
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
    parser.add_argument('--testset',
                        dest='testset',
                        action='store_true',
                        help='print test reads')

    args = parser.parse_args()
    if not args.genome_file or \
            not args.snp_fname or \
            not args.base_fname:
        parser.print_help()
        exit(1)
    main(args.genome_file,
         args.snp_fname,
         args.base_fname,
         args.inter_gap,
         args.intra_gap,
         args.verbose,
         args.testset)
