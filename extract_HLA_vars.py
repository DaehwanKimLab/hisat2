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


import os, sys, re
from collections import defaultdict as dd, Counter
from argparse import ArgumentParser, FileType



"""
"""
def extract_HLA_vars(base_fname, gap, split, verbose = False):
    # Samples of HLA_MSA_file are found in
    #    ftp://ftp.ebi.ac.uk/pub/databases/ipd/imgt/hla/msf/
    HLA_genes = {
        "A" : "A*01:01:01:01",
        "B" : "B*07:02:01",
        "C" : "C*01:02:01",       # ref: C*01:02:01 
        "DQA" : "",
        "DQB" : "",
        "DRB1" : "DRB1*01:01:01"  # ref: DRB1*15:01:01:01
        }

    # Write the backbone sequences into a fasta file
    backbone_file = open(base_fname + "_backbone.fa", 'w')        
    # variants w.r.t the backbone sequences into a SNP file
    var_file = open(base_fname + ".snp", 'w')
    # haplotypes
    haplotype_file = open(base_fname + ".haplotype", 'w')
    # pairs of a variant and the corresponding HLA allels into a LINK file    
    link_file = open(base_fname + ".link", 'w')
    # Write all the sequences with dots removed into a file
    input_file = open(base_fname + "_sequences.fa", 'w')
    num_vars, num_haplotypes = 0, 0
    for HLA_gene, HLA_ref_gene in HLA_genes.items():
        if HLA_ref_gene == "":
            continue
        HLA_MSA_fname = "IMGTHLA/msf/%s_gen.msf" % HLA_gene
        if not os.path.exists(HLA_MSA_fname):
            print >> sys.stderr, "Warning: %s does not exist" % HLA_MSA_fname
            continue
        
        HLA_names = {} # HLA allele names to numeric IDs
        HLA_seqs = []  # HLA multiple alignment sequences
        for line in open(HLA_MSA_fname):
            line = line.strip()
            if not line or \
                    not line[0].isalnum():
                continue

            if line.startswith("MSF"):
                continue

            if line.startswith("Name"):
                try:
                    name = line.split('\t')[0]
                    name = name.split()[1]
                except ValueError:
                    continue

                if name in HLA_names:
                    print >> sys.stderr, "Warning: %s is found more than once in Names" % (name)
                    continue

                HLA_names[name] = len(HLA_names)
            else:
                if len(HLA_seqs) == 0:
                    HLA_seqs = ["" for i in range(len(HLA_names))]
                try:
                    name, five1, five2, five3, five4, five5 = line.split()
                except ValueError:
                    continue

                if name not in HLA_names:
                    print >> sys.stderr, "Warning: %s is not present in Names" % (name)
                    continue

                id = HLA_names[name]
                HLA_seqs[id] += (five1 + five2 + five3 + five4 + five5)

        # sanity check -
        #    Assert the lengths of the input MSF are the same
        assert len(HLA_seqs) > 0
        seq_len = len(HLA_seqs[0])
        for i in range(1, len(HLA_seqs)):
            assert seq_len == len(HLA_seqs[i])

        print >> sys.stderr, "%s: number of HLA genes is %d." % (HLA_gene, len(HLA_names))

        Vars = {}
        backbone_name = HLA_ref_gene
        backbone_id = HLA_names[backbone_name]
        backbone_seq = HLA_seqs[backbone_id]
        for cmp_name, id in HLA_names.items():
            if cmp_name == backbone_name:
                continue
            assert id < len(HLA_seqs)
            cmp_seq = HLA_seqs[id]

            """
            for s in range(0, seq_len, 100):
                print s, backbone_seq[s:s+100]
                print s, cmp_seq[s:s+100]
            """

            def insertVar(indel, type):
                if type in "MI":
                    varKey = "%d-%s-%s" % (indel[0], type, indel[1])
                else:
                    varKey = "%d-%s-%d" % (indel[0], type, indel[1])
                if varKey not in Vars:
                    Vars[varKey] = [cmp_name]
                else:
                    Vars[varKey].append(cmp_name)

            insertion, deletion = [], []
            ndots = 0
            for s in range(seq_len):
                assert not (insertion and deletion)
                bc = backbone_seq[s]
                cc = cmp_seq[s]
                if bc != '.' and cc != '.':
                    if insertion:
                        insertVar(insertion, 'I')
                        insertion = []
                    elif deletion:
                        insertVar(deletion, 'D')
                        deletion = []
                    if bc != cc:
                        mismatch = [s - ndots, cc]
                        insertVar(mismatch, 'M')
                elif bc == '.' and cc != '.':
                    if deletion:
                        insertVar(deletion, 'D')
                        deletion = []
                    if insertion:
                        insertion[1] += cc
                    else:
                        insertion = [s - ndots, cc]
                elif bc != '.' and cc == '.':
                    if insertion:
                        insertVar(insertion, 'I')
                        insertion = []
                    if deletion:
                        deletion[1] += 1
                    else:
                        deletion = [s - ndots, 1]

                if bc == '.':
                    ndots += 1


                """
                if backbone_seq[s] != cmp_seq[s]:
                    print "%s is different %s at %d: %s vs. %s" % \
                        (backbone_name, cmp_name, s+1, backbone_seq[s], cmp_seq[s])
                """

            if insertion:
                insertVar(insertion, 'I')
            elif deletion:
                insertVar(deletion, 'D')


        print >> sys.stderr, "Number of variants is %d." % (len(Vars.keys()))

        # Compare variants
        def cmp_varKey(a, b):
            a_locus, a_type, a_data = a.split('-')
            b_locus, b_type, b_data = b.split('-')
            a_locus, b_locus = int(a_locus), int(b_locus)
            if a_locus != b_locus:
                return a_locus - b_locus
            if a_type != b_type:
                if a_type == 'I':
                    return -1
                elif b_type == 'I':
                    return 1
                elif a_type == 'M':
                    return -1
                else:
                    assert b_type == 'M'
                    return 1
            assert a_data != b_data
            if a_type in "MI":
                if a_data < b_data:
                    return -1
                else:
                    return 1
            else:
                assert a_type == 'D'
                return int(a_data) - int(b_data)            

        HLA_Vars = {}
        for key, names in Vars.items():
            for name in names:
                if not name in HLA_Vars:
                    HLA_Vars[name] = [key]
                else:
                    HLA_Vars[name].append(key)
        for name, vars in HLA_Vars.items():
            HLA_Vars[name] = sorted(vars, cmp=cmp_varKey)

        # Sanity check -
        #    (1) Reconstruct the other sequences from the backbone sequence and variants and
        #    (2) Confirm these constructed sequences are the same as those input sequences.
        for cmp_name, id in HLA_names.items():
            if cmp_name == backbone_name:
                continue

            constr_seq = backbone_seq.replace('.', '')
            constr_seq = list(constr_seq)
            locus_diff = 0
            for var in HLA_Vars[cmp_name]:
                try:
                    locus, type, data = var.split('-')
                    locus = int(locus)
                except ValueError:
                    continue

                if type == 'M':
                    assert len(data) == 1
                    constr_seq[locus + locus_diff] = data[0]
                elif type == 'I':
                    assert locus + locus_diff >= 0
                    assert locus + locus_diff <= len(constr_seq)
                    constr_seq = constr_seq[:locus + locus_diff] + list(data) + constr_seq[locus + locus_diff:]
                    locus_diff += len(data)
                else:
                    assert type == 'D'
                    assert locus + locus_diff + len(data) <= len(constr_seq)
                    assert locus + locus_diff >= 0
                    del_len = int(data)
                    constr_seq = constr_seq[:locus + locus_diff] + constr_seq[locus + locus_diff + del_len:]
                    locus_diff -= del_len

            constr_seq = "".join(constr_seq)
            assert id < len(HLA_seqs)
            cmp_seq = HLA_seqs[id].replace('.', '')
            if len(constr_seq) != len(cmp_seq):
                print >> sys.stderr, "Error: reconstruction fails (%s)! Lengths different: %d vs. %d" % \
                    (cmp_name, len(constr_seq), len(cmp_seq))
                assert False

            # Sanity check
            for s in range(len(constr_seq)):
                if constr_seq[s] != cmp_seq[s]:
                    print >> sys.stderr, "Differ at %d: %s vs. %s (reconstruction vs. original)" % \
                        (s, constr_seq[s], cmp_seq[s])
                    print "%s:%s vs. %s:%s" % \
                        (constr_seq[s-10:s], constr_seq[s:s+10], cmp_seq[s-10:s], cmp_seq[s:s+10])

            if constr_seq != cmp_seq.replace('.', ''):
                print >> sys.stderr, "Error: reconstruction fails for %s" % (cmp_name)
                assert False

        # Write the backbone sequences into a fasta file
        print >> backbone_file, ">%s" % (backbone_name)
        backbone_seq_ = backbone_seq.replace('.', '')
        for s in range(0, len(backbone_seq_), 60):
            print >> backbone_file, backbone_seq_[s:s+60]

        # Write
        #       (1) variants w.r.t the backbone sequences into a SNP file
        #       (2) pairs of a variant and the corresponding HLA allels into a LINK file    
        keys = sorted(Vars.keys(), cmp=cmp_varKey)
        var2ID = {}
        for k in range(len(keys)):
            locus, type, data = keys[k].split('-')
            locus = int(locus)
            if type == 'M':
                type_str = "single"
            elif type == 'I':
                type_str = "insertion"
            else:
                assert type == 'D'
                type_str = "deletion"

            varID = "hv%d" % (num_vars)
            print >> var_file, "%s\t%s\t%s\t%d\t%s" % \
                (varID, type_str, backbone_name, locus, data)
            names = sorted(Vars[keys[k]])
            print >> link_file, "%s\t%s" % (varID, ' '.join(names))
            var2ID[keys[k]] = num_vars
            num_vars += 1

        add_seq_len = 0
        # Write haplotypes
        i = 0
        while i < len(keys):
            key_i = keys[i]
            locus, type, data = key_i.split('-')
            locus = int(locus)
            if type == 'D':
                locus += (int(data)- 1)
            prev_locus = locus
            j = i + 1
            while j < len(keys):
                key_j = keys[j]
                locus2, type2, data2 = key_j.split('-')
                locus2 = int(locus2)
                if prev_locus + gap < locus2:
                    break
                prev_locus = locus2
                if type == 'D':
                    prev_locus += (int(data)- 1)                
                j += 1

            alleles = set()
            for k in range(i, j):
                key_k = keys[k]
                add_alleles = set(Vars[key_k])
                alleles |= add_alleles

            haplotypes = set()
            cur_vars = set(keys[i:j])
            for allele in alleles:
                allele_vars = set(HLA_Vars[allele])
                allele_cur_vars = '#'.join(sorted(list(cur_vars & allele_vars), cmp=cmp_varKey))
                haplotypes.add(allele_cur_vars)

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
                        prev_locus, prev_type, prev_data = haplotype[s-1].split('-')
                        locus, type, data = haplotype[s].split('-')
                        prev_locus, locus = int(prev_locus), int(locus)
                        if prev_type == 'D':
                            prev_locus += (int(prev_data) - 1)
                        if prev_locus + split < locus:
                            split_haplotypes.add('#'.join(haplotype[prev_s:s]))
                            prev_s = s
                        s += 1
                        if s == len(haplotype):
                            split_haplotypes.add('#'.join(haplotype[prev_s:s]))
                return split_haplotypes

            haplotypes2 = split_haplotypes(haplotypes)

            def cmp_haplotype(a, b):
                a = a.split('#')
                a1_locus, _, _ = a[0].split('-')
                a2_locus, a2_type, a2_data = a[-1].split('-')
                a_begin, a_end = int(a1_locus), int(a2_locus)
                if a2_type == 'D':
                    a_end += (int(a2_data) - 1)
                b = b.split('#')
                b1_locus, _, _ = b[0].split('-')
                b2_locus, b2_type, b2_data = b[-1].split('-')
                b_begin, b_end = int(b1_locus), int(b2_locus)
                if b2_type == 'D':
                    b_end += (int(b2_data) - 1)
                if a_begin != b_begin:
                    return a_begin - b_begin
                return a_end - b_end

            haplotypes = sorted(list(haplotypes), cmp=cmp_haplotype)
            haplotypes2 = sorted(list(haplotypes2), cmp=cmp_haplotype)
            
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
            sanity_vars = set()
            for h_i in range(len(haplotypes2)):
                h = haplotypes2[h_i].split('#')
                h1_locus, _, _ = h[0].split('-')
                h2_locus, h2_type, h2_data = h[-1].split('-')
                h_begin, h_end = int(h1_locus), int(h2_locus)
                if h2_type == 'D':
                    h_end += (int(h2_data) - 1)
                assert h_begin <= h_end
                varIDs = []
                for var in h:
                    varIDs.append(str(var2ID[var]))
                    # daehwan - for debugging purposes
                    # varIDs.append(var)
                    sanity_vars.add(var2ID[var])
                h_new_begin = h_begin
                for h_j in reversed(range(0, h_i)):
                    hc = haplotypes2[h_j].split('#')
                    hc_begin, hc_type, hc_data = hc[-1].split('-')
                    hc_begin = int(hc_begin)
                    hc_end = hc_begin
                    if hc_type == 'D':
                        hc_end += (int(hc_data) - 1)
                    if hc_end + gap < h_begin:
                        break
                    if h_new_begin > hc_end:
                        h_new_begin = hc_end
                assert h_new_begin <= h_begin
                print >> haplotype_file, "ht%d\t%s\t%d\t%d\t%s" % \
                    (num_haplotypes, backbone_name, h_new_begin, h_end, ','.join(varIDs))
                num_haplotypes += 1
                add_seq_len += (h_end - h_new_begin + 1)
            assert len(sanity_vars) == len(cur_vars)
                    
            i = j

        print >> sys.stderr, "Length of additional sequences for haplotypes:", add_seq_len
                    
        # Write all the sequences with dots removed into a file
        for name, ID in HLA_names.items():
            print >> input_file, ">%s" % (name)
            assert ID < len(HLA_seqs)
            seq = HLA_seqs[ID].replace('.', '')
            for s in range(0, len(seq), 60):
                print >> input_file, seq[s:s+60]

    backbone_file.close()
    var_file.close()
    haplotype_file.close()
    link_file.close()
    input_file.close()
    
        
"""
"""
if __name__ == '__main__':
    parser = ArgumentParser(
        description="Extract HLA variants from HLA multiple sequence alignments")
    parser.add_argument("-b", "--base",
                        dest="base_fname",
                        type=str,
                        default="hla",
                        help="base filename for backbone HLA sequence, HLA variants, and HLA linking info.")
    parser.add_argument("-g", "--gap",
                        dest="gap",
                        type=int,
                        default=30,
                        help="Maximum distance for variants to be in the same haplotype.")
    parser.add_argument("-s", "--split",
                        dest="split",
                        type=int,
                        default=50,
                        help="Break a haplotype into several haplotypes.")
    parser.add_argument("-v", "--verbose",
                        dest="verbose",
                        action="store_true",
                        help="also print some statistics to stderr")

    args = parser.parse_args()
    if args.gap > args.split:
        print >> sys.stderr, "Error: -g/--gap (%d) must be smaller than -s/--split (%d)" % (args.gap, args.split)
        sys.exit(1)
        
    extract_HLA_vars(args.base_fname, args.gap, args.split, args.verbose)
