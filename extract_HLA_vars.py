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


import sys, re
from collections import defaultdict as dd, Counter
from argparse import ArgumentParser, FileType



"""
"""
def extract_HLA_vars(HLA_MSA_file, base_fname, verbose = False):
    # Samples of HLA_MSA_file are found in
    #    ftp://ftp.ebi.ac.uk/pub/databases/ipd/imgt/hla/msf/

    HLA_names = {} # HLA allele names to numeric IDs
    HLA_seqs = []  # HLA multiple alignment sequences
    for line in HLA_MSA_file:
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

    print >> sys.stderr, "Number of HLA genes is %d." % (len(HLA_names))

    Vars = {}
    backbone_name = "A*01:01:01:01"
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
            if bc == cc:
                if bc != ".":
                    if insertion:
                        insertVar(insertion, 'I')
                        insertion = []
                    elif deletion:
                        insertVar(deletion, 'D')
                        deletion = []
            else:
                if bc != '.' and cc != '.':
                    mismatch = [s - ndots, cc]
                    insertVar(mismatch, 'M')
                else:
                    if bc == '.':
                        if insertion:
                            insertion[1] += cc
                        else:
                            insertion = [s - ndots, cc]
                    else:
                        assert cc == '.'
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
            if a_type == 'M':
                return -1
            elif b_type == 'M':
                return 1
            elif a_type == 'I':
                return -1
            else:
                assert b_type == 'I'
                return 1
        assert a_data != b_data
        if len(a_data) != len(b_data):
            return len(a_data) - len(b_data)            
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

    # sanity check -
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

        # daehwan - for debugging purposes
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
    backbone_file = open(base_fname + "_backbone.fa", 'w')
    print >> backbone_file, ">%s" % (backbone_name)
    backbone_seq_ = backbone_seq.replace('.', '')
    for s in range(0, len(backbone_seq_), 60):
        print >> backbone_file, backbone_seq_[s:s+60]
    backbone_file.close()
    
    # Write
    #       (1) variants w.r.t the backbone sequences into a SNP file
    #       (2) pairs of a variant and the corresponding HLA allels into a LINK file    
    keys = sorted(Vars.keys(), cmp=cmp_varKey)
    var_file = open(base_fname + ".snp", 'w')
    link_file = open(base_fname + ".link", 'w')
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
        print >> var_file, "hv%d\t%s\t%s\t%d\t%s" % \
            (k+1, type_str, backbone_name, locus, data)

        names = sorted(Vars[keys[k]])
        print >> link_file, "hv%d\t%s" % (k+1, ' '.join(names))
    var_file.close()
    link_file.close()

    # Write all the sequences with dots removed into a file
    input_file = open(base_fname + "_sequences.fa", 'w')
    for name, ID in HLA_names.items():
        print >> input_file, ">%s" % (name)
        assert ID < len(HLA_seqs)
        seq = HLA_seqs[ID].replace('.', '')
        for s in range(0, len(seq), 60):
            print >> input_file, seq[s:s+60]
    input_file.close()

    
    
        
"""
"""
if __name__ == '__main__':
    parser = ArgumentParser(
        description='Extract HLA variants from HLA multiple sequence alignments')
    parser.add_argument('HLA_MSA_file',
        nargs='?',
        type=FileType('r'),
        help='input snp file')
    parser.add_argument('-b', '--base',
        dest='base_fname',
        type=str,
        default="hla",
        help='base filename for backbone HLA sequence, HLA variants, and HLA linking info.')
    parser.add_argument('-v', '--verbose',
        dest='verbose',
        action='store_true',
        help='also print some statistics to stderr')

    args = parser.parse_args()
    if not args.HLA_MSA_file:
        parser.print_help()
        exit(1)
    extract_HLA_vars(args.HLA_MSA_file, args.base_fname, args.verbose)
