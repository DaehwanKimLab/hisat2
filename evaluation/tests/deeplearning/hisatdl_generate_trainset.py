#!/usr/bin/env python
#
# Copyright 2017, Daehwan Kim <infphilo@gmail.com>
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


import os, sys, subprocess, re
import inspect
from argparse import ArgumentParser, FileType
import hisatgenotype_typing_common as typing_common


"""
"""
def generate_sequences(base_fname,
                       locus_list,
                       verbose):
    base_fullpath_name = base_fname
    if base_dname != "" and not os.path.exists(base_dname):
        os.mkdir(base_dname)
        base_fullpath_name = "%s/%s" % (base_dname, base_fname)

    # Download human genome and HISAT2 index
    HISAT2_fnames = ["grch38",
                     "genome.fa",
                     "genome.fa.fai"]

    if not typing_common.check_files(HISAT2_fnames):
        typing_common.download_genome_and_index()

    chr_dic, _, _ = typing_common.read_genome(open("genome.fa"))

    local_size = 1 << 16
    chr, chr_pos = '1', 1000000
    local_seqs = [chr_dic[chr][chr_pos:chr_pos+local_size], \
                      chr_dic[chr][chr_pos+local_size:chr_pos+local_size*2]]

    seq_len = 20
    seq_dic = {}
    for l in range(len(local_seqs)):
        local_seq = local_seqs[l]
        assert len(local_seq) == local_size
        for p in range(0, local_size - seq_len):
            seq = local_seq[p:p+seq_len]
            if 'N' in seq:
                continue
            if seq not in seq_dic:
                seq_dic[seq] = [l]
            else:
                if seq_dic[seq][-1] != l:
                    seq_dic[seq].append(l)

    seq_file = open("seq.train", 'w')
    label_file = open("label.train", 'w')
    for seq, locals in seq_dic.items():
        assert len(locals) > 0
        def seq_to_binary(seq):
            bseq = []
            for nt in seq:
                if nt == 'A':
                    bseq += ['0', '0']
                elif nt == 'C':
                    bseq += ['0', '1']
                elif nt == 'G':
                    bseq += ['1', '0']
                else:
                    assert nt == 'T'
                    bseq += ['1', '1']
            return bseq
        bseq = seq_to_binary(seq)
        print >> seq_file, ' '.join(bseq)
        if len(locals) == 1:
            if locals[0] == 0:
                label = "1 0"
            else:
                label = "0 1"
        else:
            assert len(locals) == 2
            label = "0.5 0.5"
        print >> label_file, label
    seq_file.close()
    label_file.close()
    
    
    
        
"""
"""
if __name__ == '__main__':
    parser = ArgumentParser(
        description="Generate train sequences")
    parser.add_argument("-b", "--base",
                        dest="base_fname",
                        type=str,
                        default="hla",
                        help="base filename for backbone sequence, variants, and linking info (Default: hla)")
    parser.add_argument("--locus-list",
                        dest="locus_list",
                        type=str,
                        default="",
                        help="A comma-separated list of gene names (default: empty, all genes)")
    parser.add_argument("-v", "--verbose",
                        dest="verbose",
                        action="store_true",
                        help="also print some statistics to stderr")

    args = parser.parse_args()
    if args.locus_list == "":
        locus_list = []
    else:
        locus_list = args.locus_list.split(',')
             
    if args.base_fname.find('/') != -1:
        elems = args.base_fname.split('/')
        base_fname = elems[-1]
        base_dname = '/'.join(elems[:-1])
    else:
        base_fname = args.base_fname
        base_dname = ""
        
    generate_sequences(base_fname,
                       locus_list,
                       args.verbose)

