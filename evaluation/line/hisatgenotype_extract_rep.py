#!/usr/bin/env python
#
# Copyright 2017, Daehwan Kim <infphilo@gmail.com> and Hyunmin Kim <human.gim@gmail.com>
#
# This file is part of HISAT-rep.
#
# HISAT-rep is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# HISAT-bisulfite is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with HISAT-rep.  If not, see <http://www.gnu.org/licenses/>.
#


import os, sys, subprocess, re
import heapq
import inspect
from argparse import ArgumentParser, FileType
import hisatgenotype_typing_common as typing_common

"""
Extract human specific LINEs from RepBase database such as RepBase22.08.fasta.tar.gz
"""
def extract_L1HS(rep_fname,
                 l1hs_fname):
    l1hs_file = open(l1hs_fname, 'w')
    for line in open(rep_fname):
        line = line.rstrip()
        if line.startswith('>'):
            if line.startswith(">L1") and line.find("Homo sapiens") != -1:
                extract = True
            else:
                extract = False
        if extract:
            print >> l1hs_file, line
    l1hs_file.close()


"""
"""
def extract_rep(base_fname,
                verbose):
    rep_fname = "RepBase22.08"
    if not os.path.exists(rep_fname + ".fasta") and not os.path.exists(rep_fname + ".fasta.tar.gz"):
        print >> sys.stderr, "Error: please have %s.fasta.tar.gz copied in the current directory." % rep_fname
        sys.exit(1)

    if not os.path.exists(rep_fname + ".fasta"):
        os.system("tar xvzf %s.fasta.tar.gz" % rep_fname)
        assert os.path.exists(rep_fname + ".fasta")
    
    l1hs_fname = "%s.l1hs.fasta" % rep_fname
    if not os.path.exists(l1hs_fname):
        extract_L1HS("%s.fasta/humrep.ref" % rep_fname,
                     l1hs_fname)
        assert os.path.exists(l1hs_fname)

    clustalw_fname = "%s.l1hs.aln" % rep_fname
    # Run ClustalW2 to perform multiple sequencing alignment
    if not os.path.exists(clustalw_fname):
        def which(file):
            for path in os.environ["PATH"].split(os.pathsep):
                if os.path.exists(os.path.join(path, file)):
                    return os.path.join(path, file)
            return None
        if not which("clustalw2"):
            print >> sys.stderr, "Error: clustalw2 is required, please have it installed on your system."
            sys.exit(1)
        clustalw2_cmd = ["clustalw2", "%s.l1hs.fasta" % rep_fname]
        print >> sys.stderr, "Running ClustalW:", ' '.join(clustalw2_cmd)
        alignview_proc = subprocess.Popen(clustalw2_cmd,
                                          stdout=open("/dev/null", 'w'),
                                          stderr=open("/dev/null", 'w'))
        alignview_proc.communicate()
        assert os.path.exists(clustalw_fname)

    allele_dic = {}
    for line in open(clustalw_fname):
        if not line.startswith("L1"):
            continue
        l1_name, _, _, aln = line.strip().split()
        aln = aln.replace('-', '.')
        if l1_name not in allele_dic:
            allele_dic[l1_name] = aln
        else:
            allele_dic[l1_name] += aln

    # Make sure every alignment is of the same length and
    #    replace 'N' to 'A' for now
    msf_len = None
    allele_dic2 = {}
    for l1_name, aln in allele_dic.items():
        # Make sure the length of allele ID is short, less than 20 characters
        assert len(l1_name) < 20
        if msf_len == None:
            msf_len = len(aln)
        else:
            assert msf_len == len(aln)
        aln2, non_nt = "", {}
        for nt in aln:
            if nt in "ACGT.":
                aln2 += nt
                continue
            else:
                if nt not in non_nt:
                    non_nt[nt] = 1
                else:
                    non_nt[nt] += 1
                aln2 += 'A'
        if non_nt:
            print >> sys.stderr, "Warning: %s has" % (l1_name),
            for nt, num in non_nt.items():
                print >> sys.stderr, "%d %s" % (num, nt),
            print >> sys.stderr
        
        allele_dic2[l1_name] = aln2
    allele_dic = allele_dic2
                
    # Write MSF (multiple sequence alignment file)
    msf_fname = "L1HS_gen.msf"
    msf_file = open(msf_fname, 'w')
    for s in range(0, msf_len, 50):
        for allele_id, msf in allele_dic.items():
            assert len(msf) == msf_len
            allele_name = "L1HS*%s" % (allele_id)
            print >> msf_file, "%20s" % allele_name,
            for s2 in range(s, min(msf_len, s + 50), 10):
                print >> msf_file, " %s" % msf[s2:s2+10],
            print >> msf_file
        if s + 50 >= msf_len:
            break
        print >> msf_file
    msf_file.close()

    # Write FASTA file
    fasta_fname = "L1HS_gen.fasta"
    fasta_file = open(fasta_fname, 'w')
    for allele_id, msf in allele_dic.items():
        gen_seq = msf.replace('.', '')
        print >> fasta_file, ">L1HS*%s %d bp" % (allele_id, len(gen_seq))
        for s in range(0, len(gen_seq), 60):
            print >> fasta_file, gen_seq[s:s+60]
    fasta_file.close()

    
"""
"""
if __name__ == '__main__':
    parser = ArgumentParser(
        description="Build epigenome")
    parser.add_argument("--base", "--base-fname",
                        dest="base_fname",
                        type=str,
                        default="l1hs",
                        help="base filename (default: rep)")
    parser.add_argument("-v", "--verbose",
                        dest="verbose",
                        action="store_true",
                        help="also print some statistics to stderr")

    args = parser.parse_args()
    extract_rep(args.base_fname,
                args.verbose)
    
