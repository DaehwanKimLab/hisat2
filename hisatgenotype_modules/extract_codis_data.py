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
import inspect, operator
from argparse import ArgumentParser, FileType


# 13 CODIS sequences downloaded from http://www.cstl.nist.gov/biotech/strbase/seq_ref.htm 
CODIS_seq = {"CSF1PO" : "CAGCTGGGATGTGGAGTGGTGTGAGGAGTGGCCACAGGGGAGCAGAGGAGGTGGCAGAAGCCGGAGGTAAAGGTGTCTTAAAGTGAGAAAGAATAACTGCATCTTAACCTATTGGGAGGTCATTGTAAAGAGGAGAGTGATGGGGTCAGATTGTACAGAGGAGGCACTTCGTGGTGGTCAGGAGCACACACTCCAGGGCAGTGTTCCAACCTGAGTCTGCCAAGGACTAGCAGGTTGCTAACCACCCTGTGTCTCAGTTTTCCTACCTGTAAAATGAAGATATTAACAGTAACTGCCTTCATAGATAGAAGATAGATAGATTAGATAGATAGATAGATAGATAGATAGATAGATAGATAGATAGATAGATAGGAAGTACTTAGAACAGGGTCTGACACAGGAAATGCTGTCCAAGTGTGCACCAGGAGATAGTATCTGAGAAGGCTCAGTCTGGCACCATGTGGGTTGGGTGGGAACCTGGAGGCTGGAGAATGGGCTGAAGATGGCCAGTGGTGTGTGGAAGAGTCTGAGATGCAGGGATGAGGAAGAGAAAGGAGATAAGGATGACCTCCAGGTCTCTGGCTATGGTGATTGGGTGCA",
             "FGA" : "ACACCTTTAAAATTCCAAAGAAAGTTCTTCTTCTATATTTCTTTGGGATTACTAATTGCTATTAGGACATCTTAACTGGCATTCATGGAAGGCTGCAGGGCATAACATTATCCAAAAGTCAAATGCCCCATAGGTTTTGAACTCACAGATTAAACTGTAACCAAAATAAAATTAGGCATATTTACAAGCTAGTTTCTTTCTTTCTTTTTTCTCTTTCTTTCTTTCTTTCTTTCTTTCTTTCTTTCTTTCTTTCTTTCTTTCTTTCTCCTTCCTTCCTTTCTTCCTTTCTTTTTTGCTGGCAATTACAGACAAATCACTCAGCAGCTACTTCAATAACCATATTTTCGATTTCAGACCGTGATAATACCTACAACCGAGTGTCAGAGGATCTGAGAAGCAGAATTGAAGTCCTGAAGCGCAAAGTCATAGAAAAAGTACAGCATATCCAGCTTCTGCAGAAAAATGTTAGAGCTCAGTTGGTTGATATGAAACGACTGGAGGTAAGTATGTGGCTGTGGTCCCGAGTGTCCTTGTTTTTGA",
             "TH01" : "TAAATAAAAACCCATTTTAAATGTGCCAGGGAGCCCAAGGTTCTGAGTGCCCAAGGAGGCACCGAAGACCCCTCCTGTGGGCTGAAAAGCTCCCGATTATCCAGCCTGGCCCACACAGTCCCCTGTACACAGGGCTTCCGAGTGCAGGTCACAGGGAACACAGACTCCATGGTGAATGAATGAATGAATGAATGAATGAATGAATGAATGAGGGAAATAAGGGAGGAACAGGCCAATGGGAATCACCCCAGAGCCCAGATACCCTTTGAATTTTGCCCCCTATTTGCCCAGGACCCCCCACCATGAGCTGCTGCTAGAGCCTGGGAAGGGCCTTGGGGCTGCCTCCCCAAGCAGGCAGGCTGGTTGGGGTGCTGACTAGGGCAGCTGGGGCAGAGGGAGGCAGGGGCAGGTGGGAGTAGG",
             "TPOX" : "CCACTGCTGACTTCCCCAAGCTAACTGTGCCACAGAGTGGGGACCCCCTCCCAGCTCTCACAACCCCCACCTTCCTCTGCTTCACTTTTCACCAACTGAAATATGGCCAAAGGCAAAAACCCATGTTCCCACTGGCCTGTGGGTCCCCCCATAGATCGTAAGCCCAGGAGGAAGGGCTGTGTTTCAGGGCTGTGATCACTAGCACCCAGAACCGTCGACTGGCACAGAACAGGCACTTAGGGAACCCTCACTGAATGAATGAATGAATGAATGAATGAATGAATGAATGAATGAATGTTTGGGCAAATAAACGCTGACAAGGACAGAAGGGCCTAGCGGGAAGGGAACAGGAGTAAGACCAGCGCACAGCCCGACTTGTGTTCAGAAGACCTGGGATTGGACCTGAGGAGTTCAATTTTGGATGAATCTCTTAATTAACCTGTGTGGTTCCCAGTTCCTCCCCTGAGCGCCCAGGACAGTAGAGTCAACCTCACGTTTGAGCGTTGGGGACGCAAACACGAGAGTGCTTGGTGTGAGCAC",
             "VWA" : "AGCAGTAGAGAGAACTAGAGGGATCATTTACTTCAAGCCCCTCATTTTATAGACATTACTAGTCTCCTACAATGTGCCGGGCACTTTGCCCTTATTATTTTGTGAACTCCTCAGACTGATCCTATAAGGTAGAGTTCCCACCTTCCAGAAGAAGAAACAGGTCTAGAGGATCCAAGTTGACTTGGCTGAGATGTGAAAGCCCTAGTGGATGATAAGAATAATCAGTATGTGACTTGGATTGATCTATCTGTCTGTCTGTCTGTCTATCTATCTATCTATCTATCTATCTATCTATCTATCTATCTATCTATCTATCCATCTATCCATCCATCCTATGTATTTATCATCTGTCCTATCTCTATCTAACCTATGTATCTATTTATCATCTATCCTGTCTCTATCTATCCTTTGTATCTATCA",
             "D3S1358" : "AGTGGAAAAGCTATTCCCAGGTGAGGACTGCAGCTGCCAGGGCACTGCTCCAGAATGGGCATGCTGGCCATATTCACTTGCCCACTTCTGCCCAGGGATCTATTTTTCTGTGGTGTGTATTCCCTGTGCCTTTGGGGGCATCTCTTATACTCATGAAATCAACAGAGGCTTGCATGTATCTATCTGTCTGTCTATCTATCTATCTATCTATCTATCTATCTATCTATCTATCTATCTATCTATCTATCTATGAGACAGGGTCTTGCTCTGTCACCCAGATTGGACTGCAGTGGGGGAATCATAGCTCACTACAGCCTCAAACTCCTGGGCTCAAGCAGTCCTCCTGCCTCAGCCTCCCAAGTACCTGGGATTATAGGCATGAGCCACCATGTCCGGCTAATTTTTTTTTTTAAGAGATGG",
             "D5S818" : "CATTTGTAATTTTCTAATTCAAAGGAGTATATAATTATGTAATAATTTTAAAATTAAATACTGAGACATGCATATGCTTTTAAAGCTTCTAATTAAAGTGGTGTCCCAGATAATCTGTACTAATAAAAGTATATTTTAATAGCAAGTATGTGACAAGGGTGATTTTCCTCTTTGGTATCCTTATGTAATATTTTGAAGATAGATAGATAGATAGATAGATAGATAGATAGATAGATAGATAGAGGTATAAATAAGGATACAGATAAAGATACAAATGTTGTAAACTGTGGCTATGATTGGAATCACTTGGCTAAAAAGCACTAAAGCATTCCTCTGAGAGAGACAATTACTTTTTTGCTTAGGAAACTACCTCAACAGCCTATTAGCATCTGAAATATGAGGTCCACTATCCAGATGGGAGAGGTTTAGAAAAAGAAGACTTATATTACTCTGTATAATGAAATGATGGAGTATTTGGAG",
             "D7S820" : "ACTGCAACCTCCGCTTCTTGGGTCAAGTGGTTCTCCTGCCCCAGCCTCCTGAGTAGCTGGGACTACAGGCATGTGCTACTGCATCCAGCTAATTTTTGTATTTTTTTTAGAGACGGGGTTTCACCATGTTGGTCAGGCTGACTATGGAGTTATTTTAAGGTTAATATATATAAAGGGTATGATAGAACACTTGTCATAGTTTAGAACGAACTAACGATAGATAGATAGATAGATAGATAGATAGATAGATAGATAGATAGATAGATAGACAGATTGATAGTTTTTTTTAATCTCACTAAATAGTCTATAGTAAACATTTAATTACCAATATTTGGTGCAATTCTGTCAATGAGGATAAATGTGGAATCGTTATAATTCTTAAGAATATATATTCCCTCTGAGTTTTTGATACCTCAGATTTTAAGACCTCACAATTATCTCATAAGGCTTAAAATCAATCATATTTTGAGGATCAACTTATGGTATTTTTGCCTGTTTTTATTCCTTCTGGTGTGAAAACTGATGCCTTCCATCGTGTAA",
             "D8S1179" : "TCTGAAGTAAGTAAAACATTGATACAGGATCCTTGGGGTGTCGCTTTTTTGGCCAGAAACCTCTGTAGCCAGTGGCGCCTTTGCCTGAGTTTTGCTCAGGCCCACTGGGCTCTTTTTGCCCACACGGCCGGGCAACTTATATGTATTTTTGTATTTCATGTGTACATTCGTATCTATCTGTCTATCTATCTATCTATCTATCTATCTATCTATCTATCTATCTATTCCCCACAGTGAAAATAATCTACAGGATAGGTAAATAAATTAAGGCATATTCACGCAATGGGATACGATACAGTGATGAAAATGAACTAATTATAGCTACGTGAAACTATACTCATGAACACAATTTGGTAAAAGAAACTGGAAACAAGAATACATACGGTTTTTGACAGCTGTACTATTTTACATTCCCAACAA",
             "D13S317" : "GTATTGCAAGCACTTAGTTACATTTCTAGCATATAACACATGATCAATAAATATTTTGACATGAACAAATGGTAATTCTGCCTACAGCCAATGTGAATATTGGGATGGGTTGCTGGACATGGTATCACAGAAGTCTGGGATGTGGAGGAGAGTTCATTTCTTTAGTGGGCATCCGTGACTCTCTGGACTCTGACCCATCTAACGCCTATCTGTATTTACAAATACATTATCTATCTATCTATCTATCTATCTATCTATCTATCTATCTATCAATCAATCATCTATCTATCTTTCTGTCTGTCTTTTTGGGCTGCCTATGGCTCAACCCAAGTTGAAGGAGGAGATTTGACCAACAATTCAAGCTCTCTGAATATGTTTTGAAAATAATGTATATTAATGAATGTACAAATTTCCCCACTTGTACTTTCAGACTGTTATCTGTGAGTTAAAACTCCTCCACTCTTTTTCCTACCCAAATAA",
             "D16S539" : "GGAGGATGACTGTGTTCCCACTCTCAGTCCTGCCGAGGTGCCTGACAGCCCTGCACCCAGGAGCTGGGGGGTCTAAGAGCTTGTAAAAAGTGTACAAGTGCCAGATGCTCGTTGTGCACAAATCTAAATGCAGAAAAGCACTGAAAGAAGAATCCCGAAAACCACAGTTCCCATTTTTATATGGGAGCAAACAAAGCAGATCCCAAGCTCTTCCTCTTCCCTAGATCAATACAGACAGACAGACAGGTGGATAGATAGATAGATAGATAGATAGATAGATAGATAGATAGATATCATTGAAAGACAAAACAGAGATGGATGATAGATACATGCTTACAGATGCACACACAAACGCTAAATGGTATAAAAATGGAATCACTCTGTAGGCTGTTTTACCACCTACTTTACTAAATTAATGAGTTATTGAGTATAATTTAATTTTATATACTAATTTGAAACTGTGTCATTAGGTTTTTAAGT",
             "D18S51" : "TACAAAAAAATACAAAAATTAGTTGGGCATGGTGGCACGTGCCTGTAGTCTCAGCTACTTGCAGGGCTGAGGCAGGAGGAGTTCTTGAGCCCAGAAGGTTAAGGCTGCAGTGAGCCATGTTCATGCCACTGCACTTCACTCTGAGTGACAAATTGAGACCTTGTCTCAGAAAGAAAGAAAGAAAGAAAGAAAGAAAGAAAGAAAGAAAGAAAGAAAGAAAGAAAGAAAGAAAGAAAGAAAAAGAGAGAGGAAAGAAAGAGAAAAAGAAAAGAAATAGTAGCAACTGTTATTGTAAGACATCTCCACACACCAGAGAAGTTAATTTTAATTTTAACATGTTAAGAACAGAGAGAAGCCAACATGTCCACCTTAGGCTGACGGTTTGTTTATTTGTGTTGTTGCTGGTAGTCGGGTTTGTTATTTTTAAAGTAGCTTATCCAATACTTCATTAACAATTTCAGTAAGTTATTTCATCTTTCAACATAAATACGCACAAGGATTTCTTCTGGTCAAGACCAAACTAATATTAGTCCATAGTAG",
             "D21S11" : "AACTGAAGGTTACACTATCAGCTTCCGTTGTTCTAAGGGCTTCAGACTTGGACAGCCACACTGCCAGCTTCCCTGATTCTTCAGCTTGTAGATGGTCTGTTATGGGACTTTTCTCAGTCTCCATAAATATGTGAGTCAATTCCCCAAGTGAATTGCCTTCTATCTATCTATCTATCTGTCTGTCTGTCTGTCTGTCTGTCTATCTATCTATATCTATCTATCTATCATCTATCTATCCATATCTATCTATCTATCTATCTATCTATCTATCTATCTATCTATCTATCGTCTATCTATCCAGTCTATCTACCTCCTATTAGTCTGTCTCTGGAGAACATTGACTAATACAACATCTTTAATATATCACAGTTTAATTTCAAGTTATATCATACCACTTCATACATTATATAAAACCTTACAGTGTTTCTCCCTTCTCAGTGTTTATGGCTAGTAATTTTTTACTGGGTGCCAGACACTAAT"}

CODIS_ref_name = {}


"""
## Download variant information from website
"""
def get_html(url):
    download_cmd = ["wget",
                    "-O", "-",
                    url]
    proc = subprocess.Popen(download_cmd,
                            stdout=subprocess.PIPE,
                            stderr=open("/dev/null", 'w'))

    output = ""
    for line in proc.stdout:
        output += line

    return output


"""
Extract multiple sequence alignments
"""
def extract_msa(base_dname,
                base_fname,
                verbose):
    # CODIS database base URL
    base_url = "http://www.cstl.nist.gov/biotech/strbase"
    
    # Refer to Python's regular expression at https://docs.python.org/2/library/re.html
    #   <td width="16%" align="center"><font size="4">47.2 </font> </td>
    allele_re = re.compile('>(\d+\.?\d?\"?\'*\(?\d*\.?\d?\"?\'*\)?\*?)</')
    #   <td width="35%"><font size="2">[TTTC]<sub>4</sub>TTTT TT<span style="mso-spacerun: yes"> </span>[CTTT]<sub>14</sub>[CTGT]<sub>3</sub>[CTTT]<sub>14 </sub>[CTTC]<sub>4</sub>[CTTT]<sub>3</sub>CTCC[TTCC]<sub>4</sub></font> </td>
    # repeat_re = re.compile('^(\[[ACGT]+\]\d+|[ACGT]+)+$')
    repeat_re = re.compile('^(\[[ACGT]+\]\d+|[ACGT]+)+$')
    # Remove extra tags
    tag_re = re.compile('(<[^>]*>)')
    nbsp_re = re.compile('&nbsp;')
    quot_re = re.compile('&quot;')
    for locus_name in CODIS_seq.keys():
        url = "%s/str_%s.htm" % (base_url, locus_name)
        content = get_html(url).split("\r\n")
        content = map(lambda x: x.strip().replace(' ',''), content)
        content2 = []
        for line in content:
            if line.startswith("<t") or \
               line.startswith("</tr") or \
               len(content2) == 0:
                content2.append(line)
            else:
                content2[-1] += line

        content = content2
        alleles = []
        l = 0
        while l < len(content):
            line = content[l]
            if line.startswith("<tr"):
                l += 1
                if l < len(content):
                    line = content[l]
                    line = re.sub(nbsp_re, '', line)
                    line = re.sub(quot_re, "''", line)
                    allele_match = allele_re.search(line)
                    if not allele_match:
                        continue
                    allele_id = allele_match.group(1)
                    l += 1
                    repeat_match = None
                    while l < len(content):
                        line = content[l]                        
                        if not line.startswith("<td"):
                            break
                        line = re.sub(tag_re, '', line)
                        line = re.sub(nbsp_re, '', line)
                        repeat_match = repeat_re.search(line)
                        if repeat_match:
                            break
                        l += 1
                        
                    if not repeat_match:
                        continue

                    repeat_st = line
                    alleles.append([allele_id, repeat_st])
            else:
                l += 1

        if len(alleles) <= 0:
            continue

        # From   [TTTC]3TTTTTTCT[CTTT]20CTCC[TTCC]2
        # To     [['TTTC', 3], ['TTTTTTCT', 1], ['CTTT', 20], ['CTCC', 1], ['TTCC', 2]]
        def read_allele(repeat_st):
            allele = []
            s = 0
            while s < len(repeat_st):
                ch = repeat_st[s]
                assert ch in "[ACGT"
                if ch == '[':
                    s += 1
                    repeat = ""
                    while s < len(repeat_st):
                        nt = repeat_st[s]
                        if nt in "ACGT":
                            repeat += nt
                            s += 1
                        else:
                            assert nt == ']'
                            s += 1
                            break
                    assert s < len(repeat_st)
                    num = 0
                    while s < len(repeat_st):
                        digit = repeat_st[s]
                        if digit.isdigit():
                            num = num * 10 + int(digit)
                            s += 1
                        else:
                            break
                    assert num > 0
                    allele.append([repeat, num])
                else:
                    repeat = ""
                    while s < len(repeat_st):
                        nt = repeat_st[s]
                        if nt in "ACGT":
                            repeat += nt
                            s += 1
                        else:
                            assert nt == '['
                            break
                    allele.append([repeat, 1])

            # Sanity check
            cmp_repeat_st = ""
            for repeat, repeat_num in allele:
                if repeat_num > 1 or locus_name == "D8S1179":
                    cmp_repeat_st += "["
                cmp_repeat_st += repeat
                if repeat_num > 1 or locus_name == "D8S1179":
                    cmp_repeat_st += "]%d" % repeat_num
                        
            assert repeat_st == cmp_repeat_st                        
            return allele

        alleles = [[allele_id, read_allele(repeat_st)] for allele_id, repeat_st in alleles]

        def to_sequence(repeat_st):
            sequence = ""
            for repeat, repeat_num in repeat_st:
                sequence += (repeat * repeat_num)
            return sequence

        allele_seqs = [[allele_id, to_sequence(repeat_st)] for allele_id, repeat_st in alleles]
        extended_ref_allele_seq = CODIS_seq[locus_name]
        for allele_id, allele_seq in allele_seqs:
            if extended_ref_allele_seq.find(allele_seq) != -1:
                CODIS_ref_name[locus_name] = allele_id
        if locus_name not in CODIS_ref_name:
            CODIS_ref_name[locus_name] = "undefined"

        print >> sys.stderr, "%s: %d alleles with reference allele as %s" % (locus_name, len(alleles), CODIS_ref_name[locus_name])
        if verbose:
            print >> sys.stderr, "\t", extended_ref_allele_seq
            for allele_id, allele in alleles:
                print >> sys.stderr, allele_id, "\t", allele

        ### Perform ClustalW multiple sequence alignment

        # Create a temporary allele sequence file
        seq_fname = "%s.tmp.fa" % locus_name
        msf_fname = "%s.tmp.aln" % locus_name
        dnd_fname = "%s.tmp.dnd" % locus_name
        seq_file = open(seq_fname, 'w')
        for allele_id, allele_seq in allele_seqs:
            print >> seq_file, ">%s" % allele_id
            print >> seq_file, allele_seq
        seq_file.close()

        # Run ClustalW - see http://www.clustal.org/clustal2 for more details
        clustalw_cmd = ["clustalw2", seq_fname]
        try:
            clustalw_proc = subprocess.Popen(clustalw_cmd,
                                             stdout=open("/dev/null", 'w'),
                                             stderr=open("/dev/null", 'w'))
            clustalw_proc.communicate()
            if not os.path.exists(msf_fname):
                print >> sys.stderr, "Error: running ClustalW failed."
        except:
            print >> sys.stderr, "Error: please install the latest version of ClustalW."

        allele_dic = {}
        for allele_id, allele_seq in allele_seqs:
            allele_dic[allele_id] = allele_seq

        allele_repeat_msf = {}
        for line in open(msf_fname):
            line = line.strip()
            if len(line) == 0 or \
               line.startswith("CLUSTAL") or \
               line.startswith("*"):
                continue
            
            tmp_allele_id, repeat = line.split()
            # ClustalW sometimes changes allele names
            allele_id = ""
            parenthesis_open = True
            for ch in tmp_allele_id:
                if ch == '_':
                    if parenthesis_open:
                        allele_id += '('
                    else:
                        allele_id += ')'
                    parenthesis_open = not parenthesis_open
                else:
                    allele_id += ch
            assert parenthesis_open
            repeat = repeat.replace('-', '.')
            if allele_id not in allele_repeat_msf:
                allele_repeat_msf[allele_id] = repeat
            else:
                allele_repeat_msf[allele_id] += repeat

        # Sanity check
        assert len(allele_dic) == len(allele_repeat_msf)
        repeat_len = -1
        for repeat_msf in allele_repeat_msf.values():
            if repeat_len < 0:
                repeat_len = len(repeat_msf)
            else:
                assert repeat_len == len(repeat_msf)

        # Creat full multiple sequence alignment
        ref_allele_id = CODIS_ref_name[locus_name]
        ref_allele_seq = allele_dic[ref_allele_id]
        replace_left = extended_ref_allele_seq.find(ref_allele_seq)
        replace_right = replace_left + len(ref_allele_seq)
        allele_msf = {}
        for allele_id, repeat_msf in allele_repeat_msf.items():
            allele_msf[allele_id] = extended_ref_allele_seq[:replace_left] + repeat_msf + extended_ref_allele_seq[replace_right:]
            
        os.remove(seq_fname)
        os.remove(msf_fname)
        os.remove(dnd_fname)

        # Make sure the length of allele ID is short, less than 20 characters
        max_allele_id_len = max([len(allele_id) for allele_id in allele_dic.keys()])
        assert max_allele_id_len < 20

        # Write MSF (multiple sequence alignment file)
        msf_len = len(extended_ref_allele_seq) - len(ref_allele_seq) + repeat_len
        msf_fname = "%s_gen.msf" % locus_name
        msf_file = open(msf_fname, 'w')
        for s in range(0, msf_len, 50):
            for allele_id, msf in allele_msf.items():
                assert len(msf) == msf_len
                allele_name = "%s*%s" % (locus_name, allele_id)
                print >> msf_file, "%20s" % allele_name,
                for s2 in range(s, min(msf_len, s + 50), 10):
                    print >> msf_file, " %s" % msf[s2:s2+10],
                print >> msf_file
                
            if s + 50 >= msf_len:
                break
            print >> msf_file

        msf_file.close()


        # Write FASTA file
        fasta_fname = "%s_gen.fasta" % locus_name
        fasta_file = open(fasta_fname, 'w')
        for allele_id, allele_seq in allele_seqs:
            gen_seq = extended_ref_allele_seq[:replace_left] + allele_seq + extended_ref_allele_seq[replace_right:]
            print >> fasta_file, ">%s*%s %d bp" % (locus_name, allele_id, len(gen_seq))
            for s in range(0, len(gen_seq), 60):
                print >> fasta_file, gen_seq[s:s+60]
        fasta_file.close()


"""
"""
if __name__ == '__main__':
    parser = ArgumentParser(
        description="Extract multiple sequence alignments for DNA Fingerprinting loci")
    parser.add_argument("-b", "--base",
                        dest="base_fname",
                        type=str,
                        default="codis",
                        help="base filename (default: codis)")
    parser.add_argument("-v", "--verbose",
                        dest="verbose",
                        action="store_true",
                        help="also print some statistics to stderr")

    args = parser.parse_args()
    if args.base_fname.find('/') != -1:
        elems = args.base_fname.split('/')
        base_fname = elems[-1]
        base_dname = '/'.join(elems[:-1])
    else:
        base_fname = args.base_fname
        base_dname = ""
        
    extract_msa(base_dname,
                base_fname,                
                args.verbose)

