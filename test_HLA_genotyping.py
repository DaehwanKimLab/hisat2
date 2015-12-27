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


import sys, os, subprocess, re
import inspect
import random
from argparse import ArgumentParser, FileType


"""
"""
def test_HLA_genotyping(base_fname, verbose = False):
    # Current directory
    curr_script = os.path.realpath(inspect.getsourcefile(test_HLA_genotyping))
    ex_path = os.path.dirname(curr_script)
    
    # Extract HLA variants, backbone sequence, and other sequeces
    HLA_fnames = ["hla_backbone.fa",
                  "hla_sequences.fa",
                  "hla.snp",
                  "hla.link"]

    def check_files(fnames):
        for fname in fnames:
            if not os.path.exists(fname):
                return False
        return True

    if not check_files(HLA_fnames):
        extract_hla_script = os.path.join(ex_path, "extract_HLA_vars.py")
        extract_cmd = [extract_hla_script,
                       "IMGTHLA/msf/A_gen.msf"]
        proc = subprocess.Popen(extract_cmd, stderr=open("/dev/null", 'w'))
        proc.communicate()
        if not check_files(HLA_fnames):
            print >> sys.stderr, "Error: extract_HLA_vars failed!"
            sys.exit(1)

    # Build HISAT2 indexes based on the above information
    HLA_index_fnames = ["hla.%d.ht2" % (i+1) for i in range(8)]
    if not check_files(HLA_index_fnames):
        hisat2_build = os.path.join(ex_path, "hisat2-build")
        build_cmd = [hisat2_build,
                     "--snp",
                     "hla.snp",
                     "hla_backbone.fa",
                     "hla"]
        proc = subprocess.Popen(build_cmd, stderr=open("/dev/null", 'w'))
        proc.communicate()        
        if not check_files(HLA_index_fnames):
            print >> sys.stderr, "Error: indexing HLA failed!"
            sys.exit(1)   

    # Read HLA alleles (names and sequences)
    HLAs = {}
    for line in open("hla_sequences.fa"):
        if line.startswith(">"):
            HLA_name = line.strip().split()[0][1:]
        else:
            if not HLA_name in HLAs:
                HLAs[HLA_name] = ""
            HLAs[HLA_name] += line.strip()
    HLA_names = list(HLAs.keys())

    # Read HLA variants, and link information
    Vars, Var_list = {}, []
    for line in open("hla.snp"):
        var_id, var_type, allele, pos, data = line.strip().split('\t')
        assert not var_id in Vars
        Vars[var_id] = [var_type, pos, data]
        Var_list.append([pos, var_id])
            
    Links = {}
    for line in open("hla.link"):
        var_id, alleles = line.strip().split('\t')
        alleles = alleles.split()
        assert not var_id in Links
        Links[var_id] = alleles

    # Test HLA genotyping
    test_size = 20
    allele_count = 1
    for test_i in range(test_size):
        print >> sys.stderr, "Test %d" % (test_i + 1)
        nums = [i for i in range(len(HLA_names))]
        random.shuffle(nums)
        test_HLA_names = [HLA_names[nums[i]] for i in range(allele_count)]
        for test_HLA_name in test_HLA_names:
            print >> sys.stderr, "\t%s" % (test_HLA_name)

        # Simulate reads from two HLA alleles
        def simulate_reads(seq, read_len = 100):
            reads = []
            for i in range(0, len(seq) - read_len + 1):
                reads.append(seq[i:i+read_len])
            return reads
        HLA_reads = []
        for test_HLA_name in test_HLA_names:
            HLA_seq = HLAs[test_HLA_name]
            HLA_reads += simulate_reads(HLA_seq)

        # Write reads into a fasta read file
        read_file = open('hla_input.fa', 'w')
        for read_i in range(len(HLA_reads)):
            print >> read_file, ">%d" % (read_i + 1)
            print >> read_file, HLA_reads[read_i]
        read_file.close()

        # Align reads, and sort the alignments into a BAM file
        hisat2 = os.path.join(ex_path, "hisat2")
        hisat2_cmd = [hisat2,
                      "-x",
                      "hla",
                      "-f",
                      "hla_input.fa"]
        align_proc = subprocess.Popen(hisat2_cmd,
                                      stdout=subprocess.PIPE,
                                      stderr=open("/dev/null", 'w'))
        sambam_cmd = ["samtools",
                      "view",
                      "-bS",
                      "-"]

        sambam_proc = subprocess.Popen(sambam_cmd,
                                       stdin=align_proc.stdout,
                                       stdout=open("hla_input_unsorted.bam", 'w'),
                                       stderr=open("/dev/null", 'w'))
        sambam_proc.communicate()

        bamsort_cmd = ["samtools",
                       "sort",
                       "hla_input_unsorted.bam",
                       "hla_input"]
        bamsort_proc = subprocess.Popen(bamsort_cmd,
                                        stderr=open("/dev/null", 'w'))
        bamsort_proc.communicate()

        bamindex_cmd = ["samtools",
                        "index",
                        "hla_input.bam"]
        bamindex_proc = subprocess.Popen(bamindex_cmd,
                                         stderr=open("/dev/null", 'w'))
        bamindex_proc.communicate()

        os.system("rm hla_input_unsorted.bam")

        # Read alignments
        alignview_cmd = ["samtools",
                         "view",
                         "hla_input.bam",
                         "A*01:01:01:01:0-10000"]
        alignview_proc = subprocess.Popen(alignview_cmd,
                                          stdout=subprocess.PIPE,
                                          stderr=open("/dev/null", 'w'))

        # Cigar regular expression
        cigar_re = re.compile('\d+\w')

        # Count alleles
        HLA_counts = {}
        for line in alignview_proc.stdout:
            cols = line[:-1].split()
            read_id, flag, chr, pos, mapQ, cigar_str = cols[:6]
            read_seq = cols[9]

            flag = int(flag)
            pos = int(pos)

            if flag & 0x4 != 0:
                continue

            Zs, MD = "", ""
            for i in range(11, len(cols)):
                col = cols[i]
                if col.startswith("Zs"):
                    Zs = col[5:]
                elif col.startswith("MD"):
                    MD = col[5:]

            assert MD != ""
            MD_pos, MD_len = 0, 0
            read_pos, right_pos = 0, pos - 1
            cigars = cigar_re.findall(cigar_str)
            cigars = [[cigar[-1], int(cigar[:-1])] for cigar in cigars]
            for i in range(len(cigars)):
                cigar_op, length = cigars[i]
                if cigar_op == 'M':
                    first = True
                    while True:
                        if not first or MD_len == 0:
                            assert MD[MD_pos].isdigit()
                            num = int(MD[MD_pos])
                            MD_pos += 1
                            while MD_pos < len(MD):
                                if MD[MD_pos].isdigit():
                                    num = num * 10 + int(MD[MD_pos])
                                    MD_pos += 1
                                else:
                                    break
                            MD_len += num
                        # Insertion or full match followed
                        if MD_len >= length:
                            MD_len -= length
                            break
                        first = False
                        MD_ref_base = MD[MD_pos]
                        MD_pos += 1
                        assert MD_ref_base in "ACGT"
                        MD_len += 1
                        # Full match
                        if MD_len == length:
                            MD_len = 0
                            break
                elif cigar_op == 'I':
                    None
                elif cigar_op == 'D':
                    if MD[MD_pos] == '0':
                        MD_pos += 1
                    assert MD[MD_pos] == '^'
                    MD_pos += 1
                    while MD_pos < len(MD):
                        if not MD[MD_pos] in "ACGT":
                            break
                        MD_pos += 1                            
                elif cigar_op == 'S':
                    None
                else:
                    assert cigar_op == 'N'
                
                if cigar_op in "MND":
                    right_pos += length

                if cigar_op in "MIS":
                    read_pos += length

            vars = []
            HLA_count_per_read = {}
            max_count = 0
            if Zs:
                vars = Zs.split(',')
            for var in vars:
                var = var.split('|')[2]
                assert var in Links
                alleles = Links[var]
                for allele in alleles:
                    if not allele in HLA_count_per_read:
                        HLA_count_per_read[allele] = 1
                    else:
                        HLA_count_per_read[allele] += 1
                    if HLA_count_per_read[allele] > max_count:
                        max_count = HLA_count_per_read[allele]

            for allele, count in HLA_count_per_read.items():
                if count < max_count:
                    continue
                if not allele in HLA_counts:
                    HLA_counts[allele] = 1
                else:
                    HLA_counts[allele] += 1

        HLA_counts = [[allele, count] for allele, count in HLA_counts.items()]
        def HLA_count_cmp(a, b):
            if a[1] != b[1]:
                return b[1] - a[1]
            assert a[0] != b[0]
            if a[0] < b[0]:
                return -1
            else:
                return 1

        HLA_counts = sorted(HLA_counts, cmp=HLA_count_cmp)
        for count_i in range(len(HLA_counts)):
            count = HLA_counts[count_i]
            for test_HLA_name in test_HLA_names:
                if count[0] == test_HLA_name:
                    print >> sys.stderr, "%s ranked %d (count: %d)" % (test_HLA_name, count_i + 1, count[1])
            if count_i < 1:
                print >> sys.stderr, "%d %s (count: %d)" % (count_i + 1, count[0], count[1])
               
        
"""
"""
if __name__ == '__main__':
    parser = ArgumentParser(
        description='test HLA genotyping')
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
    """
    if not args.HLA_MSA_file:
        parser.print_help()
        exit(1)
    """
    random.seed(1)
    test_HLA_genotyping(args.base_fname, args.verbose)
