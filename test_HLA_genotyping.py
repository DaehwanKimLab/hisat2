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

    # Backbone sequence
    ref_seq = HLAs["A*01:01:01:01"]

    # Read HLA variants, and link information
    Vars, Var_list = {}, []
    for line in open("hla.snp"):
        var_id, var_type, allele, pos, data = line.strip().split('\t')
        assert not var_id in Vars
        Vars[var_id] = [var_type, int(pos), data]
        Var_list.append([int(pos), var_id])
    Var_list = sorted(Var_list)
    def lower_bound(Var_list, pos):
        low, high = 0, len(Var_list)
        while low < high:
            m = (low + high) / 2
            m_pos = Var_list[m][0]
            if m_pos < pos:
                low = m + 1
            elif m_pos > pos:
                high = m
            else:
                assert m_pos == pos
                while m > 0:
                    if Var_list[m-1][0] < pos:
                        break
                    m -= 1
                return m
        return low        
            
    Links = {}
    for line in open("hla.link"):
        var_id, alleles = line.strip().split('\t')
        alleles = alleles.split()
        assert not var_id in Links
        Links[var_id] = alleles

    # Test HLA genotyping
    basic_test, random_test = False, True
    test_passed = 0
    test_list = []
    if basic_test:
        for i in range(len(HLA_names)):
            test_list.append([HLA_names[i]])
    if random_test:
        test_size = 100
        allele_count = 2
        for test_i in range(test_size):
            nums = [i for i in range(len(HLA_names))]
            random.shuffle(nums)
            test_HLA_names = [HLA_names[nums[i]] for i in range(allele_count)]
            test_list.append(test_HLA_names)
    for test_i in range(len(test_list)):
        print >> sys.stderr, "Test %d" % (test_i + 1)
        test_HLA_names = test_list[test_i]
        
        # daehwan - for debugging purposes
        # test_HLA_names = ["A*80:01:01:01"]
        
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
        HLA_counts, HLA_cmpt = {}, {}
        for line in alignview_proc.stdout:
            cols = line[:-1].split()
            read_id, flag, chr, pos, mapQ, cigar_str = cols[:6]
            read_seq = cols[9]

            flag = int(flag)
            pos = int(pos)

            # daehwan - for debugging purposes
            debug = False
            if pos - 1 == 545 and False:
                debug = True

            if flag & 0x4 != 0:
                continue

            Zs, MD = "", ""
            for i in range(11, len(cols)):
                col = cols[i]
                if col.startswith("Zs"):
                    Zs = col[5:]
                elif col.startswith("MD"):
                    MD = col[5:]

            vars = []
            if Zs:
                vars = Zs.split(',')
                
            assert MD != ""
            MD_str_pos, MD_len = 0, 0
            read_pos, left_pos = 0, pos - 1
            right_pos = left_pos
            cigars = cigar_re.findall(cigar_str)
            cigars = [[cigar[-1], int(cigar[:-1])] for cigar in cigars]
            cmp_list = []
            for i in range(len(cigars)):
                cigar_op, length = cigars[i]
                if cigar_op == 'M':
                    first = True
                    MD_len_used = 0
                    while True:
                        if not first or MD_len == 0:
                            if MD[MD_str_pos].isdigit():
                                num = int(MD[MD_str_pos])
                                MD_str_pos += 1
                                while MD_str_pos < len(MD):
                                    if MD[MD_str_pos].isdigit():
                                        num = num * 10 + int(MD[MD_str_pos])
                                        MD_str_pos += 1
                                    else:
                                        break
                                MD_len += num
                        # Insertion or full match followed
                        if MD_len >= length:
                            MD_len -= length
                            cmp_list.append(["match", right_pos + MD_len_used, length - MD_len_used])
                            break
                        first = False
                        read_base = read_seq[read_pos + MD_len]
                        MD_ref_base = MD[MD_str_pos]
                        MD_str_pos += 1
                        assert MD_ref_base in "ACGT"
                        assert MD_ref_base != read_base
                        cmp_list.append(["match", right_pos + MD_len_used, MD_len - MD_len_used])
                        cmp_list.append(["mismatch", right_pos + MD_len, 1])
                        MD_len_used = MD_len + 1
                        MD_len += 1
                        # Full match
                        if MD_len == length:
                            MD_len = 0
                            break
                elif cigar_op == 'I':
                    cmp_list.append(["insertion", right_pos, length])
                elif cigar_op == 'D':
                    if MD[MD_str_pos] == '0':
                        MD_str_pos += 1
                    assert MD[MD_str_pos] == '^'
                    MD_str_pos += 1
                    while MD_str_pos < len(MD):
                        if not MD[MD_str_pos] in "ACGT":
                            break
                        MD_str_pos += 1
                    cmp_list.append(["deletion", right_pos, length])
                elif cigar_op == 'S':
                    cmp_list.append(["soft", right_pos, length])
                else:                    
                    assert cigar_op == 'N'
                    cmp_list.append(["intron", right_pos, length])
                
                if cigar_op in "MND":
                    right_pos += length

                if cigar_op in "MIS":
                    read_pos += length

            HLA_count_per_read = {}
            for HLA_name in HLA_names:
                HLA_count_per_read[HLA_name] = 0
               
            def add_count(var_id, add):
                assert var_id in Links
                alleles = Links[var_id]
                for allele in alleles:
                    HLA_count_per_read[allele] += add
            
            # Sanity check - read length, cigar string, and MD string
            # Decide which allele(s) a read most likely came from
            ref_pos, read_pos, cmp_cigar_str, cmp_MD = left_pos, 0, "", ""
            cigar_match_len, MD_match_len = 0, 0            
            for cmp in cmp_list:
                type = cmp[0]
                length = cmp[2]
                if type == "match":
                    var_idx = lower_bound(Var_list, ref_pos)
                    while var_idx < len(Var_list):
                        var_pos, var_id = Var_list[var_idx]
                        if ref_pos + length <= var_pos:
                            break
                        if ref_pos <= var_pos:
                            var_type, _, var_data = Vars[var_id]
                            if var_type == "insertion":
                                if ref_pos + length > var_pos + len(var_data):
                                    add_count(var_id, -1)
                                # daehwan - for debugging purposes
                                if debug:
                                    print cmp, var_id
                                    
                            else:
                                add_count(var_id, -1)
                        var_idx += 1
                    
                    read_pos += length
                    ref_pos += length
                    cigar_match_len += length
                    MD_match_len += length
                elif type == "mismatch":
                    read_base = read_seq[read_pos]
                    var_idx = lower_bound(Var_list, ref_pos)
                    while var_idx < len(Var_list):
                        var_pos, var_id = Var_list[var_idx]
                        if ref_pos < var_pos:
                            break
                        if ref_pos == var_pos:
                            var_type, _, var_data = Vars[var_id]
                            if var_type == "single":
                                # daehwan - for debugging purposes
                                if debug:
                                    print cmp, var_id, var_data, read_base
                                if var_data == read_base:
                                    add_count(var_id, 1)
                                else:
                                    add_count(var_id, -1)
                        var_idx += 1
                    
                    cmp_MD += ("%d%s" % (MD_match_len, ref_seq[ref_pos]))
                    MD_match_len = 0
                    cigar_match_len += 1
                    read_pos += 1
                    ref_pos += 1
                elif type == "insertion":
                    ins_seq = read_seq[read_pos:read_pos+length]
                    var_idx = lower_bound(Var_list, ref_pos)
                    # daehwan - for debugging purposes
                    if debug:
                        print left_pos, cigar_str, MD, vars
                        print ref_pos, ins_seq, Var_list[var_idx], Vars[Var_list[var_idx][1]]
                        # sys.exit(1)
                    while var_idx < len(Var_list):
                        var_pos, var_id = Var_list[var_idx]
                        if ref_pos < var_pos:
                            break
                        if ref_pos == var_pos:
                            var_type, _, var_data = Vars[var_id]
                            if var_type == "insertion":                                
                                if var_data == ins_seq:
                                    # daehwan - for debugging purposes
                                    if debug:
                                        print cmp, var_id
                                    add_count(var_id, 1)
                        var_idx += 1

                    if cigar_match_len > 0:
                        cmp_cigar_str += ("%dM" % cigar_match_len)
                        cigar_match_len = 0
                    read_pos += length
                    cmp_cigar_str += ("%dI" % length)
                elif type == "deletion":
                    var_idx = lower_bound(Var_list, ref_pos)
                    while var_idx < len(Var_list):
                        var_pos, var_id = Var_list[var_idx]
                        if ref_pos < var_pos:
                            break
                        if ref_pos == var_pos:
                            var_type, _, var_data = Vars[var_id]
                            if var_type == "deletion":
                                var_len = int(var_data)
                                if var_len == length:
                                    # daehwan - for debugging purposes
                                    if debug:
                                        print cmp, var_id
                                    add_count(var_id, 1)
                        var_idx += 1

                    if cigar_match_len > 0:
                        cmp_cigar_str += ("%dM" % cigar_match_len)
                        cigar_match_len = 0
                    cmp_MD += ("%d" % MD_match_len)
                    MD_match_len = 0
                    cmp_cigar_str += ("%dD" % length)
                    cmp_MD += ("^%s" % ref_seq[ref_pos:ref_pos+length])
                    ref_pos += length
                elif type == "soft":
                    if cigar_match_len > 0:
                        cmp_cigar_str += ("%dM" % cigar_match_len)
                        cigar_match_len = 0
                    read_pos += length
                    cmp_cigar_str += ("%dS" % length)
                else:
                    assert type == "intron"
                    if cigar_match_len > 0:
                        cmp_cigar_str += ("%dM" % cigar_match_len)
                        cigar_match_len = 0
                    cmp_cigar_str += ("%dN" % length)
                    ref_pos += length                    
            if cigar_match_len > 0:
                cmp_cigar_str += ("%dM" % cigar_match_len)
            cmp_MD += ("%d" % MD_match_len)
            if read_pos != len(read_seq) or \
                    cmp_cigar_str != cigar_str or \
                    cmp_MD != MD:
                print >> sys.stderr, "Error:", cigar_str, MD
                print >> sys.stderr, "\tcomputed:", cmp_cigar_str, cmp_MD
                print >> sys.stderr, "\tcmp list:", cmp_list
                assert False            

            max_count = None
            for allele, count in HLA_count_per_read.items():
                if not max_count:
                    max_count = count
                elif count > max_count:
                    max_count = count
            cur_cmpt = set()
            # daehwan - for debugging purposes
            allele1_found, allele2_found = False, False
            for allele, count in HLA_count_per_read.items():
                if count < max_count:
                    continue
                if allele == "A*68:02:01:01":
                    allele1_found = True
                elif allele == "A*68:18N":
                    allele2_found = True
                cur_cmpt.add(allele)                    
                if not allele in HLA_counts:
                    HLA_counts[allele] = 1
                else:
                    HLA_counts[allele] += 1

            # daehwan - for debugging purposes
            """
            if allele1_found != allele2_found:
                if allele1_found:
                    print ("read_id %s]" % read_id), left_pos, cigar_str, MD, Zs
                    print read_seq
            """

            cur_cmpt = sorted(list(cur_cmpt))
            cur_cmpt = '-'.join(cur_cmpt)
            if not cur_cmpt in HLA_cmpt:
                HLA_cmpt[cur_cmpt] = 1
            else:
                HLA_cmpt[cur_cmpt] += 1

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
                    print >> sys.stderr, "*** %d ranked %s (count: %d)" % (count_i + 1, test_HLA_name, count[1])
                    """
                    if count_i > 0 and HLA_counts[0][1] > count[1]:
                        print >> sys.stderr, "Warning: %s ranked first (count: %d)" % (HLA_counts[0][0], HLA_counts[0][1])
                        assert False
                    else:
                        test_passed += 1
                    """
            if count_i < 2 and False:
                print >> sys.stderr, "%d %s (count: %d)" % (count_i + 1, count[0], count[1])

        HLA_prob, HLA_prob_next = {}, {}
        for cmpt, count in HLA_cmpt.items():
            alleles = cmpt.split('-')
            for allele in alleles:
                if allele not in HLA_prob:
                    HLA_prob[allele] = 0.0
                HLA_prob[allele] += (float(count) / len(alleles))

        def normalize(prob):
            total = sum(prob.values())
            for allele, mass in prob.items():
                prob[allele] = mass / total

        normalize(HLA_prob)
        def prob_diff(prob1, prob2):
            diff = 0.0
            for allele in prob1.keys():
                assert allele in prob2
                diff += abs(prob1[allele] - prob2[allele])
            return diff

        def next_prob(HLA_cmpt, HLA_prob):
            HLA_prob_next = {}
            for cmpt, count in HLA_cmpt.items():
                alleles = cmpt.split('-')
                alleles_prob = 0.0
                for allele in alleles:
                    assert allele in HLA_prob
                    alleles_prob += HLA_prob[allele]
                for allele in alleles:
                    if allele not in HLA_prob_next:
                        HLA_prob_next[allele] = 0.0
                    HLA_prob_next[allele] += (float(count) * HLA_prob[allele] / alleles_prob)
            normalize(HLA_prob_next)
            return HLA_prob_next

        diff, iter = 1.0, 0
        while diff > 0.0001 and iter < 1000:
            HLA_prob_next = next_prob(HLA_cmpt, HLA_prob)
            diff = prob_diff(HLA_prob, HLA_prob_next)
            HLA_prob = HLA_prob_next
            iter += 1
        HLA_prob = [[allele, prob] for allele, prob in HLA_prob.items()]
        def HLA_prob_cmp(a, b):
            if a[1] != b[1]:
                if a[1] < b[1]:
                    return 1
                else:
                    return -1
            assert a[0] != b[0]
            if a[0] < b[0]:
                return -1
            else:
                return 1
            
        HLA_prob = sorted(HLA_prob, cmp=HLA_prob_cmp)
        success = [False for i in range(len(test_HLA_names))]
        found_list = [False for i in range(len(test_HLA_names))]
        for prob_i in range(len(HLA_prob)):
            prob = HLA_prob[prob_i]
            found = False
            for test_i in range(len(test_HLA_names)):
                test_HLA_name = test_HLA_names[test_i]
                if prob[0] == test_HLA_name:
                    print >> sys.stderr, "*** %d ranked %s (abundance: %.2f%%)" % (prob_i + 1, test_HLA_name, prob[1] * 100.0)
                    if prob_i < len(success):
                        success[prob_i] = True
                    found_list[test_i] = True
                    found = True                        
            if not False in found_list:
                break
            if not found:
                print >> sys.stderr, "%d ranked %s (abundance: %.2f%%)" % (prob_i + 1, prob[0], prob[1] * 100.0)
        if not False in success:
            test_passed += 1

    print >> sys.stderr, "%d/%d passed (%.2f)" % (test_passed, len(test_list), test_passed * 100.0 / len(test_list))
    
        
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
