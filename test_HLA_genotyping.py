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
def test_HLA_genotyping(base_fname, sequence_type, verbose = False):
    # Current directory
    curr_script = os.path.realpath(inspect.getsourcefile(test_HLA_genotyping))
    ex_path = os.path.dirname(curr_script)

    # Clone a git repository, IMGTHLA
    if not os.path.exists("IMGTHLA"):
        os.system("git clone https://github.com/jrob119/IMGTHLA.git")
    
    # Extract HLA variants, backbone sequence, and other sequeces
    HLA_fnames = ["hla_backbone.fa",
                  "hla_sequences.fa",
                  "hla.snp",
                  "hla.haplotype",
                  "hla.link"]

    def check_files(fnames):
        for fname in fnames:
            if not os.path.exists(fname):
                return False
        return True

    if not check_files(HLA_fnames):
        extract_hla_script = os.path.join(ex_path, "extract_HLA_vars.py")
        extract_cmd = [extract_hla_script,
                       "--gap", "30",
                       "--split", "50",
                       "--sequence-type", sequence_type]
        proc = subprocess.Popen(extract_cmd, stdout=open("/dev/null", 'w'), stderr=open("/dev/null", 'w'))
        proc.communicate()
        if not check_files(HLA_fnames):
            print >> sys.stderr, "Error: extract_HLA_vars failed!"
            sys.exit(1)

    # Build HISAT2 graph indexes based on the above information
    HLA_hisat2_graph_index_fnames = ["hla.graph.%d.ht2" % (i+1) for i in range(8)]
    if not check_files(HLA_hisat2_graph_index_fnames):
        hisat2_build = os.path.join(ex_path, "hisat2-build")
        build_cmd = [hisat2_build,
                     "--snp", "hla.snp",
                     "--haplotype", "hla.haplotype",
                     "hla_backbone.fa",
                     "hla.graph"]
        proc = subprocess.Popen(build_cmd, stdout=open("/dev/null", 'w'), stderr=open("/dev/null", 'w'))
        proc.communicate()        
        if not check_files(HLA_hisat2_graph_index_fnames):
            print >> sys.stderr, "Error: indexing HLA failed!"
            sys.exit(1)

    # Build HISAT2 linear indexes based on the above information
    HLA_hisat2_linear_index_fnames = ["hla.linear.%d.ht2" % (i+1) for i in range(8)]
    if not check_files(HLA_hisat2_linear_index_fnames):
        hisat2_build = os.path.join(ex_path, "hisat2-build")
        build_cmd = [hisat2_build,
                     "hla_sequences.fa",
                     "hla.linear"]
        proc = subprocess.Popen(build_cmd, stdout=open("/dev/null", 'w'), stderr=open("/dev/null", 'w'))
        proc.communicate()        
        if not check_files(HLA_hisat2_graph_index_fnames):
            print >> sys.stderr, "Error: indexing HLA failed!"
            sys.exit(1)

    # Build Bowtie2 indexes based on the above information
    HLA_bowtie2_index_fnames = ["hla.%d.bt2" % (i+1) for i in range(4)]
    HLA_bowtie2_index_fnames += ["hla.rev.%d.bt2" % (i+1) for i in range(2)]
    if not check_files(HLA_bowtie2_index_fnames):
        build_cmd = ["bowtie2-build",
                     "hla_sequences.fa",
                     "hla"]
        proc = subprocess.Popen(build_cmd, stdout=open("/dev/null", 'w'))
        proc.communicate()        
        if not check_files(HLA_bowtie2_index_fnames):
            print >> sys.stderr, "Error: indexing HLA failed!"
            sys.exit(1)

    # Read HLA alleles (names and sequences)
    refHLAs = {}
    for line in open("hla_backbone.fa"):
        if line.startswith('>'):
            HLA_name = line.strip().split()[0][1:]
            HLA_gene = HLA_name.split('*')[0]
            assert not HLA_gene in refHLAs
            refHLAs[HLA_gene] = HLA_name
    HLAs = {}
    for line in open("hla_sequences.fa"):
        if line.startswith(">"):
            HLA_name = line.strip().split()[0][1:]
            HLA_gene = HLA_name.split('*')[0]
            if not HLA_gene in HLAs:
                HLAs[HLA_gene] = {}
            if not HLA_name in HLAs[HLA_gene]:
                HLAs[HLA_gene][HLA_name] = ""
        else:
            HLAs[HLA_gene][HLA_name] += line.strip()
    HLA_names = {}
    for HLA_gene, data in HLAs.items():
        HLA_names[HLA_gene] = list(data.keys())

    # Read HLA variants, and link information
    Vars, Var_list = {}, {}
    for line in open("hla.snp"):
        var_id, var_type, allele, pos, data = line.strip().split('\t')
        gene = allele.split('*')[0]
        if not gene in Vars:
            Vars[gene] = {}
            assert not gene in Var_list
            Var_list[gene] = []
            
        assert not var_id in Vars[gene]
        Vars[gene][var_id] = [var_type, int(pos), data]
        Var_list[gene].append([int(pos), var_id])
    for gene, in_var_list in Var_list.items():
        Var_list[gene] = sorted(in_var_list)
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
    aligners = [
        ["hisat2", "graph"],
        # ["hisat2", "linear"],
        # ["bowtie2", "linear"]
        ]
    basic_test, random_test = False, True
    test_passed = {}
    test_list = []
    if basic_test:
        for HLA_gene, HLA_gene_alleles in HLA_names.items():
            for HLA_name in HLA_gene_alleles:
                test_list.append([HLA_name])
    if random_test:
        test_size = 500
        allele_count = 2
        for test_i in range(test_size):
            genes = HLA_names.keys()
            random.shuffle(genes)
            gene = genes[0]

            # daehwan - for debugging purposes
            # gene = "B"
            
            HLA_gene_alleles = HLA_names[gene]
            nums = [i for i in range(len(HLA_gene_alleles))]
            random.shuffle(nums)
            test_HLA_names = [HLA_gene_alleles[nums[i]] for i in range(allele_count)]
            test_HLA_names = sorted(test_HLA_names)
            test_list.append(test_HLA_names)
    for test_i in range(len(test_list)):
        # daehwan - for debugging purposes
        # if test_i + 1 not in [187, 195, 266, 346]:
        #    continue
        # two allele test (#266, #346)
        
        print >> sys.stderr, "Test %d" % (test_i + 1)
        test_HLA_names = test_list[test_i]

        # daehwan - for debugging purposes
        # test_HLA_names = ["A*31:14N"]
        
        for test_HLA_name in test_HLA_names:
            print >> sys.stderr, "\t%s" % (test_HLA_name)

        gene = test_HLA_names[0].split('*')[0]
        ref_allele = refHLAs[gene]
        ref_seq = HLAs[gene][ref_allele]
    
        # Simulate reads from two HLA alleles
        def simulate_reads(seq, read_len = 100):
            reads = []
            for i in range(0, len(seq) - read_len + 1):
                reads.append(seq[i:i+read_len])
            return reads
        HLA_reads = []
        for test_HLA_name in test_HLA_names:
            HLA_seq = HLAs[gene][test_HLA_name]
            HLA_reads += simulate_reads(HLA_seq)

        # Write reads into a fasta read file
        read_file = open('hla_input.fa', 'w')
        for read_i in range(len(HLA_reads)):
            print >> read_file, ">%d" % (read_i + 1)
            print >> read_file, HLA_reads[read_i]
        read_file.close()

        for aligner, index_type in aligners:
            print >> sys.stderr, "\n\t\t%s %s" % (aligner, index_type)
            # Align reads, and sort the alignments into a BAM file
            if aligner == "hisat2":
                hisat2 = os.path.join(ex_path, "hisat2")
                aligner_cmd = [hisat2]
                aligner_cmd += ["--no-unal"]
                if index_type == "linear":
                    aligner_cmd += ["-k", "10"]
                aligner_cmd += ["-x", "hla.%s" % index_type]
                # aligner_cmd += ["-x", "test"]
                aligner_cmd += ["-f", "hla_input.fa"]                
            elif aligner == "bowtie2":
                aligner_cmd = [aligner,
                               "--no-unal",
                               "-k", "10",
                               "-x", "hla",
                               "-f", "hla_input.fa"]
            else:
                assert False
            align_proc = subprocess.Popen(aligner_cmd,
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

            if index_type == "graph":
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
            else:
                os.system("mv hla_input_unsorted.bam hla_input.bam")

            # Read alignments
            alignview_cmd = ["samtools",
                             "view",
                             "hla_input.bam",
                             # "A*01:01:01:01:0-10000",
                             ]
            alignview_proc = subprocess.Popen(alignview_cmd,
                                              stdout=subprocess.PIPE,
                                              stderr=open("/dev/null", 'w'))

            # Count alleles
            HLA_counts, HLA_cmpt = {}, {}
            coverage = [0 for i in range(len(ref_seq) + 1)]
            num_reads, total_read_len = 0, 0
            if index_type == "graph":
                # Cigar regular expression
                cigar_re = re.compile('\d+\w')
                for line in alignview_proc.stdout:
                    cols = line[:-1].split()
                    read_id, flag, chr, pos, mapQ, cigar_str = cols[:6]
                    read_seq = cols[9]
                    if not chr.startswith(gene):
                        continue

                    num_reads += 1
                    total_read_len += len(read_seq)

                    flag = int(flag)
                    pos = int(pos)

                    # daehwan - for debugging purposes
                    debug = False
                    if pos - 1 == 1809 and False:
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
                            # Update coverage
                            assert right_pos + length < len(coverage)
                            coverage[right_pos] += 1
                            coverage[right_pos + length] -= 1

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
                    for HLA_name in HLA_names[gene]:
                        HLA_count_per_read[HLA_name] = 0

                    def add_count(var_id, add):
                        assert var_id in Links
                        alleles = Links[var_id]
                        for allele in alleles:
                            HLA_count_per_read[allele] += add

                    # Decide which allele(s) a read most likely came from
                    # also sanity check - read length, cigar string, and MD string
                    for var_id, data in Vars[gene].items():
                        var_type, var_pos, var_data = data
                        if var_type != "deletion":
                            continue
                        if left_pos >= var_pos and right_pos <= var_pos + int(var_data):
                            add_count(var_id, -1)                            
                    ref_pos, read_pos, cmp_cigar_str, cmp_MD = left_pos, 0, "", ""
                    cigar_match_len, MD_match_len = 0, 0            
                    for cmp in cmp_list:
                        type = cmp[0]
                        length = cmp[2]
                        if type == "match":
                            var_idx = lower_bound(Var_list[gene], ref_pos)
                            while var_idx < len(Var_list[gene]):
                                var_pos, var_id = Var_list[gene][var_idx]
                                if ref_pos + length <= var_pos:
                                    break
                                if ref_pos <= var_pos:
                                    var_type, _, var_data = Vars[gene][var_id]
                                    if var_type == "insertion":
                                        if ref_pos < var_pos and ref_pos + length > var_pos + len(var_data):
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
                            var_idx = lower_bound(Var_list[gene], ref_pos)
                            while var_idx < len(Var_list[gene]):
                                var_pos, var_id = Var_list[gene][var_idx]
                                if ref_pos < var_pos:
                                    break
                                if ref_pos == var_pos:
                                    var_type, _, var_data = Vars[gene][var_id]
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
                            var_idx = lower_bound(Var_list[gene], ref_pos)
                            # daehwan - for debugging purposes
                            if debug:
                                print left_pos, cigar_str, MD, vars
                                print ref_pos, ins_seq, Var_list[gene][var_idx], Vars[gene][Var_list[gene][var_idx][1]]
                                # sys.exit(1)
                            while var_idx < len(Var_list[gene]):
                                var_pos, var_id = Var_list[gene][var_idx]
                                if ref_pos < var_pos:
                                    break
                                if ref_pos == var_pos:
                                    var_type, _, var_data = Vars[gene][var_id]
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
                            var_idx = lower_bound(Var_list[gene], ref_pos)
                            while var_idx < len(Var_list[gene]):
                                var_pos, var_id = Var_list[gene][var_idx]
                                if ref_pos < var_pos:
                                    break
                                if ref_pos == var_pos:
                                    var_type, _, var_data = Vars[gene][var_id]
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
                        if max_count == None:
                            max_count = count
                        elif count > max_count:
                            max_count = count
                    cur_cmpt = set()
                    # daehwan - for debugging purposes
                    allele1_found, allele2_found = False, False
                    for allele, count in HLA_count_per_read.items():
                        if count < max_count:
                            continue
                        if allele == "A*24:02:01:01" and False:
                            allele1_found = True
                        elif allele == "A*24:11N" and False:
                            allele2_found = True
                        cur_cmpt.add(allele)                    
                        if not allele in HLA_counts:
                            HLA_counts[allele] = 1
                        else:
                            HLA_counts[allele] += 1
                    if allele1_found != allele2_found:
                        if allele1_found:
                            print ("A*24:02:01:01\tread_id %s]" % read_id), left_pos, cigar_str, MD, Zs
                        else:
                            print ("A*24:11N\tread_id %s]" % read_id), left_pos, cigar_str, MD, Zs
                        print read_seq

                    cur_cmpt = sorted(list(cur_cmpt))
                    cur_cmpt = '-'.join(cur_cmpt)
                    if not cur_cmpt in HLA_cmpt:
                        HLA_cmpt[cur_cmpt] = 1
                    else:
                        HLA_cmpt[cur_cmpt] += 1

                # Coverage
                assert num_reads > 0
                read_len = int(total_read_len / float(num_reads))
                coverage_sum = 0
                for i in range(len(coverage)):
                    if i > 0:
                        coverage[i] += coverage[i-1]
                    coverage_sum += coverage[i]
                coverage_avg = coverage_sum / float(len(coverage))
                assert len(ref_seq) < len(coverage)
                for i in range(len(ref_seq)):
                    coverage_threshold = 1.0 * coverage_avg
                    if i < read_len:
                        coverage_threshold *= ((i+1) / float(read_len))
                    elif i + read_len > len(ref_seq):
                        coverage_threshold *= ((len(ref_seq) - i) / float(read_len))
                    if coverage[i] >= coverage_threshold:
                        continue
                    pseudo_num_reads = (coverage_threshold - coverage[i]) / read_len
                    var_idx = lower_bound(Var_list[gene], i + 1)
                    if var_idx >= len(Var_list[gene]):
                        var_idx = len(Var_list[gene]) - 1
                    cur_cmpt = set()
                    while var_idx >= 0:
                        var_pos, var_id = Var_list[gene][var_idx]
                        var_type, _, var_data = Vars[gene][var_id]
                        if var_type == "deletion":
                            del_len = int(var_data)
                            if i < var_pos:
                                break
                            if i + read_len < var_pos + int(var_data):
                                assert var_id in Links
                                cur_cmpt = cur_cmpt.union(set(Links[var_id]))
                        var_idx -= 1
                    if cur_cmpt:
                        cur_cmpt = '-'.join(list(cur_cmpt))
                        if not cur_cmpt in HLA_cmpt:
                            HLA_cmpt[cur_cmpt] = 0
                        HLA_cmpt[cur_cmpt] += pseudo_num_reads
            else:
                assert index_type == "linear"
                def add_alleles(alleles):
                    if not allele in HLA_counts:
                        HLA_counts[allele] = 1
                    else:
                        HLA_counts[allele] += 1

                    cur_cmpt = sorted(list(alleles))
                    cur_cmpt = '-'.join(cur_cmpt)
                    if not cur_cmpt in HLA_cmpt:
                        HLA_cmpt[cur_cmpt] = 1
                    else:
                        HLA_cmpt[cur_cmpt] += 1

                prev_read_id, prev_AS = None, None
                alleles = set()
                for line in alignview_proc.stdout:
                    cols = line[:-1].split()
                    read_id, flag, allele = cols[:3]
                    flag = int(flag)
                    if flag & 0x4 != 0:
                        continue

                    if not allele.startswith(gene):
                        continue

                    AS = None
                    for i in range(11, len(cols)):
                        col = cols[i]
                        if col.startswith("AS"):
                            AS = int(col[5:])
                    assert AS != None
                    if read_id != prev_read_id:
                        if alleles:
                            if aligner == "hisat2" or \
                                    (aligner == "bowtie2" and len(alleles) < 10):
                                add_alleles(alleles)
                            alleles = set()
                        prev_AS = None
                    if prev_AS != None and AS < prev_AS:
                        continue
                    prev_read_id = read_id
                    prev_AS = AS
                    alleles.add(allele)

                if alleles:
                    add_alleles(alleles)
                    
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
                        print >> sys.stderr, "\t\t\t*** %d ranked %s (count: %d)" % (count_i + 1, test_HLA_name, count[1])
                        """
                        if count_i > 0 and HLA_counts[0][1] > count[1]:
                            print >> sys.stderr, "Warning: %s ranked first (count: %d)" % (HLA_counts[0][0], HLA_counts[0][1])
                            assert False
                        else:
                            test_passed += 1
                        """
                if count_i < 10 and False:
                    print >> sys.stderr, "\t\t\t\t%d %s (count: %d)" % (count_i + 1, count[0], count[1])

            def normalize(prob):
                total = sum(prob.values())
                for allele, mass in prob.items():
                    prob[allele] = mass / total

            def prob_diff(prob1, prob2):
                diff = 0.0
                for allele in prob1.keys():
                    if allele in prob2:
                        diff += abs(prob1[allele] - prob2[allele])
                    else:
                        diff += prob1[allele]
                return diff

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

            # daehwan - for debugging purposes
            if len(test_HLA_names) != 2:
                HLA_prob, HLA_prob_next = {}, {}
                for cmpt, count in HLA_cmpt.items():
                    alleles = cmpt.split('-')
                    for allele in alleles:
                        if allele not in HLA_prob:
                            HLA_prob[allele] = 0.0
                        HLA_prob[allele] += (float(count) / len(alleles))

                normalize(HLA_prob)
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
                for allele, prob in HLA_prob.items():
                    allele_len = len(HLAs[gene][allele])
                    HLA_prob[allele] /= float(allele_len)
                normalize(HLA_prob)
                HLA_prob = [[allele, prob] for allele, prob in HLA_prob.items()]
                
                HLA_prob = sorted(HLA_prob, cmp=HLA_prob_cmp)
                success = [False for i in range(len(test_HLA_names))]
                found_list = [False for i in range(len(test_HLA_names))]
                for prob_i in range(len(HLA_prob)):
                    prob = HLA_prob[prob_i]
                    found = False
                    for test_i in range(len(test_HLA_names)):
                        test_HLA_name = test_HLA_names[test_i]
                        if prob[0] == test_HLA_name:
                            print >> sys.stderr, "\t\t\t*** %d ranked %s (abundance: %.2f%%)" % (prob_i + 1, test_HLA_name, prob[1] * 100.0)
                            if prob_i < len(success):
                                success[prob_i] = True
                            found_list[test_i] = True
                            found = True                        
                    if not False in found_list:
                        break
                    if not found:
                        print >> sys.stderr, "\t\t\t\t%d ranked %s (abundance: %.2f%%)" % (prob_i + 1, prob[0], prob[1] * 100.0)

            else:
                assert len(test_HLA_names) == 2
                HLA_prob, HLA_prob_next = {}, {}
                for cmpt, count in HLA_cmpt.items():
                    alleles = cmpt.split('-')
                    for allele1 in alleles:
                        for allele2 in HLA_names[gene]:
                            if allele1 < allele2:
                                allele_pair = "%s-%s" % (allele1, allele2)
                            else:
                                allele_pair = "%s-%s" % (allele2, allele1)
                            if not allele_pair in HLA_prob:
                                HLA_prob[allele_pair] = 0.0
                            HLA_prob[allele_pair] += (float(count) / len(alleles))

                # Choose top allele pairs
                def choose_top_alleles(HLA_prob):
                    HLA_prob_list = [[allele_pair, prob] for allele_pair, prob in HLA_prob.items()]
                    HLA_prob_list = sorted(HLA_prob_list, cmp=HLA_prob_cmp)
                    HLA_prob = {}
                    best_prob = HLA_prob_list[0][1]
                    for i in range(len(HLA_prob_list)):
                        allele_pair, prob = HLA_prob_list[i]
                        if prob * 2 <= best_prob:
                            break                        
                        HLA_prob[allele_pair] = prob
                    normalize(HLA_prob)
                    return HLA_prob
                HLA_prob = choose_top_alleles(HLA_prob)

                def next_prob(HLA_cmpt, HLA_prob):
                    HLA_prob_next = {}
                    for cmpt, count in HLA_cmpt.items():
                        alleles = cmpt.split('-')
                        prob = 0.0
                        for allele in alleles:
                            for allele_pair in HLA_prob.keys():
                                if allele in allele_pair:
                                    prob += HLA_prob[allele_pair]
                        for allele in alleles:
                            for allele_pair in HLA_prob.keys():
                                if not allele in allele_pair:
                                    continue
                                if allele_pair not in HLA_prob_next:
                                    HLA_prob_next[allele_pair] = 0.0
                                HLA_prob_next[allele_pair] += (float(count) * HLA_prob[allele_pair] / prob)
                    normalize(HLA_prob_next)
                    return HLA_prob_next

                diff, iter = 1.0, 0
                while diff > 0.0001 and iter < 1000:
                    HLA_prob_next = next_prob(HLA_cmpt, HLA_prob)
                    diff = prob_diff(HLA_prob, HLA_prob_next)
                    HLA_prob = HLA_prob_next
                    HLA_prob = choose_top_alleles(HLA_prob)
                    iter += 1

                HLA_prob = [[allele_pair, prob] for allele_pair, prob in HLA_prob.items()]
                HLA_prob = sorted(HLA_prob, cmp=HLA_prob_cmp)

                success = [False]
                for prob_i in range(len(HLA_prob)):
                    allele_pair, prob = HLA_prob[prob_i]
                    allele1, allele2 = allele_pair.split('-')
                    if allele1 in test_HLA_names and allele2 in test_HLA_names:
                        print >> sys.stderr, "\t\t\t*** %d ranked %s (abundance: %.2f%%)" % (prob_i + 1, allele_pair, prob * 100.0)
                        if prob_i == 0:
                            success[0] = True
                        break
                    print >> sys.stderr, "\t\t\t\t%d ranked %s (abundance: %.2f%%)" % (prob_i + 1, allele_pair, prob * 100.0)
                        
            if not False in success:
                aligner_type = "%s %s" % (aligner, index_type)
                if not aligner_type in test_passed:
                    test_passed[aligner_type] = 1
                else:
                    test_passed[aligner_type] += 1


    for aligner_type, passed in test_passed.items():
        print >> sys.stderr, "%s:\t%d/%d passed (%.2f%%)" % (aligner_type, passed, len(test_list), passed * 100.0 / len(test_list))
    
        
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
    parser.add_argument("--sequence-type",
                        dest="sequence_type",
                        type=str,
                        default="gene",
                        help="Sequence type: (1) gene, (2) chromosome, and (3) genome.")
    parser.add_argument('-v', '--verbose',
        dest='verbose',
        action='store_true',
        help='also print some statistics to stderr')

    args = parser.parse_args()
    if not args.sequence_type in ["gene", "chromosome", "genome"]:
        print >> sys.stderr, "Error: --seqeuence-type (%s) must be one of gene, chromosome, and genome" % (args.sequence_type)
        sys.exit(1)
    random.seed(1)
    test_HLA_genotyping(args.base_fname, args.sequence_type, args.verbose)
